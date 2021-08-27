#' Produces the new info matrix from Simpati's classification results
#'
#' After Simpati has classified patients, this function can take Simpati's classification results
#' and produce the new info matrix related to the original input patients
#'
#' @param cons_res Output of the get_perfs function
#' @param vars_list Variable list of the Simpati's run that has produced the cons_res oject
#' @return Return new info matrix of the original patients
#' @export
#'
get_new_info = function(cons_res,vars_list){
  #Produce the new clinical information matrix -----
  #Get the testing patients and their predicted labels
  pred_te_df=data.frame(te_name=cons_res$testing_prof,
                        pred=cons_res$predicted_classes, stringsAsFactors = F)
  #Map the labels to their orignal status
  pred_te_df$class=ifelse(pred_te_df$pred==vars_list$Aletter,
                          vars_list$tab_status[1,1],
                          vars_list$tab_status[2,1])
  #Update the original clinical matrix with the new data
  info=vars_list$info
  new_info=info[-match(pred_te_df$te_name,info$shortID),]
  old_te_info=info[match(pred_te_df$te_name,info$shortID),]
  old_te_info[,2]=pred_te_df$class
  new_info=rbind(new_info[,c(1,2)],old_te_info[,c(1,2)])
  new_info=new_info[order(new_info[,2]),]
  return(new_info)
}

#' Extracts the pathways detected as signature during a Simpati's run of patient classification
#'
#' After Simpati has classified patients, this function can take Simpati's classification results
#' and provide the pathways detected as signatures
#'
#' @param PSNs_info_te Output of the predict_pats function
#' @param vars_list Variable list of the Simpati's run that has produced the cons_res oject
#' @return Return the feature list of the signature pathways used in the classification
#' @export
#'
get_used_pathways = function(PSNs_info_te,vars_list){
  #Get the name of the pathways which have been most frequently used to predict
  pathways_freq=as.data.frame(table(PSNs_info_te$pathway_name), stringsAsFactors = F)
  pathways_used=pathways_freq[pathways_freq[,2]>=quantile(pathways_freq[,2],probs = 0.4),1]
  #Get the original set of pathways
  pathways_l=vars_list$pathways_l
  pathways_used_l=pathways_l[pathways_used]
  return(pathways_used_l)
}

#' Detects signature pathways related to two patient classes
#'
#' Detects the pathways that are signature between two patient classes.
#' It does not depend by a step of classification. It computes a patient similarity network
#' per pathway and finds the PSNs which significantly separate the two classes.
#' It tests if a signature pathway specific PSN is signficant based on a permutation approach.
#'
#' @param geno_p Matrix of patient's profiles
#' @param pathways_l1 List of feature's sets
#' @param vars_list Variable list of the Simpati's run created with the set_variables function
#' @param th_pred Default 3, Power threshold to understand if Simpati should try to perform the Best Friend Connector algorithm
#' @param n_test Default 200, number of permutations to determine the p.value of significancy
#' @param Pv_th Default 0.05, value of threshold to assess a pathway-specific PSN as significant
#' @param n_cores Default 2, value which sets the number of cores to use for the parallel computation
#' @return Return PSN_info_df: matrix of information related to the tested pathway specific PSNs
#' PSN_comp_l: list of vectorized pathway specific PSNs
#' outlier_df: matrix which summerizes the likelihood for each patient to be outlier of its a priori own class
#' @import parallel
#' @import doParallel
#' @import matrixStats
#' @import stats
#' @import foreach
#' @export
#'
get_sign_pathways = function(geno_p,pathways_l1,vars_list,th_pred=3,n_test=200,Pv_th=0.05,n_cores=2){
  # library("parallel")
  # library("doParallel")
  # library("matrixStats")

  #matrix-wise 0-1 standardization for downstream analysis
  #the propagation gene value is standardized to reflect how much is high/low with respect
  #all the other genes in the matrix both row and column wise
  geno_p=scaling_matrix(geno_p)

  #Ranked version for the result phase
  #1-100 column wise stand. to get the importance of a gene in patient
  geno_pR=column_rank(geno_p,par=T)

  #Split the pathway list into chunks ----
  n_pathways=length(pathways_l1)
  pathways_indxs=seq(1,n_pathways)
  max <- length(pathways_l1)/n_cores
  x <- seq_along(pathways_indxs)
  pathway_sets <- split(pathways_indxs, ceiling(x/max))

  #Create the sample permuations for test ----
  #Generate permutations of genetic profiles for testing a signature pathway
  a=seq(1,ncol(geno_p))
  perms <- t(as.data.frame(replicate(n_test, sample(a))))
  perms = rbind(a,perms)

  #Backup of the original matrix for directional similarity ----
  dirs=c("up-inv","down-inv")
  orig_geno_p=geno_p

  if(.Platform$OS.type == "unix") {
    cat("Parallel for linux \n")
    cl <- makeCluster(n_cores,type="FORK");
  } else {
    cat("Parallel for windows \n")
    cl <- makeCluster(n_cores);
  }
  registerDoParallel(cl);data_list=list();

  PSNs_l=list();
  vars_list=vars_list[c("Aletter","Bletter","min_th_cl","max_th_cl")]
  PSNs_l=foreach(k_chunk=1:length(pathway_sets),.inorder=FALSE,.noexport=c(),
                 .packages = c("matrixStats"),.export = c("get_genes_desc",
                                                          "Adj_Weighted_Jaccard",
                                                          "Adj_Weighted_Jaccard_gene",
                                                          "autri",
                                                          "get_power",
                                                          "vectorizing_matrix",
                                                          "strong.comm.BB.cor",
                                                          "strong.comm.AA.cor",
                                                          "corr_BB_mat",
                                                          "corr_AA_mat",
                                                          "area_quadr",
                                                          "twoDist",
                                                          "strong.comm",
                                                          "th_pred",
                                                          "get_A_prob",
                                                          "get_B_prob",
                                                          "get_strength_sim",
                                                          "sim_WJ")) %dopar%
    {

      #Get the indexes of the chunck of pathways that the core has to process
      pathway_set=pathway_sets[[k_chunk]]
      PSNs_chunk_l=list()

      for(k_path in pathway_set){
        for(dir in dirs){
          if(dir=="up-inv"){geno_p=orig_geno_p}else{geno_p=1-geno_p}
          cat(k_path,"on", length(pathway_set), "\n")
          #Get the indexes of the genes belonging to the pathway
          pathway_indexes=pathways_l1[[k_path]]
          pathway_name=names(pathways_l1)[k_path]
          pathway_name=paste(pathway_name,dir,sep=" ")
          #Compute gene profiling
          ps_geno_p=geno_p[pathway_indexes,] #with only the training patients for the analysis
          gs_p=get_genes_desc(ps_geno_p,vars_list$Aletter,vars_list$Bletter)

          #Find the subset of genes in the pathway that maximize the PSN power due to two factors -----
          #how much separate the classes based on propagation value
          g_ranks=unique(gs_p$rank_dist)
          g_ranks=g_ranks[order(g_ranks,decreasing = F)]
          max_lv=0;max_g_dist_rank=0;best_m_sim_tr=NA;best_gs=NA;best_direction=c("FAIR","0")
          for(k_gp_test in 1:length(g_ranks)){
            g_rank=g_ranks[k_gp_test]
            gs_test=gs_p$rank_dist>=g_rank

            if(sum(gs_test)>=3){
              m_sim_tr=Adj_Weighted_Jaccard(ps_geno_p[gs_test,])
              direction=get_power(m_sim_tr,vars_list$Aletter,vars_list$Bletter)
              direction_lv=as.numeric(direction[2]);direction_cl=direction[1];
              if(direction_lv>max_lv){
                max_lv=direction_lv
                max_g_dist_rank=g_rank
                best_m_sim_tr=m_sim_tr
                best_gs=gs_test
                best_direction=direction
              }
            }
          }

          #how much separate the classes based on WJ value
          g_sds=quantile(gs_p$diff_sd,probs = c(0.25,0.50,0.75))
          max_g_sd=0;
          for(k_gp_test in 1:length(g_sds)){
            g_sd=g_sds[k_gp_test]
            gs_test=gs_p$diff_sd>=g_sd & gs_p$rank_dist>=max_g_dist_rank

            if(sum(gs_test)>=3){
              m_sim_tr=Adj_Weighted_Jaccard(ps_geno_p[gs_test,])
              direction=get_power(m_sim_tr,vars_list$Aletter,vars_list$Bletter)
              direction_lv=as.numeric(direction[2]);direction_cl=direction[1];
              if(direction_lv>max_lv){
                max_lv=direction_lv
                max_g_sd=g_sd
                best_m_sim_tr=m_sim_tr
                best_gs=gs_test
                best_direction=direction
              }
            }
          }

          #Clean
          rm(direction,direction_lv,direction_cl,m_sim_tr,gs_test,g_sd,g_sds,g_ranks,g_rank)

          #Prediction step -----
          if(max_lv>0){

            for(t in 1:nrow(perms)){
              #Get the submatrix of the pathway from geno matrix
              m_profiles=geno_p[pathway_indexes[best_gs],perms[t,]] #with all the patients for the class prediction
              m_profilesR=geno_pR[pathway_indexes[best_gs],perms[t,]] #with all the ranked patient data for the result processing
              colnames(m_profiles)=colnames(geno_p);colnames(m_profilesR)=colnames(geno_pR)
              #Compute the related PSNs
              m_sim=Adj_Weighted_Jaccard(m_profiles)

              if(t==1){
                m_simXgene=Adj_Weighted_Jaccard_gene(m_profiles)
                #Decompose and compact each matrix and PSN
                m_prof_l=vectorizing_matrix(m_profiles);
                m_profR_l=vectorizing_matrix(m_profilesR);
                m_sim_l=vectorizing_matrix(m_sim);
              }

              #Clean
              rm(m_profiles,m_profilesR)

              #if AA edges are stronger than BB
              if(best_direction[1]==vars_list$Aletter){
                nBs_in=length(grep(vars_list$Bletter,colnames(m_sim)))
                if(nBs_in>=4){
                  mB_weak_cor=strong.comm.BB.cor(mat=m_sim,nAs=length(grep(vars_list$Aletter,colnames(m_sim))),thW=0.4,thGroup=vars_list$min_th_cl,strong = F)
                }else{
                  mB_weak_cor=m_sim
                }
                m_sim_cl=strong.comm.AA.cor(mB_weak_cor,nAs=length(grep(vars_list$Aletter,colnames(mB_weak_cor))),thW=0.8,thGroup=vars_list$max_th_cl,strong = T)
                # m_sim_cl=strong.comm.AA.cor(mB_weak_cor,nAs=length(grep(vars_list$Aletter,colnames(mB_weak_cor))),thW=vars_list$max_th_cl,thGroup=vars_list$max_th_cl,strong = T)
                rm(mB_weak_cor)

                #if BB edges are stronger than AA
              }else{
                nAs_in=length(grep(vars_list$Aletter,colnames(m_sim)))
                if(nAs_in>=4){
                  mA_weak_cor=strong.comm.AA.cor(mat=m_sim,nAs=length(grep(vars_list$Aletter,colnames(m_sim))),thW=0.4,thGroup=vars_list$min_th_cl,strong = F)
                }else{
                  mA_weak_cor=m_sim
                }
                m_sim_cl=strong.comm.BB.cor(mat=mA_weak_cor,nAs=length(grep(vars_list$Aletter,colnames(mA_weak_cor))),thW=0.8,thGroup=vars_list$max_th_cl,strong = T)
                #m_sim_cl=strong.comm.BB.cor(mat=mA_weak_cor,nAs=length(grep(vars_list$Aletter,colnames(mA_weak_cor))),thW=vars_list$max_th_cl,thGroup=vars_list$max_th_cl,strong = T)
                rm(mA_weak_cor)
              }

              power_info=get_power(m_sim_cl,vars_list$Aletter,vars_list$Bletter)
              psn_lv=as.numeric(power_info[2]);psn_cl=power_info[1]
              if(t==1){
                psn_lvs=psn_lv
                power_info1=power_info
                m_sim_cl1=m_sim_cl
              }else{
                psn_lvs=c(psn_lvs,psn_lv)
              }
            }

            orig_psn_lv=psn_lvs[1]
            perm_psn_lvs=psn_lvs[-1]
            pv=sum(perm_psn_lvs>=orig_psn_lv)/length(perm_psn_lvs)
            outliers=setdiff(colnames(geno_p),colnames(m_sim_cl1))
          }else{
            power_info1=c("NA",0)
            outliers="NA"
            pv=999
            m_simXgene=m_prof_l=m_profR_l=m_sim_l=NULL
          }

          PSN_info=data.frame(pathway_name=pathway_name,
                              sign_class=power_info1[1],power=power_info1[2],
                              direction=dir,
                              Pvalue=pv,stringsAsFactors = F)
          PSNs_chunk_l[[pathway_name]]=list(PSN_info=PSN_info,
                                            m_simXgene=m_simXgene,
                                            m_prof_l=m_prof_l,
                                            m_profR_l=m_profR_l,
                                            m_sim_l=m_sim_l,
                                            outliers=outliers)

          #Reset
          power_info1=c("NA",0)
          outliers="NA"
          pv=999
          m_simXgene=m_prof_l=m_profR_l=m_sim_l=NULL
        }
      }
      return(PSNs_chunk_l)
    }


  cat("Preparing the final results \n")
  #Merge the results of each chunk
  PSNs_l=do.call(c, PSNs_l)

  #Format and clean the matrix with the PSN data
  PSN_info_df = as.data.frame(rbindlist(sapply(PSNs_l, "[",1),use.names=FALSE),stringsAsFactors = F);
  PSN_info_df$power=as.numeric(PSN_info_df$power)
  PSN_comp_l = lapply(PSNs_l, "[",-1);

  #Remove data about not significant pathways
  PSN_info_df=PSN_info_df[PSN_info_df$Pvalue<=Pv_th & PSN_info_df$power>=th_pred,]
  PSN_comp_l=PSN_comp_l[PSN_info_df$pathway_name]

  #Extract outliers and compute their frequency in significant pathwyas
  outl_df=as.data.frame(table(unlist(sapply(PSN_comp_l, "[",5))),stringsAsFactors = F)
  outl_df$Freq=(outl_df$Freq*100)/(length(PSN_comp_l))
  outl_df=outl_df[order(outl_df$Freq,decreasing = T),]
  PSN_comp_l = lapply(PSN_comp_l, "[",-5);

  #Prepare resulting object
  PSN_data_l=list(PSN_info_df=PSN_info_df,PSN_comp_l=PSN_comp_l,outlier_df=outl_df)

  #Stop
  stopCluster(cl)
  return(PSN_data_l)
}

#' Finds "disase type" --> "gene" associations with disgnet database
#'
#' Detects "disase type" --> "gene" associations with disgnet database and return the ID of the publications
#' confirming the associations for each gene with the user-provided key type
#'
#' @param genes Vector of gene names that are possible to query in disgnet database
#' @param key_type_disgnet character string defining the disease type
#' @param db_disgnet disgnet database
#' @return Return db_lit: list of matrices, each one reports the publications found supporting the association
#' between a gene and the disease type defined with by the key word
#' count: percentage of records supporting the gene-disease type association
#' @import stats
#' @export
#'
get_SemType_associations=function(genes,key_type_disgnet,db_disgnet){
  cat("Key semantic types that you can type in:\n")
  SemTypes=unique(db_disgnet$diseaseSemanticType)[1:2]
  print(SemTypes);cat("\n")

  db_list=list()
  count=0
  for(g in genes){
    #Filter by gene
    indxs=which(db_disgnet$geneSymbol==g)
    if(length(indxs)==0){db_list[[g]]=list(count=count);next;}
    db1=db_disgnet[indxs,]
    #Filter by type
    indxs=which(db1$diseaseSemanticType==key_type_disgnet)
    if(length(indxs)==0){db_list[[g]]=list(count=count);next;}
    count=length(indxs)
    db_lit=db1[indxs,]
    #Aggregate
    #db_lit=stats::aggregate(pmid ~ geneSymbol+diseaseName+diseaseSemanticType, db_lit, c)

    db_lit=tryCatch(
      expr = {
        stats::aggregate(pmid ~ geneSymbol+diseaseName+diseaseSemanticType, db_lit, c)
      },
      error = function(e){
        db_lit
      }
    )

    #Add it to the result
    db_list[[g]]=list(db_lit=db_lit,count=count)
    count=0
  }
  return(db_list)
}

#' Finds "disase" --> "gene" associations with disgnet database
#'
#' Detects "disase" --> "gene" associations with disgnet database and return the ID of the publications
#' confirming the associations for each gene with the user-provided key type
#'
#' @param genes Vector of gene names that are possible to query in disgnet database
#' @param key_type_disgnet character vector defining the disease of interest
#' @param db_disgnet disgnet database
#' @return Return db_lit: list of matrices, each one reports the publications found supporting the association
#' between a gene and the disease defined with by the key words
#' count: percentage of disease key words associated to the genes of interest
#' @import stats
#' @export
#'
get_KeysDis_associations=function(genes,key_names,db_disgnet){
  count=0
  db_listXgene=list()
  db_litXkey=list()
  for(g in genes){
    #Filter by gene
    indxs=which(db_disgnet$geneSymbol == g)
    if(length(indxs)==0){
      db_listXgene[[g]]=list(db_litXkey,count=0)
      next;
    }
    db1=db_disgnet[indxs,]
    #Filter by key words in cancer table
    for(key_name in key_names){
      #Filter by type
      indxs=grep(key_name,db1$diseaseName,ignore.case = T)
      if(length(indxs)!=0){
        count=count+1
        db_litXkey1=db1[indxs,]
        #db_litXkey1=stats::aggregate(pmid ~ geneSymbol+diseaseName+diseaseSemanticType, db_litXkey1, c)
        db_litXkey1=tryCatch(
          expr = {
            db_litXkey1=stats::aggregate(pmid ~ geneSymbol+diseaseName+diseaseSemanticType, db_litXkey1, c)
          },
          error = function(e){
            db_litXkey1
          }
        )
        db_litXkey[[key_name]]=db_litXkey1
      }
    }

    if(count!=0){
      #Add it to the result
      db_listXgene[[g]]=list(db_litXkey=db_litXkey,
                             count=(count/length(key_names))*100)
    }else{
      db_listXgene[[g]]=list(db_litXkey,count=0)
    }

    db_litXkey=list()
    count=0
  }
  return(db_listXgene)
}

#' Finds "cancer" --> "gene" associations with human atlas database
#'
#' Detects "cancer" --> "gene" associations with human atlas database and return the ID of the publications
#' confirming the associations for each gene with the user-provided key type
#'
#' @param genes Vector of gene names that are possible to query in disgnet database
#' @param key_names character vector defining the tumor of interest
#' @param db_cancer_ATLAS human atlas database
#' @return Return matrix with the genes information about the cancer in which are involved
#' @import stats
#' @export
#'
get_cancerAtlas_associations=function(genes,key_names,db_cancer_ATLAS){
  df=data.frame(gene="NA",cancer_name="NA",gene_class="NA",stringsAsFactors = F)
  for(g in genes){
    indxs=which(db_cancer_ATLAS$`Gene name`==g)
    if(length(indxs)==0){next;}
    db1=db_cancer_ATLAS[indxs,]
    for(k in 1:length(key_names)){
      key_name=key_names[k]
      indxs=grep(key_name,db1$Cancer,ignore.case = T)
      if(length(indxs)==0){next;}
      db2=db1[indxs,]
      for(row_k in 1:nrow(db2)){
        if(sum(!is.na(db2[row_k,c(3,4,5,6)]))==0){next;}
        gene_class=names(which.max(db2[row_k,c(3,4,5,6)]))
        cancer_name=db2$Cancer[row_k]
        df1=data.frame(gene=g,cancer_name=cancer_name,
                       gene_class=gene_class,stringsAsFactors = F)
        df=rbind(df,df1)
      }
    }
  }
  df=df[-1,]
  return(df)
}

#' Finds tissues associations to the genes of interest
#'
#' Detects "tissue" --> "gene" associations with human atlas database and return the ID of the publications
#' confirming the associations for each gene with the user-provided key type
#'
#' @param genes Vector of gene names that are possible to query in disgnet database
#' @param key_names character vector defining the tumor of interest
#' @param db_cancer_ATLAS human atlas database
#' @return Return list for each gene with the user-defined tissues to which is associated
#' @import stats
#' @export
#'
get_normalAtlas_associations=function(genes,key_names,db_normal_ATLAS){
  tissueXgene=list()
  for(g in genes){
    count=0
    keys_found="NA"
    indxs=which(db_normal_ATLAS$`Gene name`==g)
    if(length(indxs)==0){
      tissueXgene[[g]]=list(tissues_found=keys_found[-1],count=(count/length(key_names))*100)
      next
    }
    db1=db_normal_ATLAS[indxs,]
    for(k in 1:length(key_names)){
      key_name=key_names[k]
      indxs=grep(key_name,db1$Tissue,ignore.case = T)
      if(length(indxs)==0){next;}
      keys_found=c(keys_found,key_name)
      count=count+1
    }
    keys_found=keys_found[-1]
    tissueXgene[[g]]=list(tissues_found=keys_found,count=(count/length(key_names))*100)
  }
  return(tissueXgene)
}


#' Wrapper: Finds diseases, tumours and tissues associations to the genes of interest
#'
#' Detects "disase type" --> "gene" associations, Detects "disase" --> "gene" associations,
#' Detects "cancer" --> "gene" associations and Detects "tissue" --> "gene" associations
#'
#' Detects "disase type" --> "gene" associations with disgnet database and return the ID of the publications
#' confirming the associations for each gene with the user-provided key type
#'
#' Detects "disase" --> "gene" associations with disgnet database and return the ID of the publications
#' confirming the associations for each gene with the user-provided key type
#'
#' Detects "cancer" --> "gene" associations with human atlas database and return the ID of the publications
#' confirming the associations for each gene with the user-provided key type
#'
#' Detects "tissue" --> "gene" associations with human atlas database and return the ID of the publications
#' confirming the associations for each gene with the user-provided key type
#'
#' @param gene_sets_l list of genes
#' @param key_words character vector defining the tumor of interest
#' @param disease_type character vector defining the type of disease to query in disgnet
#' @param db_disgnet disgnet database
#' @param db_normal_ATLAS human atlas database
#' @param db_cancer_ATLAS human atlas database
#' @param n_cores Default 2, number of cores to use in the parallel execution of the function
#' @param records_only Default TRUE, returns only the number of associations and not the information
#' @return Return matrix with the genes information about the tissue in which are involved
#' @import parallel
#' @import doParallel
#' @import matrixStats
#' @import stats
#' @importFrom foreach %dopar%
#' @export
#'
get_pathway_enrich = function(gene_sets_l,key_words,disease_type,
                              db_disgnet,db_cancer_ATLAS,db_normal_ATLAS,
                              n_cores=2, records_only=TRUE){
  keep_genes=unique(unlist(gene_sets_l))
  db_disgnet=db_disgnet[db_disgnet$geneSymbol %in% keep_genes,]
  db_cancer_ATLAS=db_cancer_ATLAS[db_cancer_ATLAS$`Gene name` %in% keep_genes,]
  db_normal_ATLAS=db_normal_ATLAS[db_normal_ATLAS$`Gene name` %in% keep_genes,]

  n_pathways=length(gene_sets_l)
  pathways_indxs=seq(1,n_pathways)
  max <- length(gene_sets_l)/n_cores
  x <- seq_along(pathways_indxs)
  pathway_sets <- split(pathways_indxs, ceiling(x/max))

  if(grep("Neoplastic",disease_type)>=1){
    cancer=TRUE
  }else{
    cancer=FALSE
  }


  if(.Platform$OS.type == "unix") {
    cat("Parallel for linux \n")
    cl <- makeCluster(n_cores,type="FORK");
  } else {
    cat("Parallel for windows \n")
    cl <- makeCluster(n_cores);
  }
  registerDoParallel(cl);data_list=list();

  enrXpathway_l=list();
  enrXpathway_l=foreach(k_chunk=1:length(pathway_sets),.inorder=F,.noexport=c(),
                        .packages = c("matrixStats"),
                        .export = c("get_SemType_associations",
                                    "get_KeysDis_associations",
                                    "get_cancerAtlas_associations",
                                    "get_normalAtlas_associations"))%dopar% {

                                      #Get the indexes of the chunck of pathways that the core has to process
                                      pathway_set=pathway_sets[[k_chunk]]
                                      PSNs_chunk_l=list()

                                      for(kp in pathway_set){
                                        p_name=names(gene_sets_l)[kp]
                                        genesXp=gene_sets_l[[kp]]

                                        db_SemType=get_SemType_associations(genes = genesXp,
                                                                            key_type_disgnet = disease_type,
                                                                            db_disgnet = db_disgnet)
                                        records_db_SemType=max(sapply(db_SemType,function(x){x$count}))

                                        db_KeysDis=get_KeysDis_associations(genes = genesXp,
                                                                            key_names = key_words,
                                                                            db_disgnet = db_disgnet)
                                        records_db_KeysDis=max(sapply(db_KeysDis,function(x){x$count}))

                                        if(cancer==FALSE){
                                          db_normal_association=get_normalAtlas_associations(genes = genesXp,
                                                                                             key_names = key_words,
                                                                                             db_normal_ATLAS = db_normal_ATLAS)
                                          records_db_normal=max(sapply(db_normal_association,function(x){x$count}))
                                          db_cancer_association=NULL
                                          records_db_cancer=0
                                        }else{
                                          db_cancer_association=get_cancerAtlas_associations(genes = genesXp,
                                                                                             key_names = key_words,
                                                                                             db_cancer_ATLAS=db_cancer_ATLAS)
                                          records_db_cancer=(nrow(db_cancer_association)*100)/length(genesXp)
                                          db_normal_association=NULL
                                          records_db_normal=0
                                        }

                                        if(records_only){
                                          PSNs_chunk_l[[p_name]][["info"]]=data.frame(records_db_SemType=records_db_SemType,
                                                                                      records_db_KeysDis=records_db_KeysDis,
                                                                                      records_db_cancer=records_db_cancer,
                                                                                      records_db_normal=records_db_normal)
                                        }else{
                                          PSNs_chunk_l[[p_name]][["info"]]=data.frame(records_db_SemType=records_db_SemType,
                                                                                      records_db_KeysDis=records_db_KeysDis,
                                                                                      records_db_cancer=records_db_cancer,
                                                                                      records_db_normal=records_db_normal)

                                          PSNs_chunk_l[[p_name]][["DBs_results"]]=list(db_SemType=db_SemType,
                                                                                       db_KeysDis=db_KeysDis,
                                                                                       db_cancer_association=db_cancer_association,
                                                                                       db_normal_association=db_normal_association)
                                        }
                                      }

                                      return(PSNs_chunk_l)
                                    }

  enrXpathway_l=do.call(c, enrXpathway_l)

  enrXpathway_df = as.data.frame(rbindlist(sapply(enrXpathway_l, "[",1)),
                                 stringsAsFactors = F);
  enrXpathway_df$pathway_name=names(enrXpathway_l)
  enrXpathway_df=enrXpathway_df[,c(5,1,2,3,4)]

  dbXpathway_l=lapply(enrXpathway_l, "[",-1);

  stopCluster(cl)
  if(records_only){
    res=list(enrXpathway_df=enrXpathway_df)
  }else{
    res=list(enrXpathway_df=enrXpathway_df,dbXpathway_l=dbXpathway_l)
  }
  return(res)
}
