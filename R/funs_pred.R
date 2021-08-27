#' Perform the classification task to predict testing patient classes
#'
#' For each pathway, it builds a patient similarity network, it determines if it is signature due to the
#' training patient similarities and in case, it predicts the class of the testing patient.
#'
#' @param vars_l Variable list of the Simpati's run created with the set_variables function
#' @return Produce a matrix s.t. for each pathway, there are information about the PSN power and the signature class
#' @import data.table
#' @import matrixStats
#' @import scales
#' @import parallel
#' @import doParallel
#' @import matrixStats
#' @import stats
#' @import foreach
#' @export
#'
predict_pats = function(vars_l){
  #Set variables
  geno_p=vars_l$geno_p
  pathways_l1=vars_l$pathways_indxs_l
  testing_pats=vars_l$testing_pats
  n_cores=vars_l$n_cores
  n_min_genes=min(sapply(pathways_l1,length))
  vars_l=vars_l[c("th_pred","Aletter","Bletter","min_th_cl","max_th_cl")]

  #matrix-wise 0-1 standardization for downstream analysis
  geno_p=scaling_matrix(geno_p,capping = TRUE)

  #Split the pathway list into chunks ----
  n_pathways=length(pathways_l1)
  pathways_indxs=seq(1,n_pathways)
  max <- length(pathways_l1)/n_cores
  x <- seq_along(pathways_indxs)
  pathway_sets <- split(pathways_indxs, ceiling(x/max))

  #One core one chunk of pathways to process ----
  if(.Platform$OS.type == "unix") {
    cat("Parallel for linux \n")
    cl <- makeCluster(n_cores,type="FORK");
  } else {
    cat("Parallel for windows \n")
    cl <- makeCluster(n_cores);
  }
  registerDoParallel(cl);

  #Pathway analysis -----

  PSNs_l=list();
  PSNs_l=foreach(k_chunk=1:length(pathway_sets),.inorder=TRUE,.noexport=c(),.verbose = F,.errorhandling="pass",
                 .packages = c("matrixStats","Simpati"),.export = c("get_genes_desc",
                                                          "Adj_Weighted_Jaccard",
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
                                                          "get_A_prob",
                                                          "get_B_prob",
                                                          "get_strength_sim",
                                                          "sim_WJ")) %dopar%
    {
      #Get the indexes of the chunck of pathways that the core has to process
      pathway_set=pathway_sets[[k_chunk]]
      PSNs_chunk_l=list()
      for(k_path in pathway_set){

        cat(k_path,"on", length(pathway_set), "\n")
        #Get the indexes of the genes belonging to the pathway
        pathway_indexes=pathways_l1[[k_path]]
        pathway_name=names(pathways_l1)[k_path]
        #Initiliaze the vector of predicted labels for testing patients
        pred_labels="start";pred_th_A=pred_th_B=999;A_strength=B_strength=0;
        #Initiliaze the list that will contain the info regarding this running PNS
        PSN_infoXtest_l=list()
        #Initiliaze the counter of predicted testing patients
        pass=0

        #For each testing patient
        for(k_tepat in 1:length(testing_pats)){
          #Get the testing patient
          tepat=testing_pats[k_tepat]
          #Get the training patients
          keep_trpats=colnames(geno_p)!=tepat

          ##Get the submatrix of the pathway from geno matrix
          ps_tr_geno_p=geno_p[pathway_indexes,keep_trpats] #with only the training patients for the analysis

          #Compute gene profiling
          genes_desc=get_genes_desc(mp=ps_tr_geno_p, Aletter=vars_l$Aletter, Bletter=vars_l$Bletter)

          #Find the subset of genes in the pathway that maximize the initial PSN power due to two factors -----
          #PSN power: grade of separability between the two classes
          #This operation is performed without the testing patient
          g_ranks=unique(genes_desc$rank_dist)
          g_ranks=g_ranks[order(g_ranks,decreasing = F)]
          max_lv=-1;max_g_dist_rank=0;best_m_sim_tr=NA;best_gs=NA;best_direction=c("FAIR","0")
          for(k_gp_test in 1:length(g_ranks)){
            g_rank=g_ranks[k_gp_test]
            gs_test=genes_desc$rank_dist>=g_rank

            if(sum(gs_test)>=n_min_genes){
              m_sim_tr=Adj_Weighted_Jaccard(ps_tr_geno_p[gs_test,])

              m_sim_cl=m_sim_tr
              Aletter=vars_l$Aletter
              Bletter=vars_l$Bletter

              power_info=get_power(m_sim_tr,vars_l$Aletter,vars_l$Bletter)
              power_lv=power_info$power
              if(power_lv>max_lv){
                max_lv=power_lv
                max_g_dist_rank=g_rank
                best_m_sim_tr=m_sim_tr
                best_gs=gs_test
                best_direction=power_info
              }
            }
          }
          cat("0 \n")
          #how much separate the classes based on difference in classes sd values
          g_sds=quantile(genes_desc[best_gs,]$diff_sd,probs = c(0.25,0.50,0.75))
          max_g_sd="NA";
          for(k_gp_test in 1:length(g_sds)){
            name_g_sd=names(g_sds)[k_gp_test]
            g_sd=g_sds[k_gp_test]
            gs_test=genes_desc$diff_sd>=g_sd & genes_desc$rank_dist>=max_g_dist_rank

            if(sum(gs_test)>=n_min_genes){
              m_sim_tr=Adj_Weighted_Jaccard(ps_tr_geno_p[gs_test,])
              power_info=get_power(m_sim_tr,vars_l$Aletter,vars_l$Bletter)
              power_lv=power_info$power
              if(power_lv>max_lv){
                max_lv=power_lv
                max_g_sd=name_g_sd
                best_m_sim_tr=m_sim_tr
                best_gs=gs_test
                best_direction=power_info
              }
            }
          }

          #Clean
          rm(power_lv,power_info,m_sim_tr,gs_test,g_sd,g_sds,g_ranks,g_rank,
             name_g_sd,genes_desc,ps_tr_geno_p,k_gp_test)

          #Prediction step -----
          #The training PSN showed enough starting power, I consider it for the downstream analysis
          power_info=data.frame(class="NA",power=0,minSIGN=0,maxWEAK1=0,maxWEAK2=0)
          tepat_info=data.frame(pathway_name=pathway_name,power_info,testing_pat=tepat,
                                Stop=100,SSe=0,SWeS=0,DIFFe=0,assumption="NA",
                                Wtop=100,WWe=0,SWeW=0,DIFFe=0,assumption="NA",
                                A_strength=0,B_strength,DIFF_strength=0,
                                prediction="NA")
          if(max_lv>0){
            #The next two operations include the testing profile to avoid of computing
            #separately its similarities when I will have to add it for the prediction

            #Get the submatrix of the pathway from geno matrix with only the
            m_profiles=geno_p[pathway_indexes[best_gs],]
            #Compute the related PSN
            m_sim=Adj_Weighted_Jaccard(m_profiles)
            #Clean
            rm(m_profiles)

            #If the PSN doesn't have a strong starting power I try to enhance it with the algorithm
            if(max_lv<3){
              #if AA edges are stronger than BB
              if(best_direction[1]==vars_l$Aletter){
                nBs_in=length(grep(vars_l$Bletter,colnames(best_m_sim_tr)))
                if(nBs_in>=4){
                  mB_weak_cor=strong.comm.BB.cor(mat=best_m_sim_tr,nAs=length(grep(vars_l$Aletter,colnames(best_m_sim_tr))),thW=0.4,thGroup=vars_l$min_th_cl,strong = F)
                }else{
                  mB_weak_cor=best_m_sim_tr
                }
                m_sim_cl=strong.comm.AA.cor(mB_weak_cor,nAs=length(grep(vars_l$Aletter,colnames(mB_weak_cor))),thW=0.8,thGroup=vars_l$max_th_cl,strong = T)
                rm(mB_weak_cor)

                #if BB edges are stronger than AA
              }else{
                nAs_in=length(grep(vars_l$Aletter,colnames(best_m_sim_tr)))
                if(nAs_in>=4){
                  mA_weak_cor=strong.comm.AA.cor(mat=best_m_sim_tr,nAs=length(grep(vars_l$Aletter,colnames(best_m_sim_tr))),thW=0.4,thGroup=vars_l$min_th_cl,strong = F)
                }else{
                  mA_weak_cor=best_m_sim_tr
                }
                m_sim_cl=strong.comm.BB.cor(mat=mA_weak_cor,nAs=length(grep(vars_l$Aletter,colnames(mA_weak_cor))),thW=0.8,thGroup=vars_l$max_th_cl,strong = T)
                rm(mA_weak_cor)
              }
              cat("1 \n")
              if(sum(duplicated(rownames(m_sim_cl)))){
                dupl = duplicated(rownames(m_sim_cl))
                m_sim_cl = m_sim_cl[!dupl,!dupl]
              }

            }else{
              #While if the PSN is already strong with all the training profiles, I don't enhance the power
              m_sim_cl=best_m_sim_tr

              if(sum(duplicated(rownames(m_sim_cl)))){
                dupl = duplicated(rownames(m_sim_cl))
                m_sim_cl = m_sim_cl[!dupl,!dupl]
              }

            }

            power_info=get_power(m_sim_cl,vars_l$Aletter,vars_l$Bletter)
            cat("2 \n")
            if(power_info$power>=vars_l$th_pred){
              cat("**pred \n")
              #Count predicitons
              Ainfo=Ainfor=get_A_prob(tepat,m_sim,best_m_sim_tr,vars_l,power_info)
              Binfo=Binfor=get_B_prob(tepat,m_sim,best_m_sim_tr,vars_l,power_info)
              if(Ainfo$assumption=="ASign"){
                colnames(Ainfo)=c("Stop","SSe","SWeS","assumption")
                colnames(Binfo)=c("Wtop","WWe","SWeW","assumption")
                tepat_info1=cbind(pathway_name=pathway_name,power_info,
                                  Ainfo,Binfo,DIFF_top=round(abs(Ainfo$Stop-Binfo$Wtop),4))

              }else{
                colnames(Binfo)=c("Stop","SSe","SWeS","assumption")
                colnames(Ainfo)=c("Wtop","WWe","SWeW","assumption")
                tepat_info1=cbind(pathway_name=pathway_name,power_info,
                                  Binfo,Ainfo,DIFF_top=round(abs(Binfo$Stop-Ainfo$Wtop),4))
              }
              cat("3 \n")
              #Get how much the patient is similar to his BFs of As and Bs
              sim_strg=get_strength_sim(tepat,m_sim,m_sim_cl,vars_l,top_th = 0.50)
              tepat_info1$A_strength=round(sim_strg$A_strength,4)
              tepat_info1$B_strength=round(sim_strg$B_strength,4)
              tepat_info1$DIFF_strength=abs(round(sim_strg$A_strength,4)-round(sim_strg$B_strength,4))
              cat("4 \n")
              if((tepat_info1$A_strength>tepat_info1$B_strength)&
                (Ainfor$Atop<Binfor$Btop)){
                tepat_info1$MIN_top=Ainfor$Atop
                tepat_info1$MAX_strength=tepat_info1$A_strength
                tepat_info1$testing_pat=tepat
                tepat_info1$prediction=vars_l$Aletter
                tepat_info=tepat_info1
                pass=c(pass,k_tepat)
              }

              if((tepat_info1$A_strength<tepat_info1$B_strength)&
                 (Ainfor$Atop>Binfor$Btop)){
                tepat_info1$MIN_top=Binfor$Btop
                tepat_info1$MAX_strength=tepat_info1$B_strength
                tepat_info1$testing_pat=tepat
                tepat_info1$prediction=vars_l$Bletter
                tepat_info=tepat_info1
                pass=c(pass,k_tepat)
              }

            }
          }

          psnXtest_df=tepat_info
          rownames(psnXtest_df)="info";
          PSN_infoXtest_l[[k_tepat]]=psnXtest_df
        }

        #If the PSN has classfied a testing patient, then keep it
        if(length(pass[-1])>=2){
          #Keep info about only the made predictions
          keep_preds=pass[-1]
          PSN_infoXtest_l=PSN_infoXtest_l[keep_preds]
          #Compose the PSN info dataframe
          PSN_info=rbindlist(PSN_infoXtest_l,use.names=T)
          #Save into the shared list
          PSNs_chunk_l[[pathway_name]]=PSN_info
        }

      }

      PSNs_chunk_l=PSNs_chunk_l[sapply(PSNs_chunk_l,length)!=0]
      PSNs_info=rbindlist(PSNs_chunk_l)
      return(PSNs_info)
    }


  PSNs_info=rbindlist(PSNs_l)
  class(PSNs_info)="data.frame"
  stopCluster(cl)
  return(PSNs_info)
}
