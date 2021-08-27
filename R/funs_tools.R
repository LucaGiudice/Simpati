#' Preview of a matrix like object
#'
#' Allows to quickly see the first rows and columns of a matrix like object
#'
#' @param m A matrix like object
#' @param r number of rows to show, default 5
#' @param c number of columns to show, default 5
#' @return Preview of a matrix like object
#' @export
#'
see = function(m,r=5,c=5){
  if(nrow(m)<5){
    r=nrow(m)
  }
  if(ncol(m)<5){
    c=ncol(m)
  }
  m[1:r,1:c]
}

#' Checks vectors elements
#'
#' This function checks if two vectors are exactly the same
#'
#' @param x A vector
#' @param y A vector
#' @return The result of the check. 0 not match. 1 match.
#' @export
#'
check_elems = function(x,y){
  #Check if the vectors have the same length
  if(length(x)==length(y)){

    #Check if all elements in y are in x
    indx_ord=match(y,x)
    if(sum(!is.na(indx_ord))==0){
      cat("ERROR:Second vector has something that doesn't appear in the first \n")
      return(0)
    }

    #Check if all elements in x are in y
    indx_ord=match(x,y)
    if(sum(!is.na(indx_ord))==0){
      cat("ERROR:First vector has something that doesn't appear in the second \n")
      return(0)
    }

    #Check if all elements are also in the same position
    corr_ord=seq(1,length(y))
    err=0
    for(k in corr_ord){
      check=indx_ord[k]==k
      if(is.na(check)){
        err=c(err,k)
      }else{
        if(!check){
          err=c(err,k)
        }
      }
    }
    err=err[-1]

    if(length(err)==0){cat(">>Patient's genetic profiles and clinical data match \n");return(1)}
    if(length(err)!=0){cat("ERROR:Patient's genetic profiles and clinical data not match \n");return(0)}
  }else{
    cat("ERROR:Patient's genetic profiles and clinical data have different lengths \n");return(0)
  }
}

#' Harmonizer between the info matrix and the matrix of patient's profiles
#'
#' Remove patient's profiles that do not have the clinical information
#' Remove patient's info that do not have a profile in the feature matrix
#'
#' @param geno A matrix of genetic profiles: genes x samples
#' @param info A matrix of info about profiles: samples and binary classes
#' @return A list with the prepared geno and info matrix
#' @export
#'
harm_ds = function(geno,info){
  cat("*Harmonizing the info dataframe with the profiles matrix \n")
  #keep only the profiles with a matching id and class in the info matrix
  corr_indxs = match(info[, 1], colnames(geno))
  keep = !is.na(corr_indxs)
  info = info[keep, ]
  corr_indxs = corr_indxs[keep]
  geno = geno[, corr_indxs]
  cat(">>Data harmonized \n")
  res = list(info = info, geno = geno)
  return(res)
}

#' Controls the user-provided input data
#'
#' Prepare the user-provided input data to Simpati workflow
#' In case this step does not produce errors, then Simpati workflow is ready to classify and enrich the data
#'
#' @param geno A numeric matrix of patient profiles: features (e.g. genes) x samples
#' @param info A character matrix of info about the patients: patient_ID | class (e.g. clinical information)
#' @return A list with the prepared geno and info matrix
#' @import FNN
#' @import stats
#' @export
#'
data_preparation = function(geno,info,net,pathways_l,n_cores=5){
  cat("*Data preparation\n")
  indxs_gs=match(colnames(net),rownames(geno))
  indxs_gs=length(indxs_gs[!is.na(indxs_gs)])
  if(indxs_gs>0){
    cat(">>There are",indxs_gs,"genes matching between network and profiles\n")
  }else{
    cat("ERROR:There aren't genes matching between network and profiles\n")
    return(0)
  }

  indxs_gs=match(unique(unlist(pathways_l)),rownames(geno))
  indxs_gs=length(indxs_gs[!is.na(indxs_gs)])
  if(indxs_gs>0){
    cat(">>There are",indxs_gs,"genes matching between pathways and profiles\n")
  }else{
    cat("ERROR:There aren't genes matching between pathways and profiles\n")
    return(0)
  }

  #Check that IDs are equal and ordered the same as colnames(geno)
  check=check_elems(colnames(geno),info[,1])
  #If colnames(geno) and ids in info matirx don't match, harmonize
  if(check==0){
    harm_ds_l=harm_ds(geno,info)
    info=harm_ds_l$info
    geno=harm_ds_l$geno
  }

  n_classes=length(unique(info[,2]))
  if(n_classes!=2){
    cat("ERROR:Profile classes in info matrix are not two\n")
    return(0)
  }

  #Compact
  kmeans_indxs=function(x,th=3){
    x1=t(x)
    kmeans.x <- kmeans(x1, 10)
    y <- get.knnx(x1, kmeans.x$centers, round(nrow(x1)/th))
    idx1 <- sort(y$nn.index[1, ])
    return(idx1)
  }

  if(ncol(geno)>=100){
    gr1=which(info[1,2]==info[,2])
    gr2=which(info[1,2]!=info[,2])
    geno_p=rwrProp(net, geno, n_cores=n_cores);

    if(ncol(geno)>=100 & ncol(geno)<200){
      i1=kmeans_indxs(geno_p[,gr1],3);i2=kmeans_indxs(geno_p[,gr2],3)+length(gr1)
    }

    if(ncol(geno)>=200 & ncol(geno)<300){
      i1=kmeans_indxs(geno_p[,gr1],5);i2=kmeans_indxs(geno_p[,gr2],5)+length(gr1)
    }

    if(ncol(geno)>=300){
      i1=kmeans_indxs(geno_p[,gr1],7);i2=kmeans_indxs(geno_p[,gr2],7)+length(gr1)
    }

    geno=geno[,c(i1,i2)]
    info=info[c(i1,i2),]
  }

  #Order sample data s.t. one class of sample comes before the other
  neword=order(info[,2])
  info=info[neword,]
  geno=geno[,neword]
  #Save results and return
  res=list(geno=geno,info=info,net=net,pathways_l=pathways_l)
  #Visualize
  cat("*Input data are well formatted\n")
  cat(">>Class names:",names(table(info[,2])));cat("\n")
  cat(">>Class sizes:",as.vector(table(info[,2])));cat("\n")
  return(res)
}

#' Creates and sets the list of variable required to run Simpati
#'
#' Sets the variables based on the input data and required by Simpati to run any other function of the workflow
#'
#' @param data_rdy_l Simpati object containing checked and prepared input data
#' @param dataset_name The name of the dataset or project
#' @param disease_type The DISGNET semantic type of the sample condition: "Neoplastic Process" or "Disease or Syndrome"
#' @param key_words A character vector of one or more words related to the sample condition (e.g. KIRC: Kidney, Carcinoma)
#' @param n_LOO A double value lower than 1 indicating the percentage of patients to test in the LOO (e.g. 0.3 is 30 percentage)
#' @param k_tuning An integer value greater than 10 indicating the number of training patients to use for the tuning phase (e.g. 10)
#' @param n_cores An integer value greater than 0 indicating the number of cores to use for parallel computing
#' @param th_pred An integer value greater than 0 indicating the power threshold to filter signature PSNs
#' @return A list with elements describe the input data
#' @export
#'
set_variables=function(data_l,
                       dataset_name="simpati_ds",disease_type="NAN",key_words="NAN",
                       n_LOO=0.5,k_tuning=10,th_pred=3,n_cores=2){
  #Extract prepared input data
  geno=data_l$geno
  info=data_l$info
  net=data_l$net
  pathways_l=data_l$pathways_l
  #Keep track of the original genetic profile matrix
  orig_geno=geno
  #Assess the number of elements in the two classes
  tab_status=as.data.frame(table(info[,2]),stringsAsFactors = F);colnames(tab_status)[1]=c("class");
  #Get number of sample per class and their indexes in the profile matrix
  nAs=tab_status$Freq[1];nBs=tab_status$Freq[2];seqAs=seq(1,nAs);seqBs=seq(nAs+1,nrow(info))
  #Simplify the names
  Aletter=substr(tab_status$class[1],1,1);Bletter=substr(tab_status$class[2],1,1);
  #Change name of the patients into more easy2read labels
  isiAs=paste(Aletter,seqAs,sep="");isiBs=paste(Bletter,seq(1,nBs),sep="");
  #Update info
  info=cbind(info,c(isiAs,isiBs));colnames(info)[3]="shortID";colnames(geno)=info[,3];
  #Extract testing patients for the LOO approach
  pat_names=as.character(info$shortID)
  LOO_As=round((nAs*n_LOO))
  LOO_Bs=round((nBs*n_LOO))
  testing_pats=c(sample(pat_names[seqAs],LOO_As),sample(pat_names[seqBs],LOO_Bs))
  #Extract tuning patients from the training set
  left_pat_names=setdiff(pat_names,testing_pats)
  left_nAs=nAs-LOO_As
  left_nBs=nBs-LOO_Bs
  if(left_nAs<k_tuning){tuning_nAs=left_nAs}else{tuning_nAs=k_tuning}
  if(left_nBs<k_tuning){tuning_nBs=left_nBs}else{tuning_nBs=k_tuning}
  tuning_As=sample(left_pat_names[seq(1,left_nAs)],tuning_nAs)
  tuning_Bs=sample(left_pat_names[seq((left_nAs+1),length(left_pat_names))],tuning_nBs)
  tuning_pats=c(tuning_As,tuning_Bs)
  #Assess the thresholds of clustering based on the dimension of the profile matrix
  n_pats=length(pat_names)
  if(n_pats<=50){
    max_th_cl=0.80
    min_th_cl=0.20
  }
  if(n_pats>50 & n_pats<=100){
    max_th_cl=0.75
    min_th_cl=0.25
  }
  if(n_pats>100){
    max_th_cl=0.70
    min_th_cl=0.30
  }

  #Pack everything to return as result
  vars_list=list(geno=geno,info=info,net=net,pathways_l=pathways_l,tab_status=tab_status,
                 nAs=nAs,nBs=nBs,seqAs=seqAs,seqBs=seqBs,Aletter=Aletter,Bletter=Bletter,isiAs=isiAs,isiBs=isiBs,
                 dataset_name=dataset_name,disease_type=disease_type,key_words=key_words,
                 testing_pats=testing_pats,tuning_pats=tuning_pats,
                 max_th_cl=max_th_cl,min_th_cl=min_th_cl,n_cores=n_cores,th_pred=th_pred,
                 orig_geno=orig_geno)
  cat("*Variables ready\n")
  return(vars_list)
}

#' Clean the pathways or feature sets provided by the user
#'
#' Maps the feature names of the matrix of patient's profiles to the corresponding name in the feature sets/pathway
#' Plus, it removes all the feature sets mapping only few existing features (defined by: dim_pathway parameter)
#'
#' @param pathwayList List of gene sets
#' @param geno Matrix of profiles (genes x samples)
#' @param dim_pathway Number greater than 1 indicating the minimum size of the pathways
#' @param n_cores Number greater than 1 indicating the number of cores for parallel computing
#' @return Filtered and indexes pathways/feature sets
#' @import parallel
#' @import doParallel
#' @import foreach
#' @export
#'
clean_pathways = function(pathways_l,geno,size=4,n_cores=2){
  cat("*Cleaning and mapping pathway specific gene sets\n")
  #Start clusters
  if(.Platform$OS.type == "unix") {
    cat(">>Parallel for linux \n")
    cl <- makeCluster(n_cores,type="FORK");
  } else {
    cat(">>Parallel for windows \n")
    cl <- makeCluster(n_cores);
  }
  registerDoParallel(cl);

  #Obtain gene names from gene expression dataset
  gene_names=rownames(geno)

  #For each pathway
  #Compute the number of genes which are also in geno
  count_f = function(x,gene_names){
    sum(!is.na((match(gene_names,x))))
  }

  #Remove pathways with 0 or 1 genes in geno
  netList= parSapply(cl = cl, X=pathways_l, FUN=count_f, gene_names)
  pathways_l=pathways_l[netList>=size]

  #For each pathway
  #Keep the index of the genes in geno that belong to the gene set
  det_f = function(x,gene_names){
    which(gene_names %in% x)
  }
  netList= parSapply(cl = cl, X=pathways_l, FUN=det_f, gene_names)

  #Close clusters and return
  stopCluster(cl)
  cat(">>Pathway list ready to use with:",length(netList),"sets\n")
  return(netList)
}


#' 0-1 matrix standardization
#'
#' 0-1 matrix standardization, the propagation feature values are standardized to reflect
#' how much are high/low with respect all the other features in the matrix
#'
#' @param geno_p Matrix of profiles (genes x samples)
#' @param min Integer number to set the low limit of the range standardization. Default 0
#' @param max Integer number to set the high limit of the range standardization. Default 1
#' @return The same input matrix with standardized values
#' @import scales
#' @export
#'
scaling_matrix=function(geno_p,min=0,max=1,capping=F){
  if(capping==TRUE){
    geno_p=apply(geno_p,2,function(x){
      qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
      caps <- quantile(x, probs=c(.1, .90), na.rm = T)
      H <- 1.5 * IQR(x, na.rm = TRUE)
      x[x < (qnt[1] - H)] <- caps[1]
      x[x > (qnt[2] + H)] <- caps[2]
      return(x)
    })
  }

  #Convert matrix to vector
  geno_p_v=as.vector(geno_p)
  #Scale the vector
  geno_p_v_sc=scales::rescale(geno_p_v, to=c(min,max))
  #Convert vector to matrix with same original dimensions
  geno_p2=matrix(geno_p_v_sc,nrow(geno_p),ncol(geno_p))
  #Apply original labels
  colnames(geno_p2)=colnames(geno_p)
  rownames(geno_p2)=rownames(geno_p)
  #Return
  return(geno_p2)
}

#' 1-100 matrix column-wise ranking
#'
#' column-wise ranks the values and then normalizes them between 1-100
#'
#' @param geno_p Matrix of profiles (genes x samples)
#' @param par Boolean. Default TRUE. If FALSE diseable the parallel computing
#' @return Produce a new matrix s.t. each gene has the rank of its original value in the sample
#' @import scales
#' @export
#'
column_rank = function(geno_p,par=T){
  geno_pR=scales::rescale(Rfast::colRanks(geno_p,parallel = par),to=c(1,100))
  rownames(geno_pR)=rownames(geno_p);colnames(geno_pR)=colnames(geno_p);
  return(geno_pR)
}

#' Find features which separate the patient classes
#'
#' Given a matrix of patient's features, this function ranks the features based on
#' how much their distribution values are different between classes
#'
#' @param geno_p Matrix of profiles (genes x samples)
#' @param Aletter The letter signature of the labels of the first class samples
#' @param Bletter The letter signature of the labels of the second class samples
#' @return Produce a dataframe s.t.  each row correspond to a gene situation
#' the columns are: rank_dist (higher and more it differs between classes)
#' diff_sd (higher and more its variance between classes is different)
#' @import matrixStats
#' @export
#'
get_genes_desc = function(mp,Aletter,Bletter){
  #Compute gene profiling -----
  #Divide the matrix into the patient classes
  seqAs=grep(Aletter,colnames(mp))
  seqBs=grep(Bletter,colnames(mp))
  m1=mp[,seqAs];m2=mp[,seqBs];

  #Compute how much the sd differ between the two groups and rank
  g_sd_p=abs(matrixStats::rowSds(m1)-matrixStats::rowSds(m2))

  #Assess how much a gene separates the two groups based on its values
  min=0.35;max=0.65;rank_dist=rep(0,nrow(m1));power=1;
  #Check if you satisfied the number of genes to find
  for(k in 1:10){
    #Extract the distribution of each row/gene, specifically the 20 and 80Q
    d1=matrixStats::rowQuantiles(m1,probs=c(min,max))
    d2=matrixStats::rowQuantiles(m2,probs=c(min,max))
    for(row_k in 1:nrow(d1)){
      x=d1[row_k,]
      y=d2[row_k,]
      if(x[1]>y[2] | y[1]>x[2]){
        rank_dist[row_k]=power
      }
    }
    power=power+1
    min=min-0.025
    max=max+0.025
  }
  gs_p=data.frame(rank_dist=rank_dist,diff_sd=g_sd_p)
  return(gs_p)
}


#' Vectorizes a matrix
#'
#' This functions takes a matrix and returns a list of its components describing its vectorize form
#' the values, the row and col names and its dimensions
#'
#' @param geno_p Matrix
#' @return List: v = vectorized matrix content, col_names = names of the orginal columns
#' row_names = names of the orginal rows, n_col and n_row contain the number of original columns and rows
#' @export
#'
vectorizing_matrix=function(geno_p){
  #Convert matrix to vector
  geno_p_v=as.vector(geno_p)
  #Extract original labels
  col_names=colnames(geno_p)
  row_names=rownames(geno_p)
  #Extract dimensions
  n_col=ncol(geno_p)
  n_row=nrow(geno_p)
  #Return list of components
  geno_l=list(v=geno_p_v,col_names=col_names,row_names=row_names,n_col=n_col,n_row=n_row)
  return(geno_l)
}

#' Determines the POWER of a patient similarity network
#'
#' Computes how much the two classes separate based on their similarity
#'
#' @param m_sim_cl Similarity matrix
#' @param Aletter Character string. One letter to recognise the first class profiles
#' @param Bletter Character string. One letter to recognise the second class profiles
#' @return Return the power (grade of separation between classes) of the similarity matrix
#' @export
#'
get_power=function(m_sim_cl,Aletter,Bletter){
  get_label=function(m_sim_cl,Aletter,Bletter,min=0.3,max=0.7){
    seqAs=grep(Aletter,colnames(m_sim_cl))
    seqBs=grep(Bletter,colnames(m_sim_cl))
    mAA_sim_cl=m_sim_cl[seqAs,seqAs]
    mBB_sim_cl=m_sim_cl[seqBs,seqBs]
    mAB_sim_cl=m_sim_cl[seqAs,seqBs]

    AAe=quantile(mAA_sim_cl[upper.tri(mAA_sim_cl)],probs = c(min,max))
    BBe=quantile(mBB_sim_cl[upper.tri(mBB_sim_cl)],probs = c(min,max))
    ABe=quantile(mAB_sim_cl,probs = c(min,max))

    label="FAIR"
    if(AAe[1]>BBe[2]){
      if(AAe[1]>ABe[2]){
        label=c(Aletter,AAe[1],BBe[2],ABe[2])
      }
    }

    if(BBe[1]>AAe[2]){
      if(BBe[1]>ABe[2]){
        label=c(Bletter,BBe[1],AAe[2],ABe[2])
      }
    }
    return(label)
  }

  min_th=0.350
  max_th=0.650
  min_lvs=seq(min_th,0.1,-0.025)
  max_lvs=seq(max_th,0.9,0.025)
  power_lvs=cbind(min_lvs,max_lvs)
  power_lvs=power_lvs[-c(10),]
  colnames(power_lvs)=c("min","max");rownames(power_lvs)=seq(1,nrow(power_lvs))

  info_df=data.frame(class="FAIR",power=0,minSIGN=0,maxWEAK1=0,maxWEAK2=0,stringsAsFactors = F)
  power=best_distr=0;fin_label=label="FAIR";
  for(lv in 1:nrow(power_lvs)){
    v_info=get_label(m_sim_cl,Aletter,Bletter,min=power_lvs[lv,1],max=power_lvs[lv,2])
    label=v_info[1]
    e_distr=v_info[c(2,3,4)]
    if(label!="FAIR"){
      power=lv
      fin_label=label
      best_distr=as.numeric(e_distr)
    }
    label="FAIR"
  }
  info=c(fin_label,power,best_distr)
  if(power!=0){
    info_df[1,]=info
    info_df$power=as.numeric(info_df$power)
    info_df$minSIGN=as.numeric(info_df$minSIGN)
    info_df$maxWEAK1=as.numeric(info_df$maxWEAK1)
    info_df$maxWEAK2=as.numeric(info_df$maxWEAK2)
  }
  return(info_df)
}

#' Convert minutes time value in hours
#'
#' @param raw_min Double value
#' @return New double value formatted in hours
#' @export
#'
h_min=function(raw_min){
  if(raw_min>=60){
    h=floor(raw_min/60)
    mins=round(raw_min-(60*h))
    read_run_time=as.numeric(paste(h,mins,sep="."))
  }else{
    h=0
    mins=round(raw_min)
    if(mins>9){
      mins=paste(h,mins,sep=".")
    }else{
      mins=paste(h,".0",mins,sep="")
    }
    read_run_time=as.numeric(mins)
  }
  return(read_run_time)
}

get_top_lv=function(PSN_power_df){
  med_spe_lvs=0
  lvs=unique(PSN_power_df$level)
  for(lv in lvs){
    med_spe_lv=median(PSN_power_df[PSN_power_df$level>=lv,]$auroc)
    med_spe_lvs=c(med_spe_lvs,med_spe_lv)
  }

  med_spe_lvs=med_spe_lvs[-1]
  top_lv=lvs[which.max(med_spe_lvs)]
  return(top_lv)
}

#' Determines testing patient probability of belonging to the first patient class based on BFC algorithm
#'
#' Assumes that the testing patient belongs to the first class, it iterates the Best Friend Connector
#' algorithm to find how much it fits in the representative members of the class
#'
#' @param tepat String character of the name of the testing profile
#' @param m_sim Matrix: Original PSN with all the profiles both training and testing
#' @param m_sim_cl Matrix: Clustering PSN of only training profiles
#' @param vars_list Simpati object: List of variables
#' @param best_direction Power info result
#' @return A dataframe of information for the profile class prediction, precisely:
#' Atop: the top percentage in which the profile is included
#' AAe: the median of the profile's AA edges
#' DIFFe: the difference beteen the medians of profile's AA and AB edges
#' assumption: the assumption made with the function: ASign: the profile is fake of first
#' and signature class
#' @export
#'
get_A_prob=function(tepat,m_sim,m_sim_cl,vars_list,best_direction){
  #Create a fake known profile for the testing patient in the first class
  tepat_art=paste(vars_list$Aletter,"999",sep="")
  #Get the position of the running testing patient
  indx_tepat=match(tepat,colnames(m_sim))
  #Extract its values of similarity with respect the other patients
  v_tepat=as.data.frame(m_sim[,indx_tepat]);colnames(v_tepat)=rownames(v_tepat)[indx_tepat]=tepat_art
  v_tepat2=as.data.frame(v_tepat[order(rownames(v_tepat)),1])
  rownames(v_tepat2)=rownames(v_tepat)[order(rownames(v_tepat))]
  colnames(v_tepat2)=colnames(v_tepat)
  #Let's reintegrate the testing patient is the PSN as one of the first class
  diag(m_sim_cl)=1
  tmp=merge(v_tepat2,m_sim_cl,by="row.names",sort=F)
  rownames(tmp)=tmp[,1];tmp=tmp[,-1]
  tmp2=merge(v_tepat2,t(tmp),by="row.names",sort=F)
  rownames(tmp2)=tmp2[,1];tmp2=tmp2[,-1]
  xA=tmp2[,rownames(tmp2)]
  xA=as.matrix(xA)

  check_exist=F;
  th_min_test=0;th_max_test=1;best_THin=1;best_PSNin=xA

  for(run1 in 1:9){
    if(best_direction[1]==vars_list$Aletter){
      tmp1=strong.comm.AA.cor(mat=xA,nAs=length(grep(vars_list$Aletter,colnames(xA))),thW=1,thGroup=th_max_test,strong = T)
      check_exist=tepat_art %in% colnames(tmp1)
      if(check_exist){
        best_THin=th_max_test
        best_PSNin=tmp1
        th_max_test=th_max_test-0.1
        check_exist=F
      }else{
        break;
      }

    }else{
      tmp1=strong.comm.AA.cor(mat=xA,nAs=length(grep(vars_list$Aletter,colnames(xA))),thW=1,thGroup=th_min_test,strong = F)
      check_exist=tepat_art %in% colnames(tmp1)
      if(check_exist){
        best_THin=1-th_min_test
        best_PSNin=tmp1
        th_min_test=th_min_test+0.1
        check_exist=F
      }else{
        break;
      }
    }

    if(sum(table(substr(colnames(tmp1),1,1))==2)>=1){
      break;
    }
  }

  if(best_direction[1]==vars_list$Aletter){
    nAs_in=length(grep(vars_list$Aletter,colnames(best_PSNin)))

    te_AAe=best_PSNin[1:nAs_in,tepat_art]
    if(sum(te_AAe!=0)>0){
      te_AAe=round(median(te_AAe[te_AAe!=0]),4)
    }else{
      te_AAe=0
    }

    te_ABe=best_PSNin[(nAs_in+1):(nrow(best_PSNin)),tepat_art]
    if(sum(te_AAe!=0)>0){
      te_ABe=round(median(te_ABe[te_ABe!=0]),4)
    }else{
      te_ABe=0
    }

    te_ee=data.frame(Atop=best_THin,AAe=te_AAe,ABe=te_ABe,
                     assumption="ASign",stringsAsFactors = F)
  }else{
    nAs_in=length(grep(vars_list$Aletter,colnames(best_PSNin)))
    te_AAe=best_PSNin[1:nAs_in,tepat_art]
    if(sum(te_AAe!=0)>0){
      te_AAe=round(median(te_AAe[te_AAe!=0]),4)
    }else{
      te_AAe=0
    }

    te_ABe=best_PSNin[(nAs_in+1):(nrow(best_PSNin)),tepat_art]
    if(sum(te_ABe!=0)>0){
      te_ABe=round(median(te_ABe[te_ABe!=0]),4)
    }else{
      te_ABe=0
    }
    te_ee=data.frame(Atop=best_THin,AAe=te_AAe,ABe=te_ABe,
                     assumption="AWeak",stringsAsFactors = F)
  }

  return(te_ee)
}

#' Determines testing patient probability of belonging to the second patient class based on BFC algorithm
#'
#' Assumes that the testing patient belongs to the second class, it iterates the Best Friend Connector
#' algorithm to find how much it fits in the representative members of the class
#'
#' @param tepat String character of the name of the testing profile
#' @param m_sim Matrix: Original PSN with all the profiles both training and testing
#' @param m_sim_cl Matrix: Clustering PSN of only training profiles
#' @param vars_list Simpati object: List of variables
#' @param best_direction Power info result
#' @return A dataframe of information for the profile class prediction, precisely:
#' Btop: the top percentage in which the profile is included
#' BBe: the median of the profile's BB edges
#' DIFFe: the difference beteen the medians of profile's BB and AB edges
#' assumption: the assumption made with the function: BSign: the profile is fake of second
#' and signature class
#' @export
#'
get_B_prob=function(tepat,m_sim,m_sim_cl,vars_list,best_direction){
  tepat_art=paste(vars_list$Bletter,"999",sep="")
  #Get the position of the running testing patient
  indx_tepat=match(tepat,colnames(m_sim))
  #Extract its values of similarity with respect the other patients
  v_tepat=as.data.frame(m_sim[,indx_tepat]);colnames(v_tepat)=rownames(v_tepat)[indx_tepat]=tepat_art
  v_tepat2=as.data.frame(v_tepat[order(rownames(v_tepat)),1])
  rownames(v_tepat2)=rownames(v_tepat)[order(rownames(v_tepat))]
  colnames(v_tepat2)=colnames(v_tepat)

  #Let's reintegrate the testing patient is the PSN as one of the first class
  diag(m_sim_cl)=1
  tmp=merge(v_tepat2,m_sim_cl,by="row.names",sort=F)
  rownames(tmp)=tmp[,1];tmp=tmp[,-1]
  tmp2=merge(v_tepat2,t(tmp),by="row.names",sort=F)
  rownames(tmp2)=tmp2[,1];tmp2=tmp2[,-1]
  xB=tmp2[,rownames(tmp2)]
  xB=as.matrix(xB)

  check_exist=F;
  th_min_test=0;th_max_test=1;best_THin=1;best_PSNin=xB

  for(run1 in 1:9){
    if(best_direction[1]==vars_list$Aletter){
      tmp1=strong.comm.BB.cor(mat=xB,nAs=length(grep(vars_list$Aletter,colnames(xB))),thW=1,thGroup=th_min_test,strong = F)
      check_exist=tepat_art %in% colnames(tmp1)
      if(check_exist){
        best_THin=1-th_min_test
        best_PSNin=tmp1
        th_min_test=th_min_test+0.1
        check_exist=F
      }else{
        break;
      }

    }else{
      tmp1=strong.comm.BB.cor(mat=xB,nAs=length(grep(vars_list$Aletter,colnames(xB))),thW=1,thGroup=th_max_test,strong = T)
      check_exist=tepat_art %in% colnames(tmp1)
      if(check_exist){
        best_THin=th_max_test
        best_PSNin=tmp1
        th_max_test=th_max_test-0.1
        check_exist=F
      }else{
        break;
      }

    }

    if(sum(table(substr(colnames(tmp1),1,1))==2)>=1){
      break;
    }

  }

  if(best_direction[1]==vars_list$Bletter){
    nAs_in=length(grep(vars_list$Aletter,colnames(best_PSNin)))
    te_BBe=best_PSNin[(nAs_in+1):(nrow(best_PSNin)),tepat_art]
    te_BBe=round(median(te_BBe[te_BBe!=0]),4)
    te_ABe=best_PSNin[1:nAs_in,tepat_art]
    te_ABe=round(median(te_ABe[te_ABe!=0]),4)
    te_ee=data.frame(Btop=best_THin,BBe=te_BBe,ABe=te_ABe,
                     assumption="BSign",stringsAsFactors = F)
  }else{
    nAs_in=length(grep(vars_list$Aletter,colnames(best_PSNin)))
    te_BBe=best_PSNin[(nAs_in+1):(nrow(best_PSNin)),tepat_art]
    te_BBe=round(median(te_BBe[te_BBe!=0]),4)
    te_ABe=best_PSNin[1:nAs_in,tepat_art]
    te_ABe=round(median(te_ABe[te_ABe!=0]),4)
    te_ee=data.frame(Btop=best_THin,BBe=te_BBe,ABe=te_ABe,
                     assumption="BWeak",stringsAsFactors = F)
  }

  return(te_ee)
}

#' Determines testing patient class based on the strength of its similarities
#'
#' Determines testing patient class based on the strength of its similarities with the training
#' patients that populate the two classes in comparison
#'
#' @param tepat String character of the name of the testing profile
#' @param m_sim Matrix: Original PSN with all the profiles both training and testing
#' @param m_sim_cl Matrix: Clustering PSN of only training profiles
#' @param vars_list Simpati object: List of variables
#' @param top_th Percentage of profiles of each class to compare with testing profile
#' @return A dataframe of information for the profile class prediction, precisely:
#' A_strength: how much is WJ similar to A profiles similarities
#' B_strength: how much is WJ similar to B profiles similarities
#' @export
#'
get_strength_sim=function(tepat,m_sim,m_sim_cl,vars_list,top_th=0.30){
  m_sim_sim=sim_WJ(m_sim)
  #Get the position of the running testing patient
  indx_tepat=match(tepat,colnames(m_sim_sim))
  #Extract its values of similarity with respect the other patients
  v_tepat=as.data.frame(m_sim_sim[,indx_tepat]);colnames(v_tepat)=tepat
  v_tepat=v_tepat[-indx_tepat,]
  names(v_tepat)=colnames(m_sim_sim)[-indx_tepat]
  #Keep only the patients in the crews
  v_tepat=v_tepat[colnames(m_sim_cl)]
  Asims=v_tepat[grep(vars_list$Aletter,names(v_tepat))]
  Bsims=v_tepat[grep(vars_list$Bletter,names(v_tepat))]
  #Get how similar it is to his best friends As and Bs
  A_strength=mean(Asims[order(Asims, decreasing=TRUE)[1:(length(Asims)*top_th)]])
  B_strength=mean(Bsims[order(Bsims, decreasing=TRUE)[1:(length(Bsims)*top_th)]])
  #Compose result
  res=data.frame(A_strength=A_strength,B_strength=B_strength)
  return(res)
}

#' Determines testing patient class based on the strength of its similarities
#'
#' Determines testing patient class based on the strength of its similarities with the training
#' patients that populate the two classes in comparison
#'
#' @param PSNs_info String character of the name of the testing profile
#' @param k1 Matrix: Original PSN with all the profiles both training and testing
#' @param th_v Matrix: Clustering PSN of only training profiles
#' @param k_feat Simpati object: List of variables
#' @param comp Percentage of profiles of each class to compare with testing profile
#' @return A dataframe of information for the profile class prediction, precisely:
#' A_strength: how much is WJ similar to A profiles similarities
#' B_strength: how much is WJ similar to B profiles similarities
#' @export
#'
get_indxs=function(PSNs_info,k1,th_v,k_feat,comp){
  conc_th=th_v[k1];name_conc_th=names(th_v)[k1]
  feat=PSNs_info[,k_feat]

  if(comp=="=="){
    indxs=(feat==conc_th)
  }
  if(comp=="<="){
    indxs=(feat<=conc_th)
  }
  if(comp==">="){
    indxs=(feat>=conc_th)
  }

  if(sum(indxs)>1){
    PSNs_info0=PSNs_info[indxs,]
    n_classes1=dim(table(substr(PSNs_info0$testing_pats,1,1)))
    n_classes2=dim(table(PSNs_info0$pred_cl))
    if(n_classes1<2 | n_classes2<2){
      return(F)
    }else{
      return(indxs)
    }
  }else{
    return(F)
  }
}
