#' Wrapper of Simpati workflow for a human gene study
#'
#' This function performs all Simpati operations, both classification and pathway enrichment, on a study which includes
#' genes as features describing the human patient's profiles. This function is intended to help users which do not want to go
#' through all the Simpati operations manually and their study falls into the wrapper application.
#'
#' @param geno A numeric matrix of patient profiles: features (e.g. genes) x samples
#' @param info A character matrix of info about the patients: patient_ID | class (e.g. clinical information)
#' @param net Feature interaction network as adjacency matrix (e.g. gene-gene interactio network)
#' @param pathways_l List of feature sets (e.g. gene lists)
#' @param dataset_name The name of the dataset or project
#' @param disease_type The DISGNET semantic type of the sample condition: "Neoplastic Process" or "Disease or Syndrome"
#' @param key_words A character vector of one or more words related to the sample condition (e.g. KIRC: Kidney, Carcinoma)
#' @param n_cores Default 5, An integer value greater than 0 indicating the number of cores to use for parallel computing
#' @param n_LOO Default 0.3, A double value lower than 1 indicating the percentage of patients to test in the LOO (e.g. 0.3 is 30 percentage)
#' @param seed Default 0, An integer value to set the random seed generator of the run
#' @param test_run Default TRUE, boolean indicating if to run the function with a limited numbers of sets from pathways_l to check if Simpati works properly on the input data
#' @return A list of results
#' @import data.table
#' @import Matrix
#' @import scales
#' @import parallel
#' @import doParallel
#' @import Rfast
#' @import igraph
#' @import stats
#' @import foreach
#' @import FNN
#' @import ROCR
#' @import ggplot2
#' @import matrixStats
#' @import plyr
#' @export
#'
wrapper_human_mutations=function(geno,info,net,pathways_l,dataset_name,disease_type,key_words,n_cores=5,n_LOO=0.3,seed=0,test_run=TRUE){
  #Set seed for reproduce the results
  set.seed(seed)

  #Reduce pathways in case you are running a test run ----
  if(test_run){
    pathways_l=pathways_l[1:1000]
    pathways_l=pathways_l[grep("source-MSIGDB",names(pathways_l))]
  }

  #Load the databases ----
  #human protein atlas for normal cell lines needed for the tissue enrichment
  data(db_normal_ATLAS)

  #human protein atlas for cancer cell lines needed for the tumor enrichment
  data(db_cancer_ATLAS)

  #disgnet needed for the disease-feature (aka gene) enrichment
  data(db_disgnet)

  #Data preparation
  data_l=data_preparation(geno,info,net,pathways_l,n_cores)
  #Clean
  rm(geno,info)

  #Generate the object list with all the variables needed by the method
  vars_l=set_variables(data_l,dataset_name,disease_type,key_words,n_cores=n_cores,n_LOO=n_LOO)
  #Clean
  rm(disease_type,key_words)

  #Start to measure the running time
  start_time <- Sys.time()

  #Propagation
  vars_l[["geno_p"]]=rwrProp(net=vars_l$net, geno=vars_l$geno, n_cores=vars_l$n_cores);

  #Prepare pathways
  vars_l[["pathways_indxs_l"]]=clean_pathways(pathways_l=vars_l$pathways_l, geno=vars_l$geno_p, n_cores=vars_l$n_cores)

  #Predicts testing profiles
  cat("*Making pathway specific predictions of testing profiles \n")
  PSNs_info=predict_pats(vars_l);

  #Predict testing profiles and get classification performances
  res=get_perfs(PSNs_info)

  #Stop running time of classification
  end_time <- Sys.time()
  running_time=as.numeric(difftime(end_time, start_time, units = "mins")[[1]])

  #Save results of classification and pathway analysis
  cat("Ended classification in",running_time,"\n")

  cat("Starting pathway enrichment with classification results \n")

  #Get the updated info dataframe
  new_info=get_new_info(res,vars_l)
  #Get the pathways used in the classification
  pathways_used_l=get_used_pathways(PSNs_info,vars_l)

  #Data preparation
  data_l=data_preparation(vars_l$orig_geno[,new_info$patientID],new_info,vars_l$net,pathways_used_l)

  #Generate the object list with all the variables needed by the method
  pr_vars_l=set_variables(data_l,vars_l$dataset_name,vars_l$disease_type,vars_l$key_words,n_cores=vars_l$n_cores)

  #Propagation
  pr_vars_l[["geno_p"]]=rwrProp(net=pr_vars_l$net,geno=pr_vars_l$geno,n_cores=pr_vars_l$n_cores);

  #Prepare pathways
  pr_vars_l[["pathways_used_indxs_l"]]=clean_pathways(pathways_l=pr_vars_l$pathways_l,geno=pr_vars_l$geno_p,n_cores=pr_vars_l$n_cores)

  #Get signature pathways ----
  PSN_data_l=get_sign_pathways(pr_vars_l$geno_p, pr_vars_l$pathways_used_indxs_l, pr_vars_l, n_cores=pr_vars_l$n_cores)

  #Get enrichment ----
  gene_sets_l=sapply(PSN_data_l$PSN_comp_l,function(x){x$m_prof_l$row_names})
  if(length(gene_sets_l)>1){
    enrXpathway_l=get_pathway_enrich(gene_sets_l,pr_vars_l$key_words,pr_vars_l$disease_type,
                                     db_disgnet,db_cancer_ATLAS,db_normal_ATLAS,
                                     n_cores=3, records_only = T)
  }else{enrXpathway_l=list()}

  cat("Saving results and workflow data \n")

  #Save classification and enrichment results -----
  PSN_enr_df=PSN_data_l[["PSN_info_df"]]
  enr_PSN_df=enrXpathway_l[["enrXpathway_df"]]
  if(!is.null(enr_PSN_df)){
    PSN_enr_df=merge(PSN_enr_df,enr_PSN_df,by="pathway_name",all.x=T)
  }

  res_l=list(classification_res=res,PSNs_info=PSNs_info,running_time=running_time,PSN_enr_df=PSN_enr_df,
             vars_l=vars_l,pr_vars_l=pr_vars_l,PSN_data_l=PSN_data_l,enrXpathway_l=enrXpathway_l)
  return(res_l)
}

