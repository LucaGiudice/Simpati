#' Determine the consensus class of a testing patient from a run of classification
#'
#' Given the classes assigned to a testing patient by all the signature PSN during the classification phase, it
#' determines the final class based on the frequency
#'
#' @param df Matrix which list the signature pathways with their corresponding predicted class per testing patient
#' @return Produce a new matrix s.t. each gene has the rank of its original value in the sample
#' @export
#'
get_consensus_pred = function(df){
  #Get consensus predictions
  testing_patients=unique(df$testing_pat)
  for(k_pat in 1:length(testing_patients)){
    testing_patient=testing_patients[k_pat]
    freq_df=table(df$prediction[df$testing_pat==testing_patient])
    if(k_pat==1){
      cons_preds=names(freq_df)[which.max(freq_df)]
    }else{
      cons_preds=c(cons_preds,names(freq_df)[which.max(freq_df)])
    }
  }
  names(cons_preds)=testing_patients
  return(cons_preds)
}

#' Determine the consensus class of a testing patient from a run of classification
#'
#' Undefined
#'
#' @param PSNs_info Undefined
#' @param comb Undefined
#' @return Produce a new matrix s.t. each gene has the rank of its original value in the sample
#' @export
#'
wrapper_cons = function(PSNs_info,comb){
  indxs=match(names(comb)[1:(ncol(comb)-2)],colnames(PSNs_info))
  for(k in 1:length(indxs)){
    if(k==1){
      sub_df=PSNs_info[PSNs_info[,indxs[k]]>=comb[1,k],]
    }else{
      if(indxs[k]!=19){
        sub_df=sub_df[sub_df[,indxs[k]]>=comb[1,k],]
      }else{
        sub_df=sub_df[sub_df[,indxs[k]]<=comb[1,k],]
      }
    }
    if(nrow(sub_df)==0){return(NA)}
  }

  if(nrow(sub_df)>1){
    cons_pred=get_consensus_pred(sub_df)
    return(cons_pred)
  }else{
    return(NA)
  }
}

#' Determine the consensus class of a testing patient from a run of classification
#'
#' Undefined
#'
#' @param PSNs_info Undefined
#' @param vars_l Undefined
#' @return Produce a new matrix s.t. each gene has the rank of its original value in the sample
#' @import ROCR
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
#'
tune_performances = function(PSNs_info,vars_l){
  n_cores=8
  #Add the column of right/wrong prediction -----
  Aletter=names(table(substr(PSNs_info$testing_pat,1,1)))[1]
  Bletter=names(table(substr(PSNs_info$testing_pat,1,1)))[2]
  n_te_pats=length(unique(PSNs_info$testing_pat))

  corr_label=as.integer(substr(PSNs_info$testing_pat,1,1)==PSNs_info$prediction)
  cl_pred_label=corr_label
  cl_pred_label[cl_pred_label==1 & substr(PSNs_info$prediction,1,1)==Bletter]=2
  cl_pred_label[cl_pred_label==0 & substr(PSNs_info$prediction,1,1)==Bletter]=-2
  cl_pred_label[cl_pred_label==0 & substr(PSNs_info$prediction,1,1)==Aletter]=-1

  PSNs_info$cor_label=substr(PSNs_info$testing_pat,1,1)
  PSNs_info$bin_label=corr_label
  PSNs_info$cl_pred=cl_pred_label

  #The most frequent PSN class provides better prediction
  PSNsXcl=data.frame(table(PSNs_info$class),stringsAsFactors = F)
  low_freq_cl=as.character(PSNsXcl$Var1[which.min(PSNsXcl$Freq)])
  PSNs_info=PSNs_info[PSNs_info$class!=low_freq_cl,]

  #Two variables discriminate with opposite trends good predicitions for the two classes
  PSNs_info$MAX_strength=abs(scale(PSNs_info$MAX_strength)[,1])
  PSNs_info$SSe=abs(scale(PSNs_info$SSe)[,1])

  #Set variables to create the combination of selections s.t. they will be tested because
  #in one lies the best performances ----
  v_powers=unique(PSNs_info$power)
  v_powers=v_powers[order(v_powers,decreasing = F)]
  names(v_powers)=v_powers

  qs=seq(0.1,0.9,0.1)
  v_max_strength=quantile(PSNs_info$MAX_strength,probs=qs);names(v_max_strength)=qs
  v_min_top=quantile(PSNs_info$MIN_top,probs=qs);names(v_min_top)=qs
  v_diff_strength=quantile(PSNs_info$DIFF_strength,probs=qs);names(v_diff_strength)=qs
  v_diff_top=quantile(PSNs_info$DIFF_top,probs=qs);names(v_diff_top)=qs
  v_sse=quantile(PSNs_info$SSe,probs=qs);names(v_sse)=qs

  #Create combinations
  combs_df=expand.grid(v_powers,v_max_strength,v_min_top,v_diff_strength,v_diff_top,v_sse,0,0)
  colnames(combs_df)=c("power","MAX_strength","MIN_top","DIFF_strength","DIFF_top","SSe","auroc","aupr")

  name_combs_df=expand.grid(names(v_powers),names(v_max_strength),names(v_min_top),
                            names(v_diff_strength),names(v_diff_top),names(v_sse),0,0,stringsAsFactors = F)
  colnames(name_combs_df)=c("power","MAX_strength","MIN_top","DIFF_strength","DIFF_top","SSe","auroc","aupr")
  #Divide in chunks
  combs_df_l=split.data.frame(combs_df,
                              cut(seq_len(nrow(combs_df)),
                                  pretty(seq_len(nrow(combs_df)), n_cores)))
  name_combs_df_l=split.data.frame(name_combs_df,
                              cut(seq_len(nrow(name_combs_df)),
                                  pretty(seq_len(nrow(name_combs_df)), n_cores)))

  #One core one chunk of pathways to process ----
  if(.Platform$OS.type == "unix") {
    cat("Parallel for linux \n")
    cl <- makeCluster(n_cores,type="FORK");
  } else {
    cat("Parallel for windows \n")
    cl <- makeCluster(n_cores);
  }
  registerDoParallel(cl);

  combs_chunks_res_l=list();
  combs_chunks_res_l=foreach(k_chunk=1:length(combs_df_l),.inorder=FALSE,.noexport=c(),
                             .packages = c("ROCR"),.export = c("get_consensus_pred",
                                                               "wrapper_cons")) %dopar% {
    combs_df_chunk=combs_df_l[[k_chunk]]
    name_combs_df_chunk=name_combs_df_l[[k_chunk]]

    cat(k_chunk,"on",length(combs_df_l),"\n")
    for(k_comb in 1:nrow(combs_df_chunk)){
      cat(k_comb,"on",nrow(combs_df_chunk),"\n")
      comb=combs_df_chunk[k_comb,]
      name_comb=name_combs_df_chunk[k_comb,]
      cons_preds=wrapper_cons(PSNs_info,comb)
      print(cons_preds);cat("\n")

      if(length(cons_preds)!=1){
        corr_labs=as.vector(substr(names(cons_preds),1,1))
        pred_labs=as.vector(cons_preds)

        corr_labs=ifelse(corr_labs==Aletter, 0, 1)
        pred_labs=ifelse(pred_labs==Aletter, 0, 1)

        if(length(unique(cons_preds))!=1 & length(pred_labs)==n_te_pats){
          pred <- ROCR::prediction(pred_labs,corr_labs)
          auroc=ROCR::performance(pred, "auc")@y.values[[1]]
          aupr=ROCR::performance(pred, "aucpr")@y.values[[1]]
        }else{
          auroc=aupr=0
        }
      }else{
        auroc=aupr=0
      }

      name_comb$auroc=auroc;name_comb$aupr=aupr
      name_combs_df_chunk[k_comb,]=name_comb
    }

    return(name_combs_df_chunk)
  }
  combs_df=data.table::rbindlist(combs_chunks_res_l)
  class(combs_df)="data.frame"
  stopCluster(cl)
  return(combs_df)
}

#' Determine the consensus class of a testing patient from a run of classification
#'
#' Undefined
#'
#' @param PSNs_info Undefined
#' @param vars_l Undefined
#' @param best_comb Undefined
#' @return Produce a new matrix s.t. each gene has the rank of its original value in the sample
#' @import ROCR
#' @export
#'
get_consensus_classes = function(PSNs_info,vars_l,best_comb){
  #Add the column of right/wrong prediction -----
  Aletter=names(table(substr(PSNs_info$testing_pat,1,1)))[1]
  Bletter=names(table(substr(PSNs_info$testing_pat,1,1)))[2]
  n_te_pats=length(unique(PSNs_info$testing_pat))

  corr_label=as.integer(substr(PSNs_info$testing_pat,1,1)==PSNs_info$prediction)
  cl_pred_label=corr_label
  cl_pred_label[cl_pred_label==1 & substr(PSNs_info$prediction,1,1)==Bletter]=2
  cl_pred_label[cl_pred_label==0 & substr(PSNs_info$prediction,1,1)==Bletter]=-2
  cl_pred_label[cl_pred_label==0 & substr(PSNs_info$prediction,1,1)==Aletter]=-1

  PSNs_info$cor_label=substr(PSNs_info$testing_pat,1,1)
  PSNs_info$bin_label=corr_label
  PSNs_info$cl_pred=cl_pred_label

  #The most frequent PSN class provides better prediction
  PSNsXcl=data.frame(table(PSNs_info$class),stringsAsFactors = F)
  low_freq_cl=as.character(PSNsXcl$Var1[which.min(PSNsXcl$Freq)])
  PSNs_info=PSNs_info[PSNs_info$class!=low_freq_cl,]

  #Two variables discriminate with opposite trends good predicitions for the two classes
  PSNs_info$MAX_strength=abs(scale(PSNs_info$MAX_strength)[,1])
  PSNs_info$SSe=abs(scale(PSNs_info$SSe)[,1])

  best_comb$power=as.numeric(best_comb$power)
  best_comb$MAX_strength=quantile(PSNs_info$MAX_strength,probs = as.numeric(best_comb$MAX_strength))
  best_comb$MIN_top=quantile(PSNs_info$MIN_top,probs = as.numeric(best_comb$MIN_top))
  best_comb$DIFF_strength=quantile(PSNs_info$DIFF_strength,probs = as.numeric(best_comb$DIFF_strength))
  best_comb$DIFF_top=quantile(PSNs_info$DIFF_top,probs = as.numeric(best_comb$DIFF_top))
  best_comb$SSe=quantile(PSNs_info$SSe,probs = as.numeric(best_comb$SSe))

  cons_preds=wrapper_cons(PSNs_info,best_comb)

  corr_labs=as.vector(substr(names(cons_preds),1,1))
  pred_labs=as.vector(cons_preds)

  corr_labs=ifelse(corr_labs==Aletter, 0, 1)
  pred_labs=ifelse(pred_labs==Aletter, 0, 1)

  pred <- ROCR::prediction(pred_labs,corr_labs)
  auroc=ROCR::performance(pred, "auc")@y.values[[1]]
  aupr=ROCR::performance(pred, "aucpr")@y.values[[1]]

  perf=data.frame(auroc=auroc,aupr=aupr)
  res=list(consensus_classes=cons_preds,peformances=perf)
  return(res)
}








