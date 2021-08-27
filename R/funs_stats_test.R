#' Performs the last step of classification by predicting the final class of the testing profiles
#'
#' Performs the consensus prediction of the testing patients and retrieves the classification performances
#'
#' @param df matrix of PSN information per testing patient profile
#' @param node_profiles_db data frame of profiles
#' @return A list with: the testing profiles, the correct a prior classes, the predicted classes and the performances
#' @import ROCR
#' @import stats
#' @export
#'
get_perfs = function(df){
  #Get class variables -----
  Aletter=names(table(substr(df$testing_pat,1,1)))[1]
  Bletter=names(table(substr(df$testing_pat,1,1)))[2]
  freq_cls=as.data.frame(table(df$class))

  if(nrow(freq_cls)==1){
    max_freq_cl=as.character(freq_cls$Var1[which.max(freq_cls$Freq)])
    min_freq_cl=setdiff(c(Aletter,Bletter),max_freq_cl)
  }else{
    max_freq_cl=as.character(freq_cls$Var1[which.max(freq_cls$Freq)])
    min_freq_cl=as.character(freq_cls$Var1[which.min(freq_cls$Freq)])
  }

  #Get info about profiles and initiate a fake predicted label vector -----
  pats=unique(df$testing_pat)
  cls=substr(pats,1,1)
  preds=rep("C",length(cls))
  n_trs=nrow(node_profiles_db)

  #Get the behaviour of each testing profile ----
  probs_v=seq(0,1,0.01)
  for(k_pat in 1:length(pats)){
    pat=pats[k_pat]
    pat_df=df[df$testing_pat==pat,]
    Aqs=quantile(pat_df$minSIGN,probs=probs_v)
    Bqs=quantile(pat_df$maxWEAK1,probs=probs_v)
    Cqs=quantile(pat_df$maxWEAK2,probs=probs_v)
    qs=c(Aqs,Bqs,Cqs)
    if(k_pat==1){
      qs_df=qs
    }else{
      qs_df=rbind(qs_df,qs)
    }
  }
  rownames(qs_df)=pats

  #Compute how much the behaviour of the testing profile is similar to
  #randomly generated profiles of tuning -----
  colnames(qs_df)=colnames(node_profiles_db)
  qs_df2=rbind(node_profiles_db,qs_df)
  qs_df2=t(qs_df2)
  #PSN_COR=Hmisc::rcorr(qs_df2, type="pearson")$r[-seq(1,10),seq(1,10)]
  PSN_COR=cor(qs_df2, method = "pearson")[-seq(1,n_trs),seq(1,n_trs)]

  #Get the predicted class of the testing profile based on its strongest similar tuning profile ----
  labels=substr(colnames(PSN_COR)[apply(PSN_COR,1,which.max)],1,12)
  preds=ifelse(labels=="most_freq_cl", max_freq_cl, min_freq_cl)

  #Get performances -----
  corr_labs=ifelse(cls==Aletter, 0, 1)
  pred_labs=ifelse(preds==Aletter, 0, 1)

  pred <- ROCR::prediction(pred_labs,corr_labs)
  auroc=ROCR::performance(pred, "auc")@y.values[[1]]
  aupr=ROCR::performance(pred, "aucpr")@y.values[[1]]
  cat("auroc:",auroc," aupr:",aupr,"\n")

  res_perfs=list(testing_prof=pats,correct_classes=cls,predicted_classes=preds,auroc=auroc,aupr=aupr)
  return(res_perfs)
}
