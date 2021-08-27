#' Compute Trending Matching similarity
#'
#' From the two components of the similarity, it computes the overall score also called as TM in the paper
#'
#' @param v vector of two values between 0 and 1, each value is a component of the TM similarity
#' @return Trending Matching similarity
#' @export
#'
autri = function(v,w=c(1,1.4,2)){
  diff=1-(abs(v[1]-v[2]))
  v=c(v[1],v[2],diff)
  y=sum(v*w)
  return(y)
}

#' Compute the Trending Matching similarity per each gene
#'
#' Determine which is the similarity between patients using only one gene at the time to evaluate
#' how much each gene is contributing to the final patient similarities
#'
#' @param input_m2 Matrix of profiles (genes x samples)
#' @return Produce a similarity matrix per gene/row of the original matrix
#' @import data.table
#' @import matrixStats
#' @import scales
#' @import stats
#' @export
#'
Adj_Weighted_Jaccard_gene <- function(input_m2) {
  # require(data.table)
  # require(matrixStats)

  for(k in 1:nrow(input_m2)){
    input_m=input_m2[k,]
    dim(input_m)=c(1,ncol(input_m2))
    colnames(input_m)=colnames(input_m2)
    rownames(input_m)=rownames(input_m2)[k]

    input_m <- na.omit(input_m)

    col <- colnames(input_m)
    row <- rownames(input_m)
    nrow <- nrow(input_m)
    ncol <- ncol(input_m)

    samples <- 1:ncol(input_m)
    comb <- data.table::CJ(samples, samples)
    comb[, i := .I]
    comb <- data.table::melt(comb, 'i')
    data.table::setorder(comb, value)
    v2 <- paste0("V", 1:2)
    comb[, variable2 := v2 , keyby = i]
    comb2 <- data.table::dcast(comb, i ~ variable2, value.var = 'value')
    combUnique <- unique(comb2, by = c('V1', 'V2'))   #creation of all combination of genes

    XXcomb=combUnique
    XX <- apply(combUnique[, -'i'], 1, function(x) {
      x2 <- rowRanges(input_m, cols = x)
      s <- colSums2(x2)
      s[1] / s[2]
    })

    data.table::set(XXcomb, j = 'xx', value = XX)
    rez2 <- merge(comb2, XXcomb[, -'i'], by = c('V1', 'V2'), all.x = T)
    data.table::setorder(rez2, i)
    rez2 <- array(rez2$xx, dim = rep(ncol(input_m), 2))
    rownames(rez2) <- colnames(input_m)
    colnames(rez2) <- colnames(input_m)

    newCombUnique <- combUnique[, -'i']
    newCombUnique <- as.matrix(newCombUnique)
    colnames(newCombUnique) = NULL

    #-------------------------------------------------------------------------------------------
    #average similarity

    YY <- apply(newCombUnique, 1, function(x){
      y <- cbind(input_m[,x[1]],input_m[,x[2]])
      r_mean <- rowMeans(y)
      mean(r_mean)
    })

    data.table::set(combUnique, j = 'yy', value = YY)
    rez3 <- merge(comb2, combUnique[, -'i'], by = c('V1', 'V2'), all.x = T)
    data.table::setorder(rez3, i)
    rez3 <- array(rez3$yy, dim = rep(ncol(input_m), 2))
    rownames(rez3) <- col
    colnames(rez3) <- col
    #--------------------------------------------------------------------------------------------

    min=0;max=1;
    v2=as.vector(rez2);v2=scales::rescale(v2, to=c(min,max))
    v3=as.vector(rez3);v3=scales::rescale(v3, to=c(min,max))

    #Adjusted score of similarity measures
    meas_m=cbind(v2,v3)
    eigs_adj <- apply(meas_m, 1, autri)
    eigs_adj = scales::rescale(eigs_adj, to=c(min,max))

    m4=matrix(eigs_adj,ncol,ncol)
    diag(m4)=1

    rownames(m4) <- col
    colnames(m4) <- col
    m4[is.na(m4)]=0

    if(k==1){
      jac_sim=v2
      prop_sim=v3
      adj_sim=as.vector(m4)
    }else{
      jac_sim=rbind(jac_sim,v2)
      prop_sim=rbind(prop_sim,v3)
      adj_sim=rbind(adj_sim,as.vector(m4))
    }

  }

  rownames(jac_sim)=rownames(prop_sim)=rownames(adj_sim)=rownames(input_m2)
  g_sims_l=list(jac_sim=jac_sim,prop_sim=prop_sim,adj_sim=adj_sim,n_row=length(col),n_col=length(col),genes=rownames(input_m2),names=col)
  return(g_sims_l)
}

#' Compute the Trending Matching similarity
#'
#' Determine which is the similarity between patients using the TM similarity on pathway specific genes
#'
#' @param input_m2 Matrix of profiles (genes x samples)
#' @return Produce the patient similarity network
#' @import data.table
#' @import matrixStats
#' @import scales
#' @import stats
#' @export
#'
Adj_Weighted_Jaccard <- function(input_m) {
  # require(data.table)
  # require(matrixStats)

  input_m <- na.omit(input_m)

  col <- colnames(input_m)
  row <- rownames(input_m)
  nrow <- nrow(input_m)
  ncol <- ncol(input_m)

  samples <- 1:ncol(input_m)
  comb <- data.table::CJ(samples, samples)
  comb[, i := .I]
  comb <- data.table::melt(comb, 'i')
  data.table::setorder(comb, value)
  v2 <- paste0("V", 1:2)
  comb[, variable2 := v2 , keyby = i]
  comb2 <- data.table::dcast(comb, i ~ variable2, value.var = 'value')
  combUnique <- unique(comb2, by = c('V1', 'V2'))   #creation of all combination of genes

  XXcomb=combUnique
  XX <- apply(combUnique[, -'i'], 1, function(x) {
    x2 <- rowRanges(input_m, cols = x)
    s <- colSums2(x2)
    s[1] / s[2]
  })

  data.table::set(XXcomb, j = 'xx', value = XX)
  rez2 <- merge(comb2, XXcomb[, -'i'], by = c('V1', 'V2'), all.x = T)
  data.table::setorder(rez2, i)
  rez2 <- array(rez2$xx, dim = rep(ncol(input_m), 2))
  rownames(rez2) <- colnames(input_m)
  colnames(rez2) <- colnames(input_m)

  newCombUnique <- combUnique[, -'i']
  newCombUnique <- as.matrix(newCombUnique)
  colnames(newCombUnique) = NULL

  #-------------------------------------------------------------------------------------------
  #average similarity

  YY <- apply(newCombUnique, 1, function(x){
    y <- cbind(input_m[,x[1]],input_m[,x[2]])
    r_mean <- rowMeans(y)
    mean(r_mean)
  })

  data.table::set(combUnique, j = 'yy', value = YY)
  rez3 <- merge(comb2, combUnique[, -'i'], by = c('V1', 'V2'), all.x = T)
  data.table::setorder(rez3, i)
  rez3 <- array(rez3$yy, dim = rep(ncol(input_m), 2))
  rownames(rez3) <- col
  colnames(rez3) <- col
  #--------------------------------------------------------------------------------------------

  min=0;max=1;
  v2=as.vector(rez2);v2=scales::rescale(v2, to=c(min,max))
  v3=as.vector(rez3);v3=scales::rescale(v3, to=c(min,max))

  #Adjusted score of similarity measures
  meas_m=cbind(v2,v3)
  eigs_adj <- apply(meas_m, 1, autri)
  eigs_adj = scales::rescale(eigs_adj, to=c(min,max))

  m4=matrix(eigs_adj,ncol,ncol)
  diag(m4)=1

  rownames(m4) <- col
  colnames(m4) <- col
  m4[is.na(m4)]=0
  return(m4)
}

#' Determine the first component of the Trending Matching similarity
#'
#' Given the matrix of patient's profiles, it measures how much the feature information varies between patients.
#' In case the information is similar between two patients, then this similarity is high.
#'
#' @param input_m Matrix of profiles (genes x samples)
#' @param samples_at_columns Boolean. TRUE if samples to compute similarity are at column.
#' @return A similarity matrix aka PSN with the values provided with this first component measure
#' @import data.table
#' @import matrixStats
#' @import stats
#' @export
#'
sim_WJ <- function(input_m,samples_at_columns=TRUE) {
  if(!samples_at_columns){
    input_m=t(input_m)
  }

  samples <- 1:ncol(input_m)
  comb <- data.table::CJ(samples, samples)
  comb[, i := .I]
  comb <- data.table::melt(comb, 'i')
  data.table::setorder(comb, value)
  v2 <- paste0("V", 1:2)
  comb[, variable2 := v2 , keyby = i]
  comb2 <- data.table::dcast(comb, i ~ variable2, value.var = 'value')
  combUnique <- unique(comb2, by = c('V1', 'V2'))

  XX <- apply(combUnique[, -'i'], 1, function(x) {
    x2 <- matrixStats::rowRanges(input_m, cols = x)
    s <- matrixStats::colSums2(x2)
    s[1] / s[2]
  })

  data.table::set(combUnique, j = 'xx', value = XX)
  rez2 <- merge(comb2, combUnique[, -'i'], by = c('V1', 'V2'), all.x = T)
  data.table::setorder(rez2, i)
  rez2 <- array(rez2$xx, dim = rep(ncol(input_m), 2))
  rownames(rez2) <- colnames(input_m)
  colnames(rez2) <- colnames(input_m)
  rez2[is.na(rez2)]=0
  return(rez2)
}

#' Determines how much a patient is connected to its friends against the others
#'
#' Function which computes how much a node is connected to its class, its enemies and then
#' it produces the adjusted eigenscore centrality of the nodes. Higher are both the measures
#' and higher is the score
#'
#' @param m2 Matrix of profiles (genes x samples)
#' @param nAs number of elements of the first group
#' @return A similarity matrix for the column objects
#' @import scales
#' @import Rfast
#' @import stats
#' @export
#'
eig_score_adj = function(m2,nAs){
  seqAs=seq(1,nAs);seqOPP=seq(nAs+1,ncol(m2));
  #Compute how much a node is connected inside its group
  INw=c(Rfast::rowMedians(m2[seqAs,seqAs]),Rfast::rowMedians(m2[seqOPP,seqOPP]))
  #Compute how much a node is connected outside its group
  ABw=c(1-Rfast::rowMedians(m2[seqAs,seqOPP]), 1-Rfast::colMedians(m2[seqAs,seqOPP]))
  #Finalize the vectors with the names of the nodes and create a matrix
  names(ABw)=names(INw)=colnames(m2)
  eigs=cbind(INw,ABw)
  #Rank scale
  eigs=apply(eigs,2,rank)

  #Adjusted score of eigencentrality
  eigs=scaling_matrix(eigs)
  eigs_adj <- apply(eigs, 1, autri)
  eigs_adj = scales::rescale(eigs_adj, to=c(0,1))

  #Update
  eigs=cbind(eigs,eigs_adj)
  return(eigs)
}
