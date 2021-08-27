#' Converts adjacency matrix to edge list
#'
#' Converts a graph represented as adjacency matrix to edge list
#'
#' @param mat adjacency matrix of a graph
#' @param symmetric default:TRUE, the adjacency matrix is symmetric
#' @param diagonal default:FALSE, do not keep the diagonal values in the edge list
#' @param text default:FALSE, do not set the column of the edge list as character vector
#' @return Return edge list of the original network as matrix
#' @export
#'
adj2edg<-function(mat,symmetric=TRUE,diagonal=FALSE,text=FALSE){
  mat<-as.matrix(mat)
  id<-is.na(mat) # used to allow missing
  mat[id]<-"nna"
  if(symmetric){mat[lower.tri(mat)]<-"na"} # use to allow missing values
  if(!diagonal){diag(mat)<-"na"}
  obj<-reshape2::melt(mat)
  colnames(obj)<-c("source","target","value")
  obj<-obj[!obj$value=="na",]
  obj$value[obj$value=="nna"]<-NA
  if(!text){obj$value<-as.numeric(as.character(obj$value))}
  obj=obj[obj$value!=0,]
  return(obj)
}

#' Converts edge list to adjacency matrix
#'
#' Converts a graph represented as edge list to adjacency matrix
#'
#' @param edg1 edge list of a graph
#' @return Return djacency matrix of the original network as matrix
#' @export
#'
edg2adj<-function(edg1){
  mygraph <- igraph::graph.data.frame(edg1,directed = F)
  igraph::get.adjacency(mygraph, sparse = FALSE, attr=colnames(edg1)[3])
}


#' Convert graph list to adjacency matrix
#'
#' Takes a list which contains a vector of values being the
#' vectorized version of an adjacency matrix, the orginal number of columns, rows and
#' the names of the columns. Then, this function convertes the vector to the original adjacency matrix.
#'
#' @param xL list containing vectorized adjacency matrix, its original number of columns,
#' its original number of columns and node names
#' @return Return adjacency matrix of the original network as matrix
#' @export
#'
vec2m=function(xL){
  v=xL$v;nR=xL$n_row;nC=xL$n_col;labelR=xL$row_names;labelC=xL$col_names;
  #Convert vector to matrix with same original dimensions
  m=matrix(v,nR,nC)
  #Apply original labels
  colnames(m)=labelC
  rownames(m)=labelR
  #Return matrix
  return(m)
}
