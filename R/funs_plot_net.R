#' This function creates an image that visualize a Patient Similarity Network (PSN)
#' @param m a similarity matrix
#' @param label1 a string indicating the name of the first class
#' @param label2 a string indicating the name of the second class
#' @param top_connections an integer in [1,100], percentage of top connection to visualize
#' @param node_highlight a string or vector of strinvs indicating the node id(s) to highlight in the network. By default
#' all connections in top_connection are highlited 
#' @param image_type a string indicating the format of the output image
#' @param image_name a string containing the name of the output image
#' @return a list representing a qgraph object. The method saves an image with the visualization
#' of the network 
#' @import data.table
#' @import matrixStats
#' @import scales
#' @importFrom reshape2 melt
#' @import plyr
#' @import ggplot2
#' @importFrom stats quantile
#' @export
#'
plot_network = function(m,top_connections=30,node_highlight="all", image_type="png",image_name="plot_network"){
  label1=names(table(substr(colnames(m),1,1))[1])
  label2=names(table(substr(colnames(m),1,1))[2])
  same_labels = length(rownames(m)) == length(colnames(m))
  same_labels = same_labels && (sum(is.na(match(rownames(m),colnames(m)))) == 0)
  same_labels = same_labels && (sum(is.na(match(colnames(m),rownames(m)))) == 0)
  if(!same_labels){
    cat("Labels on rows doesn't match the labels on columns. Plot aborted")
    return(NA)
  }
  if(sum(duplicated(rownames(m)))){
    cat("Duplicated labels, plot aborted. Call remove_duplicates(PSN) to delete the duplicates.")
    return(NA)
  }
  if(length(unique(substring(rownames(m), 1, 1))) == 1){
    cat("Labels of the two class must start with different letters. Plot aborted")
    return(NA)
  }
  first_class_label = substring(rownames(m)[1], 1, 1)
  nAs = sum(startsWith(rownames(m),first_class_label))
  #Set the same length for both the labels
  end_l1=nchar(label1)
  end_l2=nchar(label2)
  if(end_l2<end_l1){end_l1=end_l2}
  label1=substr(label1,1,end_l1)
  label2=substr(label2,1,end_l1)
  
  #Find the nodes which have very similar interactions to all the other nodes
  #Creating groups of nodes, everytime reducing the starting set of nodes of the 20%
  #Finally, represent each group by only one node to reduce the number of future plotted nodes
  Asel=reduce_net(m,nAs=nAs,A_bol=TRUE)
  Akeep=Asel$keep_pats
  Al=Asel[[1]]
  Bsel=reduce_net(m,nAs=nAs,A_bol=FALSE)
  Bkeep=Bsel$keep_pats
  Bl=Bsel[[1]]

  #Apply reduction
  m2=m[c(Akeep,Bkeep),c(Akeep,Bkeep)]
  nAs=length(Akeep)
  memb=c(rep(label1,nAs),rep(label2,length(Bkeep)))
  #Compute adjusted eigen centrality for each node which will be rapresented
  #as size of the node in the plot
  eigs=eig_score_adj(m2,nAs)
  
  #Rescaled score for matching with node size
  eigs_adj_sc=scales::rescale(eigs[,"eigs_adj"], to=c(2,5))
  
  #Produce edge list
  el=adj2edg(m2)

  #Set the default colors for the nodes and edges due to the groups
  edge_cols=c("#CC0033","#948df7","#0066FF")
  node_cols=c("#FF3300","#00FF66","#00CCFF")
  
  #Set the colors for the nodes
  node_col=memb
  node_col[node_col==label1]=node_cols[1]
  node_col[node_col==label2]=node_cols[3]
  
  #Create edge list annotation converting edges to values and setting their colors
  el_ann=cbind(substr(el[,1],1,end_l1),substr(el[,2],1,end_l1))
  el_ann[el_ann==label1]=0
  el_ann[el_ann==label2]=1
  el_ann=as.data.frame(as.matrix(el_ann))
  el_ann$V1=as.numeric(as.character(el_ann$V1))
  el_ann$V2=as.numeric(as.character(el_ann$V2))
  rs=as.character(rowSums(el_ann))
  rs[rs=="0"]=edge_cols[1]
  rs[rs=="1"]=edge_cols[2]
  rs[rs=="2"]=edge_cols[3]
  el_ann$col_edges=rs
  
  #Prepare node annotation s.t. if a node is representing a group then it has it
  #in the annotation and future legend
  node_ann = "start"
  node_names = colnames(m2)
  for(k in 1:length(node_names)){
    node_name=node_names[k]
    A_text=paste(Al[[node_name]],collapse = ",")
    B_text=paste(Bl[[node_name]],collapse = ",")
    if(A_text!=""){
      node_ann=c(node_ann,A_text)
      next
    }
    if(B_text!=""){
      node_ann=c(node_ann,B_text)
      next
    }
    
    node_ann=c(node_ann,node_name)
  }
  node_ann=node_ann[-1]
  names(node_ann)=node_names
  
  #Network plot
  cat("Plotting \n")
  if(top_connections < 0 || top_connections > 100){
    cat("top_connections must be between 0 and 100. Plot aborted.\n")
    return(NA)
  }
  top_connections = 100 - top_connections
  top_connections=top_connections%/%10
  qs=quantile(el$value,probs = seq(0.1,1,0.10))
  if(length(node_highlight) > 1 || node_highlight != "all"){
    node_highlight_checked = c()
    for(i in 1:length(node_highlight)){
      if(node_highlight[i] %in% names(eigs_adj_sc)){
        node_highlight_checked[length(node_highlight_checked)+1] = node_highlight[i]
      }
      else{
        cat("Node to highlight with id ",node_highlight[i],
          " not found in the network.\n",sep="")
      }
    }
    for(j in 1:(dim(el)[1])){
      source = as.character(el[j,1][[1]])
      target = as.character(el[j,2][[1]])
      if(!((source %in% node_highlight_checked) || (target %in% node_highlight_checked))){
        rs[j] = "#d2d2d2"
      }
    }
  }
  threshold = if(top_connections>0) qs[top_connections] else 0
  Q <- qgraph::qgraph(el, threshold = threshold, edgelist=TRUE, weighted=TRUE, directed = FALSE,
              borders = FALSE, groups=memb, vsize = eigs_adj_sc, esize = 5,
              edge.color=rs, color=unique(node_col), label.cex=1.2,
              repulsion=0.7, layout="spring",
              legend = TRUE, legend.mode="style1", GLratio=1.5, nodeNames=node_ann,
              filetype=image_type,filename=image_name,width=14,height=12)
}