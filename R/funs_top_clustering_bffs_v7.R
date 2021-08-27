#' Best Friend Connector algorithm
#'
#' Undefined
#'
#' @param m adjacency matrix of a patient similarity network
#' @param thW Default 0.6, double value between 0 and 1 determining the percentage of best friends to considers for a root (1BFS are the 60 percentage of the patients most similar to the root)
#' @param thGroup Default 0.6, double value between 0 and 1 determining the percentage of best friends to keep in the end as final subclass
#' @return name of the nodes/columns which compose the best group with the most similar best friends and the root
#' @import Rfast
#' @export
#'
strong.comm = function(m,thW=0.6,thGroup=0.6){
  cutW=round(nrow(m)*thW)
  cutGroup=round(nrow(m)*thGroup)
  if(cutW<=2){cutW=3}
  if(cutGroup==1){cutGroup=2}

  #Ranks the connections that each node has with its neighbours
  mR=Rfast::colRanks(m,parallel = FALSE,method="average")
  colnames(mR) = rownames(mR) = colnames(m)
  #Convert the rank values into the name of the neighbours
  ordR=apply(mR,2,function(x){
    names(x)[order(x,decreasing = TRUE)]
  })

  #Cut based on the dimension of cutW, more cutW and more the connections
  #are significantly high, however this can bring to not having a subset of nodes
  #creating a complete network
  ordR=ordR[1:cutW-1,]

  l_BFFs=list()
  #For each node in the network
  for(k in 1:ncol(ordR)){
    #Set the node as root
    subj=colnames(ordR)[k]
    #Get its friends (those with which it has an edge within the 60% of its top interactions)
    friends=ordR[,k]
    #Now iterates over its friends to understand if they are bestfriends with the root
    #the root and a friend are bestfriends if they contain eachother within the 60% of their top interactions
    for(f in 1:length(friends)){
      #Investigate a friend of the root
      friend=friends[f]
      #Get the friends of the root's friend
      col2=ordR[,friend]
      #Check if bff
      bff=subj %in% col2
      #If bff then add the friend's root to the root list of BFFs
      if(bff){
        if(f==1){
          l_BFFs[[subj]]=friend
        }else{
          v=as.character(unlist(l_BFFs[[subj]]))
          l_BFFs[[subj]]=c(v,friend)
        }
      }else{
        next;
      }
    }
  }

  #Discard nodes with few bffs
  l_BFFs=l_BFFs[sapply(l_BFFs,length)>round((cutW-1)/4)]
  if(length(l_BFFs)!=0){
    #Order from the root with more BFFs to the one root with the least of BFFs
    l_BFFs=l_BFFs[order(sapply(l_BFFs,length),decreasing = T)]

    if(F){
      #For each root, measure the coesion of the group defined by its BFFs
      #A bff is integrated in the network of root + bffs if it has at least
      #1/3 of root's bffs as its own bffs
      for(i in 1:length(l_BFFs)){
        coes_bffs=0
        #Extract root and its bffs
        root=names(l_BFFs)[i]
        bffs=unlist(l_BFFs[[i]])
        #Compute how much each bffs has the other root's bffs
        for(j in 1:length(bffs)){
          bff=bffs[j]
          bffs_test=bffs;bffs_test[j]=root
          coes_bffs=c(coes_bffs,length(intersect(unlist(l_BFFs[[bff]]),bffs_test)))
        }
        #Keep only the bffs which have at least 1/3 of the root's bffs as their own bffs
        coes_bffs=coes_bffs[-1]
        keep_bffs=which(coes_bffs>round(length(bffs)/6))
        l_BFFs[[i]]=l_BFFs[[i]][keep_bffs]
      }
      l_BFFs=l_BFFs[sapply(l_BFFs,length)!=0]
    }

    #Join roots that share the 70% of the BFFs
    if(length(l_BFFs)>1){
      for(i in 1:(length(l_BFFs)-1)){
        for(j in (i+1):length(l_BFFs)){
          i_BFFs=unlist(l_BFFs[[i]]);ni=length(i_BFFs);i_root=names(l_BFFs[i]);orig_i_BFFs=i_BFFs
          j_BFFs=unlist(l_BFFs[[j]]);nj=length(j_BFFs);j_root=names(l_BFFs[j])
          intp=(length(intersect(c(i_BFFs,i_root),c(j_BFFs,j_root)))*100)/(nj+1)
          if(intp>=70){
            i_BFFs=union(i_BFFs,c(j_BFFs,j_root))
            l_BFFs[[i]]=i_BFFs
            l_BFFs[[j]]=NA
          }
        }
      }
    l_BFFs=l_BFFs[!sapply(sapply(l_BFFs,is.na),function(x){x[1]})]
    l_BFFs=l_BFFs[!sapply(l_BFFs,is.null)]
    }

    #Now creates subnetwork with 60% of nodes and compute how in average are connected
    gr_scores=c(0)
    l_groups=list()
    #Interate over the list of BFFs
    for(k in 1:length(l_BFFs)){
      #Get the root
      root=names(l_BFFs)[k]
      #Get the BFFs of the root into account
      bffs=bffs.final=unlist(l_BFFs[[k]])

      check_update=1
      #If the group of the ROOT and its BFFs is smaller than the expected resulting group
      #Get the BFFs of the ROOT and add BFFs of the root's BFFs who are BFFs with the root
      while(length(bffs)<(cutGroup-1) & check_update!=0){
        check_update=0
        #Get the BFFs of the root's BFFs who are not directly BFFs of the root
        bffs2=setdiff(unique(unlist(l_BFFs[bffs])),c(root,bffs.final))
        if(length(bffs2)==0){break}
        #Check triangularity
        bffs2=l_BFFs[bffs2]
        check_tri=sapply(bffs2,function(x){
          root %in% x
        })
        #Get the name of the BFFs of the root's BFFs having the root as their own BFF
        bffs2_tri=names(check_tri)[check_tri]
        check_update=sum(check_tri)
        #If exist, then update and iterate again, if none then stop
        if(check_update!=0){
          bffs.final=c(bffs.final,bffs2_tri)
          bffs=bffs2_tri
        }
      }
      #Add to the nodes which will create the network, both the root and its BFFs
      bffs.final=c(root,bffs.final)
      n_bffs.final=length(bffs.final)

      #If the final group is smaller than the expected resulting group then its score is 0
      if(n_bffs.final<cutGroup){
        gr_scores=c(gr_scores,0)
        l_groups[[k]]=NA
        if(k==1){
          gr_scores=gr_scores[-1]
        }
        next;
      }

      #If the final group is larger than the expected resulting group then remove the weakests
      if(n_bffs.final>cutGroup){
        n_del_friends=n_bffs.final-cutGroup
        keep_friends=seq(1,n_bffs.final-n_del_friends)

        x=Rfast::colmeans(m[bffs.final,bffs.final],parallel=F)
        y=bffs.final[order(Rfast::Rank(x,method = "average",descending = TRUE))]

        bffs.final=y[keep_friends]
      }

      final_subm=m[bffs.final,bffs.final]
      final_score=mean(final_subm)
      gr_scores=c(gr_scores,final_score)
      l_groups[[k]]=bffs.final

      if(k==1){
        gr_scores=gr_scores[-1]
      }
    }

    best_gr=l_groups[[which.max(gr_scores)]]
  }else{
    best_gr=NA
  }
  return(best_gr)
}

#' Calls the BFC algorithm over the first patient class
#'
#' Calls the BFC algorithm over the first patient class
#'
#' @param mat adjacency matrix of a patient similarity network
#' @param nAs number of patients composing the first class
#' @param thW Default 0.6, double value between 0 and 1 determining the percentage of best friends to considers for a root (1BFS are the 60 percentage of the patients most similar to the root)
#' @param thGroup Default 0.6, double value between 0 and 1 determining the percentage of best friends to keep in the end as final subclass
#' @param strong Default TRUE, if the first class is signature class of the PSN passed as matrix (mat)
#' @return filtered input PSN without the patients that have not been considered representatives of their own class
#' @export
#'
strong.comm.AA=function(mat,nAs,thW=0.6,thGroup=0.6,strong=TRUE){
  #Get the original mAs
  diag(mat)=0
  seqAs=seq(1,nAs)
  mAs=mat[seqAs,seqAs]

  #Find the stronger/weakest community inside the matrix
  best_As=strong.comm(m=mAs,thW=thW,thGroup=thGroup)
  if(sum(is.na(best_As))>0){return(mat)}
  if(!strong){best_As=setdiff(colnames(mAs),best_As)}
  top_As=match(best_As,colnames(mAs))

  #Reduce the original matrix with the stronger/weakest community
  top_nAs=length(top_As)
  top_mAs=mAs[top_As,top_As]
  mABs=mat[seq(nAs+1,nrow(mat)),c(top_As,seq(nAs+1,ncol(mat)))]
  mBs=mat[top_As,seq(nAs+1,ncol(mat))]
  new_mat=cbind(top_mAs,mBs)
  new_mat=rbind(new_mat,mABs)
  return(new_mat)
}

#' Calls the BFC algorithm over the first patient class correcting the values of the input patient similarity networks
#'
#' Calls the BFC algorithm over the first patient class correcting the values of the input patient similarity networks
#'
#' @param mat adjacency matrix of a patient similarity network
#' @param nAs number of patients composing the first class
#' @param thW Default 0.6, double value between 0 and 1 determining the percentage of best friends to considers for a root (1BFS are the 60 percentage of the patients most similar to the root)
#' @param thGroup Default 0.6, double value between 0 and 1 determining the percentage of best friends to keep in the end as final subclass
#' @param strong Default TRUE, if the first class is signature class of the PSN passed as matrix (mat)
#' @return filtered input PSN without the patients that have not been considered representatives of their own class
#' @export
#'
strong.comm.AA.cor=function(mat,nAs,thW=0.6,thGroup=0.6,strong=T){
  #Get the original mAs
  diag(mat)=0
  seqAs=seq(1,nAs)
  mAs=mat[seqAs,seqAs]
  #Get the corrected one
  mAs_corr=corr_AA_mat(mat,nAs,strong)

  #Find the stronger/weakest community inside the matrix
  best_As=strong.comm(m=mAs_corr,thW=thW,thGroup=thGroup)
  if(sum(is.na(best_As))>0){return(mat)}
  if(!strong){best_As=setdiff(colnames(mAs),best_As)}
  top_As=match(best_As,colnames(mAs))

  #Reduce the original matrix with the stronger/weakest community
  top_nAs=length(top_As)
  top_mAs=mAs[top_As,top_As]
  mABs=mat[seq(nAs+1,nrow(mat)),c(top_As,seq(nAs+1,ncol(mat)))]
  mBs=mat[top_As,seq(nAs+1,ncol(mat))]
  new_mat=cbind(top_mAs,mBs)
  new_mat=rbind(new_mat,mABs)
  return(new_mat)
}

#' Calls the BFC algorithm over the second patient class
#'
#' Calls the BFC algorithm over the second patient class
#'
#' @param mat adjacency matrix of a patient similarity network
#' @param nAs number of patients composing the first class
#' @param thW Default 0.6, double value between 0 and 1 determining the percentage of best friends to considers for a root (1BFS are the 60 percentage of the patients most similar to the root)
#' @param thGroup Default 0.6, double value between 0 and 1 determining the percentage of best friends to keep in the end as final subclass
#' @param strong Default TRUE, if the second class is signature class of the PSN passed as matrix (mat)
#' @return filtered input PSN without the patients that have not been considered representatives of their own class
#' @export
#'
strong.comm.BB=function(mat,nAs,thW=0.6,thGroup=0.6,strong=T){
  #Get the original mBs
  diag(mat)=0
  seqAs=seq(1,nAs)
  seqBs=seq(nAs+1,ncol(mat))
  mBs=mat[seqBs,seqBs]

  #Find the stronger/weakest community inside the matrix
  best_Bs=strong.comm(mBs,thW,thGroup)
  if(sum(is.na(best_Bs))>0){return(mat)}
  if(!strong){best_Bs=setdiff(colnames(mBs),best_Bs)}
  top_Bs=match(best_Bs,colnames(mBs))

  #Reduce the original matrix with the stronger/weakest community
  top_nBs=length(top_Bs)
  top_mBs=mat[c(seqAs,nAs+top_Bs),nAs+top_Bs]
  mABs=mat[nAs+top_Bs,seqAs]
  mAs=mat[seqAs,seqAs]
  new_mat=rbind(mAs,mABs)
  new_mat=cbind(new_mat,top_mBs)
  return(new_mat)
}

#' Calls the BFC algorithm over the second patient class correcting the values of the input patient similarity networks
#'
#' Calls the BFC algorithm over the second patient class correcting the values of the input patient similarity networks
#'
#' @param mat adjacency matrix of a patient similarity network
#' @param nAs number of patients composing the first class
#' @param thW Default 0.6, double value between 0 and 1 determining the percentage of best friends to considers for a root (1BFS are the 60 percentage of the patients most similar to the root)
#' @param thGroup Default 0.6, double value between 0 and 1 determining the percentage of best friends to keep in the end as final subclass
#' @param strong Default TRUE, if the second class is signature class of the PSN passed as matrix (mat)
#' @return filtered input PSN without the patients that have not been considered representatives of their own class
#' @export
#'
strong.comm.BB.cor=function(mat,nAs,thW=0.6,thGroup=0.6,strong=T){
  #Get the original mBs
  diag(mat)=0
  seqAs=seq(1,nAs)
  seqBs=seq(nAs+1,ncol(mat))
  mBs=mat[seqBs,seqBs]
  #Get the corrected mBs
  mBs_corr=corr_BB_mat(mat,nAs,strong)

  #Find the stronger/weakest community inside the matrix
  best_Bs=strong.comm(mBs_corr,thW,thGroup)
  if(sum(is.na(best_Bs))>0){return(mat)}
  if(!strong){best_Bs=setdiff(colnames(mBs),best_Bs)}
  top_Bs=match(best_Bs,colnames(mBs))

  #Reduce the original matrix with the stronger/weakest community
  top_nBs=length(top_Bs)
  top_mBs=mat[c(seqAs,nAs+top_Bs),nAs+top_Bs]
  mABs=mat[nAs+top_Bs,seqAs]
  mAs=mat[seqAs,seqAs]
  new_mat=rbind(mAs,mABs)
  new_mat=cbind(new_mat,top_mBs)
  return(new_mat)
}

#' Euclidean distance between two points
#'
#' Euclidean distance between two points
#'
#' @param P1 vector with two values defining the x and y coordinates of the first point
#' @param P2 vector with two values defining the x and y coordinates of the second point
#' @return Euclidean distance between two points
#' @export
#'
twoDist=function(P1,P2){
  x1=P1[1];y1=P1[2]
  x2=P2[1];y2=P2[2]
  dist=sqrt((x2-x1)^2 + (y2-y1)^2)
  return(dist)
}

#' Quadrilateral area defined by four points
#'
#' Quadrilateral area defined by four points
#'
#' @param A vector with two values defining the x and y coordinates of the first point
#' @param B vector with two values defining the x and y coordinates of the second point
#' @param C vector with two values defining the x and y coordinates of the third point
#' @param D vector with two values defining the x and y coordinates of the fourth point
#' @return Quadrilateral area defined by four points
#' @export
#'
area_quadr=function(A,B,C,D){
  a=twoDist(A,B);
  b=twoDist(B,C);
  c=twoDist(C,D);
  d=twoDist(D,A);
  e=twoDist(A,C);
  f=twoDist(B,D)
  Area=(sqrt((4*(e^2)*(f^2))-(((b^2)+(d^2)-(a^2)-(c^2))^2)))/4
  return(Area)
}

#' Corrects the values of similarities of patients belonging to the first class
#'
#' Given the intra and inter similarities of one patient belonging to the first class which is signature,
#' it increases the intra values more the inter are low, it decreases the intra values more the inter are high.
#' In case, the first class is not signature it corrects in the opposite way.
#'
#' @param mat adjacency matrix of a patient similarity network
#' @param nAs number of patients composing the first class
#' @param strong Default TRUE, if the first class is signature class of the PSN passed as matrix (mat)
#' @return correct input PSN with adjusted values of similarities
#' @export
#'
corr_AA_mat = function(mat,nAs,strong=T){
  diag(mat)=0
  #Find the position of As and Bs
  seqAs=seq(1,nAs)
  seqBs=seq(nAs+1,ncol(mat))

  #Compute how in average each A is connected with the Bs
  wA_cor=seqAs
  for(rowA in seqAs){
    wA_cor[rowA]=mean(mat[rowA,seqBs])
  }

  matAA=mat[seqAs,seqAs]
  #Correct the weigth A1-A2 with how weaker is their connection with Bs
  for(rowA in seqAs){
    #Extract the original weight of the A1
    wAA=matAA[rowA,seqAs];wAA[wAA==0]=0.0001;
    #Extract how A1 is connected with the Bs
    w_cor1=wA_cor[rowA]
    #Extract the other As to consider for the correction
    seqA2s=setdiff(seqAs,rowA)
    for(rowA2 in seqA2s){
      #Extract how A2 is connected with the Bs
      w_cor2=wA_cor[rowA2]

      #Now apply the correction to the weight A1-A2
      w_raw=wAA[rowA2]
      if(strong){
        neg_w_cor1=1-w_cor1
        neg_w_cor2=1-w_cor2
      }else{
        neg_w_cor1=w_cor1
        neg_w_cor2=w_cor2
      }

      #Set coordinates
      A=c(0,0);B=c(0,w_raw);C=c(0,neg_w_cor1);D=c(w_raw,neg_w_cor2)
      #Compute area quadrilat.
      wAA[rowA2]=area_quadr(A,B,C,D)
    }
    matAA[rowA,seqAs]=wAA
  }
  matAA[lower.tri(matAA)]=matAA[upper.tri(matAA)]
  return(matAA)
}

#' Corrects the values of similarities of patients belonging to the second class
#'
#' Given the intra and inter similarities of one patient belonging to the second class which is signature,
#' it increases the intra values more the inter are low, it decreases the intra values more the inter are high.
#' In case, the first class is not signature it corrects in the opposite way.
#'
#' @param mat adjacency matrix of a patient similarity network
#' @param nAs number of patients composing the second class
#' @param strong Default TRUE, if the second class is signature class of the PSN passed as matrix (mat)
#' @return correct input PSN with adjusted values of similarities
#' @export
#'
corr_BB_mat = function(mat,nAs,strong=T){
  diag(mat)=0
  #Find the position of As and Bs
  seqAs=seq(1,nAs)
  seqBs=seq(nAs+1,ncol(mat))

  #Compute how much a node B is connected to the As
  wB_cor=0
  for(rowB in seqBs){
    wB_cor=c(wB_cor,mean(mat[rowB,seqAs]))
  }
  wB_cor=wB_cor[-1]

  #Get matrix of Bs because we are going to adjust their weights
  matBB=mat[seqBs,seqBs];seqBs=seq(1,ncol(matBB))
  #Correct the weigth B1-B2 with how weak/strong is their connection with As
  #For each row of Bs
  for(rowB in seqBs){
    #Extract the original weight of the B1
    wBB=matBB[rowB,seqBs];wBB[wBB==0]=0.0001;
    #Extract how B1 is connected with the As
    w_cor1=wB_cor[rowB]
    #Skip the correction of the running B1 itself
    seqB2s=setdiff(seqBs,rowB)
    #For each of the other Bs
    for(rowB2 in seqB2s){
      #Extract how B2 is connected with the As
      w_cor2=wB_cor[rowB2]

      #Now apply the correction to the weight B1-B2
      w_raw=wBB[rowB2]
      if(strong){
        #In case you want to correct on how much the
        #B1-As and B2-As are low for detecting the
        #strongest Bs
        neg_w_cor1=1-w_cor1
        neg_w_cor2=1-w_cor2
      }else{
        #In case you want to correct on how much the
        #B1-As and B2-As are high for detecting the
        #weakest Bs
        neg_w_cor1=w_cor1
        neg_w_cor2=w_cor2
      }

      #Set coordinates
      A=c(0,0);B=c(0,w_raw);C=c(0,neg_w_cor1);D=c(w_raw,neg_w_cor2)
      #Compute area quadrilat.
      wBB[rowB2]=area_quadr(A,B,C,D)
    }
    #Replace
    matBB[rowB,seqBs]=wBB
  }
  matBB[lower.tri(matBB)]=matBB[upper.tri(matBB)]
  return(matBB)
}

#' Wrapper of the BFC algorithm to plot a pathway specific patient similarity network
#'
#' Given a patient similarity network, for each patient class, it iterates the BFC algorithm in order to
#' create small subclasses of similar patients.
#'
#' @param mat adjacency matrix of a patient similarity network
#' @param nAs number of patients composing the second class
#' @param thGroup Default 0.2, double value between 0 and 1 determining the percentage of best friends to keep in the end as final subclass
#' @param thW Default 0.2, double value between 0 and 1 determining the percentage of best friends to considers for a root (1BFS are the 60 percentage of the patients most similar to the root)
#' @param A_bol Default TRUE, if the first class is signature class of the PSN passed as matrix (mat)
#' @return The nodes of each subclass that must be kept in the plot and the nodes that must be summerized in the legend
#' @import matrixStats
#' @export
#'
reduce_net = function(mat,nAs,thGroup=0.2,thW=0.2,A_bol=T){
  #Compute how much the patient's profiles of interactions are similar
  mat=sim_WJ(mat)
  #Temp vars
  As_gr_l=list();top_A_gr="start";k=0;
  #Get the matrix of selected group of patients
  if(A_bol==T){
    cat("\n I work on the first class of patients \n")
    seqAs=seq(1,nAs);
    seqOPP=seq(nAs+1,ncol(mat))
  }else{
    cat("\n I work on the second class of patients \n")
    seqAs=seq(nAs+1,ncol(mat))
    seqOPP=seq(1,nAs);
  }
  mAs=mat[seqAs,seqAs];diag(mAs)=0;
  #dimension of the first subgroup which be summarized
  cutGroup=round(nrow(mAs)*thGroup)
  #Starting colnames
  rem_pats=colnames(mAs)
  #Iterate until the final matrix had a dimension equal to 10/11
  while(cutGroup>2){
    k=k+1
    cat(k,"with",thW,"equal to",cutGroup,"\n")
    #Compute the group of patients in the class who have the most similar profiles of interactions
    As_gr=strong.comm(mAs,thW,thGroup)
    #Save the group
    As_gr_l[[k]]=As_gr
    #Find their indexes
    As_gr_del=match(As_gr,colnames(mAs))
    #Remove them
    mAs=mAs[-As_gr_del,-As_gr_del]
    #Detect the size of the next group
    cutGroup=round(nrow(mAs)*thGroup)

    #Now we have to find the patient who will rapresent the group
    #Extract how much these patients have interactions similar to the other ones
    net_gr=mat[,As_gr]
    #Compute how much they are similar between eachother
    net_gr1=sim_WJ(net_gr)
    #If the patients of the group have very similar profile then
    if(mean(unlist(net_gr1[upper.tri(net_gr1)]))>=0.7){
      #Find the patient which in the group is the most dissimilar from being an enemy
      net_gr=mat[seqOPP,As_gr]
      #Update the final vector of patients
      top_A_gr=c(top_A_gr,colnames(net_gr)[which.min(matrixStats::colMeans2(net_gr))])
    }else{
      #Find the patient which is most similar to the other ones in the group and
      #Update the final vector of patients
      top_A_gr=c(top_A_gr,colnames(net_gr1)[which.max(matrixStats::colMeans2(net_gr1))])
    }

    #Detect the patients which are remaining ungrouped
    rem_pats=setdiff(rem_pats,As_gr)
  }
  #Update
  top_A_gr=top_A_gr[-1]
  keep_A_patients=c(top_A_gr,rem_pats)
  names(As_gr_l)=top_A_gr
  #Return
  res_l=list(top_pats_l=As_gr_l,keep_pats=keep_A_patients)
  return(res_l)
}
