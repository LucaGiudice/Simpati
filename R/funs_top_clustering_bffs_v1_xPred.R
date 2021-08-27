#' Module of the Best Friend Connector algorithm
#'
#' Undefined
#'
#' @param m Undefined
#' @param testP Undefined
#' @param thW Undefined
#' @param thGroup Undefined
#' @return Undefined
#' @import Rfast
#' @export
#'
strong.comm.xPred2 = function(m,testP,thW=0.4,thGroup=0.2){
  diag(m)=0
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

  #Discard nodes without the root
  keep=sapply(l_BFFs,function(x){
    y=(testP %in% x)
    return(y)
  })
  l_BFFs=c(l_BFFs[testP],l_BFFs[keep])

  if(is.null(l_BFFs[[1]])){
    best_gr=ordR[,testP]
    return(best_gr)
  }

  if(length(l_BFFs)!=0){
    #Join roots that share the 70% of the BFFs
    for(i in 1){
      for(j in (i+1):length(l_BFFs)){
        i_BFFs=unlist(l_BFFs[[i]]);ni=length(i_BFFs);i_root=names(l_BFFs[i]);orig_i_BFFs=i_BFFs
        j_BFFs=unlist(l_BFFs[[j]]);nj=length(j_BFFs);j_root=names(l_BFFs[j])
        intp=(length(intersect(c(i_BFFs,i_root),c(j_BFFs,j_root)))*100)/(min(c(ni,nj))+1)
        if(intp>=70){
          i_BFFs=union(i_BFFs,c(j_BFFs,j_root))
          l_BFFs[[i]]=i_BFFs
          l_BFFs[[j]]=NA
        }
      }
    }
    l_BFFs=l_BFFs[!sapply(sapply(l_BFFs,is.na),function(x){x[1]})]

    #Now creates subnetwork with 60% of nodes and compute how in average are connected
    gr_scores=c(0)
    l_groups=list()
    #Interate over the list of BFFs
    for(k in 1:1){
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
      if(!(testP %in% bffs.final)){
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

        bffs.final=unique(c(testP,y[keep_friends]))
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

#' Module of the Best Friend Connector algorithm
#'
#' Undefined
#'
#' @param m Undefined
#' @param testP Undefined
#' @param thW Undefined
#' @param thGroup Undefined
#' @return Undefined
#' @import Rfast
#' @export
#'
strong.comm.xPred3 = function(m,testP,thW=0.4,thGroup=0.2){
  diag(m)=0
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
  ordR=ordR[1:cutW,]

  best_gr=ordR[1:cutGroup,testP]
  return(best_gr)
}





