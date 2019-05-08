IndivGraph=function(Data, metric, q, status, graph=TRUE, test=FALSE){
  pal=c(brewer.pal(12,name="Set3"), brewer.pal(n=8, name="Set2"), brewer.pal(n=4, name="Set1"))
  
  n=sum(Data$Condition==status)
  p=ncol(metric)
  
  MyEdgelist_raw=data.frame(node1=row(metric)[upper.tri(metric)], node2=col(metric)[upper.tri(metric)],
                            cor=metric[upper.tri(metric)], stringsAsFactors = FALSE)
  MyEdgelist=MyEdgelist_raw[sort.list(abs(MyEdgelist_raw$cor), decreasing = TRUE),]
  
  cor_q=abs(MyEdgelist$cor[q])
  
  adjacency=matrix(0,p,p)
  adjacency[upper.tri(adjacency)]=ifelse(abs(metric[upper.tri(metric)])>=cor_q, yes=1, no=0)
  adjacency=adjacency+t(adjacency)
  colnames(adjacency)=rownames(adjacency)=colnames(Data$Data)
  
  if (test){
    z_list=zTransform(MyEdgelist$cor)
    z_pval=sapply(1:length(z_list), FUN=function(i){Zpval(z_list[i], n)})
    MyEdgelist=cbind(MyEdgelist, z_list, z_pval)
    colnames(MyEdgelist)[4:5]=c('z_cor', 'z_pvalue')
  }
  
  z_q=zTransform(cor_q)
  pval=Zpval(z_q, n)
  
  if (graph){
    gg=GetGraph(Data=Data, adjacency=adjacency)
    res=list(Graph=gg, adjacency=adjacency, Edgelist=MyEdgelist, pval=pval)
  } else {
    res=list(adjacency=adjacency, Edgelist=MyEdgelist, pval=pval)
  }
  return(res)
}


DiffGraph=function(Data, metric, q, graph=TRUE, test=FALSE){
  pal=c(brewer.pal(12,name="Set3"), brewer.pal(n=8, name="Set2"), brewer.pal(n=4, name="Set1"))
  
  n=nrow(Data$Data0)
  
  MyEdgelist_raw=data.frame(node1=row(metric)[upper.tri(metric)], node2=col(metric)[upper.tri(metric)],
                            metric=metric[upper.tri(metric)], stringsAsFactors = FALSE)
  MyEdgelist=MyEdgelist_raw[sort.list(abs(MyEdgelist_raw$metric), decreasing = TRUE),]
  
  cor_q=abs(MyEdgelist$metric[q])
  
  adjacency=matrix(0,p,p)
  adjacency[upper.tri(adjacency)]=ifelse(abs(metric[upper.tri(metric)])>=cor_q, yes=1, no=0)
  adjacency=adjacency+t(adjacency)
  colnames(adjacency)=rownames(adjacency)=colnames(Data$Data0)
  
  if (test){
    chi_pval=sapply(MyEdgelist[,3], FUN=function(x){
      pval=1-pchisq(n*x, df=1)
      pval=ifelse(pval==0,
                  yes=-pchisq(n*x, df=1, log.p = TRUE), no=pval)
      return(pval)})
    MyEdgelist=cbind(MyEdgelist, chi_pval)
    colnames(MyEdgelist)[4]='chi_pval'
  }
  
  z_q=zTransform(cor_q)
  pval=Zpval(z_q, n)
  
  if (graph){
    gg=GetGraph(Data=Data, adjacency=adjacency)
    res=list(Graph=gg, adjacency=adjacency, Edgelist=MyEdgelist, pval=pval)
  } else {
    res=list(adjacency=adjacency, Edgelist=MyEdgelist, pval=pval)
  }
  return(res)
}


SelectionMetric=function(metric, q, symmetric=TRUE){
  if (symmetric){  
    MyEdgelist_raw=data.frame(node1=row(metric)[upper.tri(metric)], node2=col(metric)[upper.tri(metric)],
                              cor=metric[upper.tri(metric)], stringsAsFactors = FALSE)
  } else {
    MyEdgelist_raw=data.frame(node1=as.vector(row(metric)), node2=as.vector(col(metric)), 
                              cor=as.vector(metric), stringsAsFactors = FALSE)
  }
  MyEdgelist=MyEdgelist_raw[sort.list(abs(MyEdgelist_raw$cor), decreasing = TRUE),]
  cor_q=abs(MyEdgelist$cor[q])
  if (symmetric){
    adjacency=matrix(0,nrow(metric),ncol(metric))
    adjacency[upper.tri(adjacency)]=ifelse(abs(metric[upper.tri(metric)])>=cor_q, yes=1, no=0)
    adjacency=adjacency+t(adjacency)
  } else {
    adjacency=ifelse(abs(metric)>=cor_q, yes=1, no=0)
  }
  return(adjacency)
}


# GetGraph=function(Data, adjacency, weighted=NULL){
#   pal=c(brewer.pal(12,name="Set3"), brewer.pal(n=8, name="Set2"), brewer.pal(n=4, name="Set1"))
#   
#   gg=graph_from_adjacency_matrix(adjacency, mode = 'undirected', weighted = weighted)
#   V(gg)$size=as.numeric(as.character(cut(degree(gg), breaks = 3, labels = 1:3)))*3
#   V(gg)$color=pal[as.numeric(Data$Annot$Vertex_color[V(gg)$name])]
#   V(gg)$frame.color=V(gg)$color
#   V(gg)$label.family="sans"
#   E(gg)$color="black"
#   V(gg)$label.cex=as.numeric(as.character(cut(degree(gg), breaks = 3, labels = c(0.4, 0.6, 0.7))))
#   V(gg)$label.color="grey30"
#   E(gg)$width=0.5
#   #  if (!is.null(weighted)){
#   #  E(gg)$color=c('red', 'blue', 'forestgreen')
#   #}
#   
#   return(gg)
# }


GetGraph=function(Data, adjacency, weighted=NULL, pal=NULL){
  if (is.null(pal)){
    pal=c(brewer.pal(12,name="Set3"), brewer.pal(n=8, name="Set2"), brewer.pal(n=4, name="Set1"))
  }
  
  gg=graph_from_adjacency_matrix(adjacency, mode = 'undirected', weighted = weighted)
  V(gg)$size=as.numeric(as.character(cut(degree(gg), breaks = 3, labels = 1:3)))*3
  V(gg)$color=pal[as.numeric(Data$Annot$Vertex_color[V(gg)$name])]
  V(gg)$frame.color=V(gg)$color
  V(gg)$label.family="sans"
  E(gg)$color="black"
  V(gg)$label.cex=as.numeric(as.character(cut(degree(gg), breaks = 3, labels = c(0.4, 0.6, 0.7))))
  V(gg)$label.color="grey30"
  E(gg)$width=0.5
  if (!is.null(weighted)){
    E(gg)$color=c('red', 'blue', 'forestgreen')[E(gg)$weight]
  }
  
  return(gg)
}


zTransform=function(r){
  return(1/2*log((1+r)/(1-r)))
}


Zpval=function(z, n){
  z=abs(z) # true because a correlation is scaled
  pval=(1-pnorm(z, mean = 0, sd=1/sqrt(n-3), log.p = FALSE))/2
  if(pval==0|pval==1){
    pval=-pnorm(z, mean = 0, sd=1/sqrt(n-3), log.p = TRUE)*2 # ln(1-x) = -x
  }
  return(pval)
}


# StabilityIndiv=function(Data, status=0, K, Q, no_cores=NULL){
#   p=ncol(Data$Data)
#   
#   no_cores <- ifelse(is.null(no_cores), yes=detectCores(), no=no_cores)
#   clust <- makeCluster(no_cores, type="FORK")
#   
#   Output=parSapply(clust, Q, FUN=function(q){
#     Cum_adj=matrix(0,p,p)
#     for (k in 1:K){
#       s=sample(1:nrow(Data$Data), size = 0.8*nrow(Data$Data))
#       MyDataIter=Data
#       MyDataIter$Data=MyDataIter$Data[s,]
#       MyDataIter$Condition=MyDataIter$Condition[s]
#       if (nrow(MyDataIter$Data)>ncol(MyDataIter$Data)){ # if n_sub>p
#         pcor1=pcor(MyDataIter$Data)$estimate
#       } else {
#         pcor1=pcor.shrink(MyDataIter$Data)
#       }
#       gg_iter=IndivGraph(Data=MyDataIter, metric=pcor1, q=q, status=status, graph=FALSE, test=FALSE)
#       Cum_adj=Cum_adj+gg_iter$adjacency
#     }
#     Cum_adj=Cum_adj/K
#     return(Cum_adj)
#   })
#   
#   stopCluster(clust)
#   
#   Output=array(Output, dim = c(p,p,length(Q)))
#   return(list(Array=Output, Q=Q, K=K))
# }


StabilityIN=function(Data, status=0, K, Q, no_cores=NULL, sub_size=0.8, type="conditional"){
  p=ncol(Data$Data)
  
  no_cores <- ifelse(is.null(no_cores), yes=detectCores(), no=no_cores)
  clust <- makeCluster(no_cores, type="FORK")
  
  Output=parSapply(clust, Q, FUN=function(q){
    Cum_adj=matrix(0,p,p)
    for (k in 1:K){
      s=sample(1:nrow(Data$Data), size = sub_size*nrow(Data$Data))
      MyDataIter=Data
      MyDataIter$Data=MyDataIter$Data[s,]
      MyDataIter$Condition=MyDataIter$Condition[s]
      if (type=="conditional"){
        if (nrow(MyDataIter$Data)>ncol(MyDataIter$Data)){ # if n_sub>p
          pcor1=pcor(MyDataIter$Data)$estimate
        } else {
          pcor1=pcor.shrink(MyDataIter$Data)
        }
        Cum_adj=Cum_adj+SelectionMetric(metric=pcor1, q=q, symmetric=TRUE)
      } else {
        cor1=cor(MyDataIter$Data)
        Cum_adj=Cum_adj+SelectionMetric(metric=cor1, q=q, symmetric=TRUE)
      }
    }
    Cum_adj=Cum_adj/K
    return(Cum_adj)
  })
  
  stopCluster(clust)
  
  Output=array(Output, dim = c(p,p,length(Q)))
  return(list(Array=Output, Q=Q, K=K))
}


StabilityDN=function(Data, K, Q, no_cores=NULL, sub_size=0.8, type="conditional"){
  p=ncol(Data$Data0)
  
  no_cores <- ifelse(is.null(no_cores), yes=detectCores(), no=no_cores)
  clust <- makeCluster(no_cores, type="FORK")
  
  Output=parSapply(clust, Q, FUN=function(q){
    Cum_adj=matrix(0,p,p)
    for (k in 1:K){
      s0=sample(1:nrow(Data$Data0), size = sub_size*nrow(Data$Data0))
      s1=sample(1:nrow(Data$Data1), size = sub_size*nrow(Data$Data1))
      MyDataIter=Data
      MyDataIter$Data0=MyDataIter$Data0[s0,]
      MyDataIter$Data1=MyDataIter$Data1[s1,]
      if (type=="conditional"){
        hat_omega0=solve(cov(MyDataIter$Data0))
        hat_omega1=solve(cov(MyDataIter$Data1))
        Delta=ComputeDelta(omegaX0=hat_omega0, omegaX1=hat_omega1)
      }
      if (type=="marginal"){
        hat_omega0=cor(MyDataIter$Data0)
        hat_omega1=cor(MyDataIter$Data1)
        Delta=hat_omega1-hat_omega0
      }
      Cum_adj=Cum_adj+SelectionMetric(metric=Delta, q=q, symmetric=TRUE)
    }
    Cum_adj=Cum_adj/K
    return(Cum_adj)
  })
  
  stopCluster(clust)
  
  Output=array(Output, dim = c(p,p,length(Q)))
  return(list(Array=Output, Q=Q, K=K))
}


ComputeN=function(p1, p2=NULL, symmetric=FALSE){
  if (symmetric){
    return(p1*(p1-1)/2)
  } else {
    return(p1*p2)
  }
}


StabilityMOIN=function(Data, status=0, K, Q_card, no_cores=NULL, sub_size=0.8, type="conditional"){ 
  p=ncol(Data$Data)
  
  L=length(unique(Data$Platform))
  l=matrix(c(1,4,3,2),ncol=2) # to generalise
  symm=matrix(c(TRUE, FALSE, FALSE, TRUE), ncol=2) # to generalise
  pk=table(Data$Platform)
  p1=matrix(rep(pk, 2), ncol=2) # to generalise
  p2=t(p1) # to generalise
  Nblock=sapply(1:L, FUN=function(i){sapply(1:L, FUN=function(j){ComputeN(p1=p1[i,j], p2=p2[i,j], symmetric=symm[i,j])})})
  # Nblock=c(sapply(pk, FUN=function(p1){ComputeN(p1)}), ComputeN(pk[1], pk[2])) # only for 2 blocks
  
  Q=lapply(1:(L+L*(L-1)/2), FUN=function(k){round(seq(1, Nblock[l==k], length.out = Q_card))})
  
  
  no_cores <- ifelse(is.null(no_cores), yes=detectCores(), no=no_cores)
  clust <- makeCluster(no_cores, type="FORK")
  
  Output=parSapply(clust, 1:Q_card, FUN=function(l){
    Cum_adj1=matrix(0,pk[1],pk[1])
    Cum_adj2=matrix(0,pk[2],pk[2])
    Cum_adjb=matrix(0,pk[1],pk[2])
    for (k in 1:K){
      s=sample(1:nrow(Data$Data), size = sub_size*nrow(Data$Data))
      MyDataIter=Data
      MyDataIter$Data=MyDataIter$Data[s,]
      MyDataIter$Condition=MyDataIter$Condition[s]
      if (type=="conditional"){
        if (nrow(MyDataIter$Data)>ncol(MyDataIter$Data)){ # if n_sub>p
          pcor1=pcor(MyDataIter$Data)$estimate
        } else {
          pcor1=pcor.shrink(MyDataIter$Data)
        }
      }
      else {
        pcor1=cor(MyDataIter$Data)
      }
      Cum_adj1=Cum_adj1+SelectionMetric(metric=pcor1[1:pk[1], 1:pk[1]], q=(Q[[1]])[l], symmetric=TRUE)
      Cum_adj2=Cum_adj2+SelectionMetric(metric=pcor1[(pk[1]+1):(pk[1]+pk[2]), (pk[1]+1):(pk[1]+pk[2])], q=(Q[[2]])[l], symmetric=TRUE)
      Cum_adjb=Cum_adjb+SelectionMetric(metric=pcor1[1:pk[1], (pk[1]+1):(pk[1]+pk[2])], q=(Q[[3]])[l], symmetric=FALSE)
    }
    Cum_adj=rbind(cbind(Cum_adj1,Cum_adjb), cbind(t(Cum_adjb), Cum_adj2))/K
    return(Cum_adj)
  })
  
  stopCluster(clust)
  
  Output=array(Output, dim = c(p,p,Q_card))
  return(list(Array=Output, Q=Q, K=K))
}


StabilityMODN=function(Data, K, Q_card, no_cores=NULL, sub_size=0.8, type="conditional"){ 
  p=ncol(Data$Data0)
  
  L=length(unique(Data$Platform))
  l=matrix(c(1,4,3,2),ncol=2) # to generalise
  symm=matrix(c(TRUE, FALSE, FALSE, TRUE), ncol=2) # to generalise
  pk=table(Data$Platform)
  p1=matrix(rep(pk, 2), ncol=2) # to generalise
  p2=t(p1) # to generalise
  Nblock=sapply(1:L, FUN=function(i){sapply(1:L, FUN=function(j){ComputeN(p1=p1[i,j], p2=p2[i,j], symmetric=symm[i,j])})})
  # Nblock=c(sapply(pk, FUN=function(p1){ComputeN(p1)}), ComputeN(pk[1], pk[2])) # only for 2 blocks
  
  Q=lapply(1:(L+L*(L-1)/2), FUN=function(k){round(seq(1, Nblock[l==k], length.out = Q_card))})
  
  
  no_cores <- ifelse(is.null(no_cores), yes=detectCores(), no=no_cores)
  clust <- makeCluster(no_cores, type="FORK")
  
  Output=parSapply(clust, 1:Q_card, FUN=function(l){
    Cum_adj1=matrix(0,pk[1],pk[1])
    Cum_adj2=matrix(0,pk[2],pk[2])
    Cum_adjb=matrix(0,pk[1],pk[2])
    for (k in 1:K){
      s0=sample(1:nrow(Data$Data0), size = sub_size*nrow(Data$Data0))
      s1=sample(1:nrow(Data$Data1), size = sub_size*nrow(Data$Data1))
      MyDataIter=Data
      MyDataIter$Data0=MyDataIter$Data0[s0,]
      MyDataIter$Data0=MyDataIter$Data1[s1,]
      if (type=="conditional"){
        hat_omega0=solve(cov(MyDataIter$Data0))
        hat_omega1=solve(cov(MyDataIter$Data1))
        Delta=ComputeDelta(omegaX0=hat_omega0, omegaX1=hat_omega1)
      } else {
        hat_omega0=cor(MyDataIter$Data0)
        hat_omega1=cor(MyDataIter$Data1)
        Delta=hat_omega1-hat_omega0
      }
      # SubDelta=NULL
      # for (sub_k in 1:10){
      #   s0=sample(1:nrow(MyDataIter$Data0), size = sub_size*nrow(MyDataIter$Data0))
      #   s1=sample(1:nrow(MyDataIter$Data1), size = sub_size*nrow(MyDataIter$Data1))
      #   MyDataSubIter=MyDataIter
      #   MyDataSubIter$Data0=MyDataSubIter$Data0[s0,]
      #   MyDataSubIter$Data0=MyDataSubIter$Data1[s1,]
      #   
      #   hat_omega0=solve(cov(MyDataSubIter$Data0))
      #   hat_omega1=solve(cov(MyDataSubIter$Data1))
      #   SubDelta=abind(SubDelta, ComputeDelta(omegaX0=hat_omega0, omegaX1=hat_omega1), along=3)
      # }
      # Delta=apply(SubDelta, c(1,2), mean)
      Cum_adj1=Cum_adj1+SelectionMetric(metric=Delta[1:pk[1], 1:pk[1]], q=(Q[[1]])[l], symmetric=TRUE)
      Cum_adj2=Cum_adj2+SelectionMetric(metric=Delta[(pk[1]+1):(pk[1]+pk[2]), (pk[1]+1):(pk[1]+pk[2])], q=(Q[[2]])[l], symmetric=TRUE)
      Cum_adjb=Cum_adjb+SelectionMetric(metric=Delta[1:pk[1], (pk[1]+1):(pk[1]+pk[2])], q=(Q[[3]])[l], symmetric=FALSE)
    }
    Cum_adj=rbind(cbind(Cum_adj1,Cum_adjb), cbind(t(Cum_adjb), Cum_adj2))/K
    return(Cum_adj)
  })
  
  stopCluster(clust)
  
  Output=array(Output, dim = c(p,p,Q_card))
  return(list(Array=Output, Q=Q, K=K))
}


GetStableEdges=function(stability, pi=0.8, q, symmetric=TRUE, blockrows=NULL, blockcols=NULL){
  if (is.null(blockrows)){
    EdgeStab=stability$Array[,,stability$Q==q]
  } else {
    EdgeStab=stability$Array[blockrows,blockcols,stability$Q==q]
  }
  Edgelist=data.frame(node1=row(EdgeStab)[EdgeStab>=pi], 
                      node2=col(EdgeStab)[EdgeStab>=pi])
  if(symmetric){
    Edgelist=Edgelist[Edgelist$node1<Edgelist$node2,]
  }
  return(Edgelist)
}


GetStableMOEdges=function(stability, pi=0.8, q, l, symmetric=TRUE, blockrows=NULL, blockcols=NULL){
  if (is.null(blockrows)){
    EdgeStab=stability$Array[,,stability$Q[[l]]==q]
  } else {
    EdgeStab=stability$Array[blockrows,blockcols,stability$Q[[l]]==q]
  }
  Edgelist=data.frame(node1=row(EdgeStab)[EdgeStab>=pi], 
                      node2=col(EdgeStab)[EdgeStab>=pi])
  if(symmetric){
    Edgelist=Edgelist[Edgelist$node1<Edgelist$node2,]
  }
  return(Edgelist)
}


HypergeomCalib=function(metric, stability, pi=0.95){
  p=ncol(metric)
  N=p*(p-1)/2
  
  pvalues_y=NULL
  for (q in stability$Q){
    hat_theta=SelectionMetric(metric=metric, q=q, symmetric=TRUE)
    
    likely_edges=GetStableEdges(stability = stability, pi = pi, q = q)
    K_1=nrow(likely_edges)
    likely_theta=matrix(0,p,p)
    likely_theta=SetToCoordinate(matrix = likely_theta, rown = likely_edges[,1], coln = likely_edges[,2], value=1)
    likely_theta=likely_theta+t(likely_theta)
    
    tmp=hat_theta+likely_theta
    y=sum(tmp[upper.tri(tmp)]==2)
    pvalue_y=phyper(q=y-1, m=K_1, n=N-K_1, k=q, lower.tail = FALSE, log.p = TRUE)
    pvalues_y=c(pvalues_y, pvalue_y)
  }
  
  return(list(hyper_pvalue=pvalues_y, hat_q=stability$Q[which.min(pvalues_y)], Q=stability$Q))
}


MOHypergeomCalib=function(Data, hat_phi, stability, pi){
  pk=unname(table(Data$Platform))
  L=length(unique(Data$Platform))
  l=matrix(c(1,4,3,2),ncol=2) # to generalise
  # l=matrix(c(1,3,3,2),ncol=2)
  l_rows_lower=t(matrix(c(1,1,pk[1]+1, pk[1]+1),ncol=2))
  l_rows_upper=t(matrix(c(pk[1],pk[1],pk[1]+pk[2], pk[1]+pk[2]),ncol=2))
  symm=matrix(c(TRUE, FALSE, FALSE, TRUE), ncol=2) # to generalise
  
  p1=matrix(rep(pk, 2), ncol=2) # to generalise
  p2=t(p1) # to generalise
  Nblock=sapply(1:L, FUN=function(i){sapply(1:L, FUN=function(j){ComputeN(p1[i,j], p2[i,j], symm[i,j])})})
  
  hyper_pvalue=NULL
  for (i in 1:L){
    for (j in i:L){ # block-matrix coordinates 
      pvalues_y=NULL
      for (q in stability$Q[[l[i,j]]]){
        hat_theta=SelectionMetric(metric=hat_phi[l_rows_lower[i,j]:l_rows_upper[i,j], l_rows_lower[j,i]:l_rows_upper[j,i]],
                                  q=q, symmetric=symm[i,j])
        # print(pi)
        likely_edges=GetStableMOEdges(stability = stability, pi = pi, q = q, l=l[i,j],
                                      blockrows = l_rows_lower[i,j]:l_rows_upper[i,j], 
                                      blockcols = l_rows_lower[j,i]:l_rows_upper[j,i], symmetric = symm[i,j])
        d_q=nrow(likely_edges)
        likely_theta=matrix(0,p1[i,j],p2[i,j])
        likely_theta=SetToCoordinate(matrix = likely_theta, rown = likely_edges[,1], coln = likely_edges[,2], value=1)
        # print(sum(likely_theta)==nrow(likely_edges))
        tmp=hat_theta+likely_theta
        y=sum(tmp==2)
        pvalue_y=phyper(q=y-1, m=d_q, n=Nblock[i,j]-d_q, k=q, lower.tail = FALSE, log.p = TRUE)
        pvalues_y=c(pvalues_y, pvalue_y)
      }
      hyper_pvalue=rbind(hyper_pvalue, pvalues_y)
    }
  }
  rownames(hyper_pvalue)=l[upper.tri(l,diag=T)]
  # do.call(rbind, stability$Q)
  hat_q=unlist(sapply(1:L, FUN=function(i){sapply(i:L, FUN=function(j){stability$Q[[l[i,j]]][which.min(hyper_pvalue[as.character(l[i,j]),])]})}))
  names(hat_q)=l[upper.tri(l,diag=T)]
  return(list(hyper_pvalue=hyper_pvalue, 
              hat_q=hat_q, 
              Q=stability$Q))
}


MODNHypergeomCalib=function(Data, metric, stability, pi){ # same as MOHypergeom but with metric instead of hat_phi
  pk=unname(table(Data$Platform))
  L=length(unique(Data$Platform))
  l=matrix(c(1,4,3,2),ncol=2) # to generalise
  # l=matrix(c(1,3,3,2),ncol=2)
  l_rows_lower=t(matrix(c(1,1,pk[1]+1, pk[1]+1),ncol=2))
  l_rows_upper=t(matrix(c(pk[1],pk[1],pk[1]+pk[2], pk[1]+pk[2]),ncol=2))
  symm=matrix(c(TRUE, FALSE, FALSE, TRUE), ncol=2) # to generalise
  
  p1=matrix(rep(pk, 2), ncol=2) # to generalise
  p2=t(p1) # to generalise
  Nblock=sapply(1:L, FUN=function(i){sapply(1:L, FUN=function(j){ComputeN(p1[i,j], p2[i,j], symm[i,j])})})
  
  hyper_pvalue=NULL
  for (i in 1:L){
    for (j in i:L){ # block-matrix coordinates 
      pvalues_y=NULL
      for (q in stability$Q[[l[i,j]]]){
        hat_theta=SelectionMetric(metric=metric[l_rows_lower[i,j]:l_rows_upper[i,j], l_rows_lower[j,i]:l_rows_upper[j,i]],
                                  q=q, symmetric=symm[i,j])
        # print(pi)
        likely_edges=GetStableMOEdges(stability = stability, pi = pi, q = q, l=l[i,j],
                                      blockrows = l_rows_lower[i,j]:l_rows_upper[i,j], 
                                      blockcols = l_rows_lower[j,i]:l_rows_upper[j,i], symmetric = symm[i,j])
        d_q=nrow(likely_edges)
        likely_theta=matrix(0,p1[i,j],p2[i,j])
        likely_theta=SetToCoordinate(matrix = likely_theta, rown = likely_edges[,1], coln = likely_edges[,2], value=1)
        # print(sum(likely_theta)==nrow(likely_edges))
        tmp=hat_theta+likely_theta
        y=sum(tmp==2)
        pvalue_y=phyper(q=y-1, m=d_q, n=Nblock[i,j]-d_q, k=q, lower.tail = FALSE, log.p = TRUE)
        pvalues_y=c(pvalues_y, pvalue_y)
      }
      hyper_pvalue=rbind(hyper_pvalue, pvalues_y)
    }
  }
  rownames(hyper_pvalue)=l[upper.tri(l,diag=T)]
  # do.call(rbind, stability$Q)
  hat_q=unlist(sapply(1:L, FUN=function(i){sapply(i:L, FUN=function(j){stability$Q[[l[i,j]]][which.min(hyper_pvalue[as.character(l[i,j]),])]})}))
  names(hat_q)=l[upper.tri(l,diag=T)]
  return(list(hyper_pvalue=hyper_pvalue, 
              hat_q=hat_q, 
              Q=stability$Q))
}


Heatmap=function(matrix, bound, filename=NA, numbers=FALSE, width=7, height=7, 
                 show_colnames=FALSE, show_rownames=FALSE){
  Tmp=matrix
  Tmp[lower.tri(Tmp, diag = T)]=0
  if (numbers){
    TmpNumbers=round(matrix, digits=2)
    TmpNumbers[lower.tri(TmpNumbers, diag = T)]=""
    pheatmap(Tmp, cluster_rows = FALSE, cluster_cols = FALSE, width = width, height = height,
             breaks=seq(-bound,bound,length.out=100), border=NA, 
             show_rownames = show_rownames, show_colnames = show_colnames, filename = filename, display_numbers = TmpNumbers,
             color=colorRampPalette(c('navy', 'blue', 'skyblue', 'white', 'gold', 'red', 'darkred'))(100))
  } else {
    pheatmap(Tmp, cluster_rows = FALSE, cluster_cols = FALSE, width = width, height = height,
             breaks=seq(-bound,bound,length.out=100), border=NA, 
             show_rownames = show_rownames, show_colnames = show_colnames, filename = filename, 
             color=colorRampPalette(c('navy', 'blue', 'skyblue', 'white', 'gold', 'red', 'darkred'))(100))
  }
}


SetToCoordinate=function(matrix, rown, coln, value){
  matrix[rown+nrow(matrix)*(coln-1)]=value
  return(matrix)
}


CalibrationPlot=function(hyperpvalues, theta_star=NULL){
  par(mar=c(4.2,4.2,3,1))
  colors=rep('grey30', length(hyperpvalues$Q))
  colors[hyperpvalues$Q==hyperpvalues$hat_q]='blue'
  if (is.null(theta_star)){
    plot(-hyperpvalues$hyper_pvalue, pch=19, xaxt='n', las=1, ylab=expression(italic(p[value])), xlab='',
         main=substitute(hat(q)~'='~q_hat1, list(q_hat1=hyperpvalues$hat_q)), col=colors)
  } else {
    tmp=theta_star
    q_star=sum(tmp[upper.tri(tmp)])
    colors[which.min(abs(hyperpvalues$Q-q_star))]='red'
    mypch=rep(19,length(hyperpvalues$hyper_pvalue))
    mypch[colors!='grey30']=8
    mycex=rep(0.7,length(hyperpvalues$hyper_pvalue))
    mycex[colors!='grey30']=1.2
    mylwd=rep(1,length(hyperpvalues$hyper_pvalue))
    mylwd[colors!='grey30']=2
    plot(-hyperpvalues$hyper_pvalue, xaxt='n', las=1, ylab=expression(italic(p[value])), xlab=expression(italic(q)), 
         cex=mycex, pch=mypch, lwd=mylwd,
         main=substitute('q*='~q_star1*';'~hat(q)~'='~q_hat1, 
                         list(q_star1=q_star, q_hat1=hyperpvalues$hat_q)), col=colors, cex.lab=1.5)
  }
  # abline(v=seq(1,length(hyperpvalues$Q),by=10),lty=3)
  # axis(side=1, at=seq(1,length(hyperpvalues$Q),by=10), 
  #      labels=(hyperpvalues$Q)[seq(1,length(hyperpvalues$Q),by=10)])
}


MOCalibrationPlot=function(Data, hyperpvalues, theta_star=NULL){
  L=length(unique(Data$Platform))
  par(mfrow=c(1,L+L*(L-1)/2))
  
  l=matrix(c(1,4,3,2),ncol=2)
  l_rows_lower=t(matrix(c(1,1,pk[1]+1, pk[1]+1),ncol=2))
  l_rows_upper=t(matrix(c(pk[1],pk[1],pk[1]+pk[2], pk[1]+pk[2]),ncol=2))
  symm=matrix(c(TRUE, FALSE, FALSE, TRUE), ncol=2)
  
  for (i in 1:L){
    for (j in i:L){
      colors=rep('grey30', length(hyperpvalues$Q[[l[i,j]]]))
      colors[hyperpvalues$Q[[l[i,j]]]==hyperpvalues$hat_q[as.character(l[i,j])]]='blue'
      if (is.null(theta_star)){
        plot(-hyperpvalues$hyper_pvalue[as.character(l[i,j]),], pch=19, xaxt='n', las=1, ylab=expression(italic(p[value])), xlab='',
             main=substitute('Block'~block*':'~hat(q)~'='~q_hat1, 
                             list(block=l[i,j], q_hat1=hyperpvalues$hat_q[as.character(l[i,j])])), col=colors, cex.axis=1.5)
      } else {
        if (symm[i,j]){
          tmp=theta_star[l_rows_lower[i,j]:l_rows_upper[i,j], l_rows_lower[j,i]:l_rows_upper[j,i]]
          q_star=sum(tmp[upper.tri(tmp)])
        } else {
          tmp=theta_star[l_rows_lower[i,j]:l_rows_upper[i,j], l_rows_lower[j,i]:l_rows_upper[j,i]]
          q_star=sum(tmp)
        }
        colors[which.min(abs(hyperpvalues$Q[[l[i,j]]]-q_star))]='red'
        plot(-hyperpvalues$hyper_pvalue[as.character(l[i,j]),], pch=19, xaxt='n', las=1, ylab=expression(italic(p[value])), xlab='',
             main=substitute('Block'~block*':'~'q*='~q_star1*';'~hat(q)~'='~q_hat1, 
                             list(block=l[i,j], q_star1=q_star, q_hat1=hyperpvalues$hat_q[as.character(l[i,j])])), col=colors)
      }
      abline(v=seq(1,length(hyperpvalues$Q[[l[i,j]]]),by=5),lty=3)
      axis(side=1, at=seq(1,length(hyperpvalues$Q[[l[i,j]]]),by=5), 
           labels=(hyperpvalues$Q[[l[i,j]]])[seq(1,length(hyperpvalues$Q[[l[i,j]]]),by=5)])
    }
  }
}


ConditionalCovariance=function(omega){
  p=ncol(omega)
  kappa=sapply(1:p, FUN=function(j){sapply(1:p, FUN=function(i){-omega[i,j]/(omega[i,i]*omega[j,j]-omega[i,j]^2)})})
  diag(kappa)=0
  return(kappa)
}


KappaTrace=function(kappa, tau){
  d=det(kappa)
  return(2/d*(kappa[1,1]*kappa[2,2]-tau*kappa[1,2]))
}


KappaLogDet=function(kappa, tau){
  d=det(kappa)
  return(log(kappa[1,1]*kappa[2,2]-tau^2)-log(d))
}


ComputeDelta=function(omegaX0, omegaX1){
  p=ncol(omegaX0)
  Delta=matrix(0,p,p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      omega0=rbind(c(omegaX0[i,i], omegaX0[i,j]), c(omegaX0[i,j], omegaX0[j,j]))
      omega1=rbind(c(omegaX1[i,i], omegaX1[i,j]), c(omegaX1[i,j], omegaX1[j,j]))
      kappa0=solve(omega0)
      kappa1=solve(omega1)
      
      tau=c(kappa0[1,2], kappa1[1,2])[which.min(abs(c(kappa0[1,2], kappa1[1,2])))]
      # tau=c(kappa0[1,2], kappa1[1,2])[which.max(abs(c(kappa0[1,2], kappa1[1,2])))]
      # tau=mean(c(kappa0[1,2], kappa1[1,2]))
      Delta[i,j]=KappaTrace(kappa=kappa0, tau=tau)+KappaTrace(kappa=kappa1, tau=tau)-
        KappaLogDet(kappa=kappa0, tau=tau)-KappaLogDet(kappa=kappa1, tau=tau)-4
    }
  }
  Delta=Delta+t(Delta)
  return(Delta)
}


GetPerformance=function(theta){
  p=ncol(theta)
  N=p*(p-1)/2
  
  TP=sum(theta[upper.tri(theta)]==3)
  FN=sum(theta[upper.tri(theta)]==2)
  FP=sum(theta[upper.tri(theta)]==1)
  TN=sum(theta[upper.tri(theta)]==0)
  
  sensitivity=TP/(TP+FN)
  specificity=TN/(TN+FP)
  precision=TP/(TP+FP)
  accuracy=(TP+TN)/N
  return(list(TP=TP, FN=FN, FP=FP, TN=TN, 
              sensitivity=sensitivity, specificity=specificity, 
              precision=precision, accuracy=accuracy))
}


CalibrationGrid=function(object, pi_list){
  hyperpvalues=NULL
  for (pi in pi_list){
    print(pi)
    hyperpvalues=cbind(hyperpvalues, -HypergeomCalib(metric=object$estimate$metric, stability = object$stability, pi = pi)$hyper_pvalue)
  }
  colnames(hyperpvalues)=pi_list
  rownames(hyperpvalues)=object$parameters$Q
  return(hyperpvalues)
}


MOIndivGraph=function(Data, metric, q, status, graph=TRUE, test=FALSE, pal=NULL){
  if (is.null(pal)){
    pal=c(brewer.pal(12,name="Set3"), brewer.pal(n=8, name="Set2"), brewer.pal(n=4, name="Set1"))
  }
  
  n=sum(Data$Condition==status)
  
  pk=unname(table(Data$Platform))
  p1=matrix(rep(pk, 2), ncol=2) # to generalise
  p2=t(p1) # to generalise
  L=length(unique(Data$Platform))
  l=matrix(c(1,4,3,2),ncol=2) # to generalise
  l_rows_lower=t(matrix(c(1,1,pk[1]+1, pk[1]+1),ncol=2))
  l_rows_upper=t(matrix(c(pk[1],pk[1],pk[1]+pk[2], pk[1]+pk[2]),ncol=2))
  symm=matrix(c(TRUE, FALSE, FALSE, TRUE), ncol=2) # to generalise
  
  adjacency_list=list()
  edgelist_list=list()
  for (i in 1:L){
    for (j in i:L){
      tmpmetric=metric[l_rows_lower[i,j]:l_rows_upper[i,j], l_rows_lower[j,i]:l_rows_upper[j,i]]
      if (symm[i,j]){
        MyEdgelist_raw=data.frame(node1=row(tmpmetric)[upper.tri(tmpmetric)], node2=col(tmpmetric)[upper.tri(tmpmetric)],
                                  cor=tmpmetric[upper.tri(tmpmetric)], stringsAsFactors = FALSE)
      } else {
        MyEdgelist_raw=data.frame(node1=as.vector(row(tmpmetric)), node2=as.vector(col(tmpmetric)), 
                                  cor=as.vector(tmpmetric), stringsAsFactors = FALSE)
      }
      MyEdgelist=MyEdgelist_raw[sort.list(abs(MyEdgelist_raw$cor), decreasing = TRUE),]
      
      cor_q=abs(MyEdgelist$cor[q[as.character(l[i,j])]])
      
      adjacency=matrix(0,p1[i,j],p2[i,j])
      if (symm[i,j]){
        adjacency[upper.tri(adjacency)]=ifelse(abs(tmpmetric[upper.tri(tmpmetric)])>=cor_q, yes=1, no=0)
        adjacency=adjacency+t(adjacency)
      } else {
        adjacency=ifelse(abs(tmpmetric)>=cor_q, yes=1, no=0)
        diag(adjacency)=0
      }
      
      if (test){
        z_list=zTransform(MyEdgelist$cor)
        z_pval=sapply(1:length(z_list), FUN=function(i){Zpval(z_list[i], n)})
        MyEdgelist=cbind(MyEdgelist, z_list, z_pval)
        colnames(MyEdgelist)[4:5]=c('z_cor', 'z_pvalue')
      }
      adjacency_list=c(adjacency_list, list(adjacency))
      edgelist_list=c(edgelist_list, list(MyEdgelist))
    }
  }
  
  adjacency=rbind(cbind(adjacency_list[[1]], adjacency_list[[2]]), cbind(t(adjacency_list[[2]]), adjacency_list[[3]]))
  colnames(adjacency)=rownames(adjacency)=colnames(Data$Data)
  
  if (graph){
    gg=GetGraph(Data=Data, adjacency=adjacency, pal=pal)
    res=list(Graph=gg, adjacency=adjacency, Edgelist=edgelist_list)
  } else {
    res=list(adjacency=adjacency, Edgelist=edgelist_list)
  }
  
  return(res)
}


MODiffGraph=function(Data, metric, q, graph=TRUE, test=FALSE, pal=NULL){
  if (is.null(pal)){
    pal=c(brewer.pal(12,name="Set3"), brewer.pal(n=8, name="Set2"), brewer.pal(n=4, name="Set1"))
  }
  
  n=nrow(Data$Data0)
  
  pk=unname(table(Data$Platform))
  p1=matrix(rep(pk, 2), ncol=2) # to generalise
  p2=t(p1) # to generalise
  L=length(unique(Data$Platform))
  l=matrix(c(1,4,3,2),ncol=2) # to generalise
  l_rows_lower=t(matrix(c(1,1,pk[1]+1, pk[1]+1),ncol=2))
  l_rows_upper=t(matrix(c(pk[1],pk[1],pk[1]+pk[2], pk[1]+pk[2]),ncol=2))
  symm=matrix(c(TRUE, FALSE, FALSE, TRUE), ncol=2) # to generalise
  
  adjacency_list=list()
  edgelist_list=list()
  for (i in 1:L){
    for (j in i:L){
      tmpmetric=metric[l_rows_lower[i,j]:l_rows_upper[i,j], l_rows_lower[j,i]:l_rows_upper[j,i]]
      if (symm[i,j]){
        MyEdgelist_raw=data.frame(node1=row(tmpmetric)[upper.tri(tmpmetric)], node2=col(tmpmetric)[upper.tri(tmpmetric)],
                                  cor=tmpmetric[upper.tri(tmpmetric)], stringsAsFactors = FALSE)
      } else {
        MyEdgelist_raw=data.frame(node1=as.vector(row(tmpmetric)), node2=as.vector(col(tmpmetric)), 
                                  cor=as.vector(tmpmetric), stringsAsFactors = FALSE)
      }
      MyEdgelist=MyEdgelist_raw[sort.list(abs(MyEdgelist_raw$cor), decreasing = TRUE),]
      
      cor_q=abs(MyEdgelist$cor[q[as.character(l[i,j])]])
      
      adjacency=matrix(0,p1[i,j],p2[i,j])
      if (symm[i,j]){
        adjacency[upper.tri(adjacency)]=ifelse(abs(tmpmetric[upper.tri(tmpmetric)])>=cor_q, yes=1, no=0)
        adjacency=adjacency+t(adjacency)
      } else {
        adjacency=ifelse(abs(tmpmetric)>=cor_q, yes=1, no=0)
        diag(adjacency)=0
      }
      
      if (test){
        chi_pval=sapply(MyEdgelist[,3], FUN=function(x){
          pval=1-pchisq(n*x, df=1)
          pval=ifelse(pval==0,
                      yes=-pchisq(n*x, df=1, log.p = TRUE), no=pval)
          return(pval)})
        MyEdgelist=cbind(MyEdgelist, chi_pval)
        colnames(MyEdgelist)[4]='chi_pval'
      }
      adjacency_list=c(adjacency_list, list(adjacency))
      edgelist_list=c(edgelist_list, list(MyEdgelist))
    }
  }
  
  adjacency=rbind(cbind(adjacency_list[[1]], adjacency_list[[2]]), cbind(t(adjacency_list[[2]]), adjacency_list[[3]]))
  colnames(adjacency)=rownames(adjacency)=colnames(Data$Data0)
  
  if (graph){
    gg=GetGraph(Data=Data, adjacency=adjacency, pal=pal)
    res=list(Graph=gg, adjacency=adjacency, Edgelist=edgelist_list)
  } else {
    res=list(adjacency=adjacency, Edgelist=edgelist_list)
  }
  
  return(res)
}

