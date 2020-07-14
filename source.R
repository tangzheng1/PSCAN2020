
ACAT<-function(Pvals,Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals<0)+sum(Pvals>1))>0){
    stop("P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(Pvals),length(Pvals))
  }else if (length(Weights)!=length(Pvals)){
    stop("The length of weights should be the same as that of the p-values")
  }else if (sum(Weights<0)>0){
    stop("All the weights must be positive!")
  }else{
    Weights<-Weights/sum(Weights)
  }
  
  
  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}



library(MASS)
library("SKAT")
library("igraph")
library(seqMeta)  

######## scan Test 3D
## test for mean
scanTest2 <- function(U, V, MAC, weight=NULL, cords.ls, nperm=NULL, signal.alpha=NULL){
  
  MAC.LB = 10
  if(sum(MAC) < MAC.LB){  
    
    Q.max = NA; pval=NA; signal.region =NA; C = NA; Q.scan2.gz = NA
    
  }else{
    
    m = length(U)
    n.dom = length(cords.ls)
    m.ls = sapply(cords.ls, nrow)
    m.mapro = sum(m.ls)
    m.imapro = m - m.mapro
    
    flag.search = as.numeric(m.ls>1)
    
    
    cluster.list = NULL
    
    m.idx = 0
    for(d in 1:n.dom){
      
      if(flag.search[d]){
        
        cords.dist = as.matrix(dist(cords.ls[[d]]))
        dist.uniq = sort(unique(  as.numeric(cords.dist[lower.tri(cords.dist)])  ))
        n.poss = length(dist.uniq)
        cluster.list.d = as.list(c(1:m.ls[d])) 
        cluster.pre = cluster.list.d
        cluster.num.pre = m.ls[d]
        
        for(j in 1:n.poss){
          tmp = clusters(graph_from_adjacency_matrix(0+(cords.dist<=dist.uniq[j]), mode="undirected"))
          cluster.num.cur = tmp$no
          
          if(cluster.num.cur == cluster.num.pre){
            next
          }else{
            cluster.cur = split(1:m.ls[d], tmp$membership)
            cluster.add = cluster.cur[!(cluster.cur %in% cluster.pre)]
            cluster.list.d = c(cluster.list.d, cluster.add)
            
            cluster.pre = cluster.cur
            cluster.num.pre = cluster.num.cur
          }
          if(cluster.num.cur ==1){
            break
          }
        }
        
        for(j in 1:length(cluster.list.d)){
          cluster.list.d[[j]] = m.idx+cluster.list.d[[j]]
        }
        
        cluster.list = c(cluster.list, cluster.list.d)
        
      }else{
        
        cluster.list = c(cluster.list, list(m.idx + 1))
      }
      
      m.idx = m.idx + m.ls[d]
      
    }
    
    
    if( m.imapro ){
      
      cluster.list = c(cluster.list, list( (1:m)[-c(1:m.mapro)] ) )
      
    }
    
    
    n.scan = length(cluster.list)
    C = matrix(0, nrow=n.scan, ncol=m)
    for(j in 1:n.scan){
      C[j, cluster.list[[j]]] = 1
    }
    
    if( length( which(rowSums(C)==m) )==0 ){
      C = rbind(C, rep(1, m))
      n.scan = n.scan + 1
      cluster.list = c(cluster.list, list(rep(1:m)) )
    }
    
    
    C.MAC = C %*% MAC
    rm.idx = which(C.MAC < MAC.LB)
    if(length(rm.idx)>0){
      C = C[-rm.idx,,drop=FALSE]
      n.scan = n.scan - length(rm.idx)
      cluster.list = cluster.list[-rm.idx]
    }
    
    
    if(!is.null(weight)){
      
      C = t(t(C) * weight)
      
    }
    
    UC = C %*% U
    VC = C %*% V %*% t(C)
    Q = UC*UC / diag(VC)
    Q.max =  max(Q)
    
    Qp = 1-pchisq(Q, 1)
    pval = ACAT( Qp )
    
    
   
    signal.region = NULL
    
    if(!is.null(signal.alpha)){
      
      set.seed(11)
      if(is.null(nperm)){
        nperm = 50/signal.alpha
      }
      Q.max.perm = NULL
      for(l in 1:nperm){
        U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        tmp = C %*% U.perm
        Q.max.perm = c(Q.max.perm, max( tmp*tmp/diag(VC) ) )
      }

  
      signal.thresh = quantile(Q.max.perm, 1-signal.alpha)
      cand.idx = which(Q >= signal.thresh)
      
      if( length(cand.idx)>0 ){
        # Q.max.idx is the idx for the cluster that gives maximum test stat. If there is tie, the cluster with smallest region will be selected
        Q.max.idx = which(Q==Q.max)
        if( length(Q.max.idx)>1 ){
          Q.max.idx = Q.max.idx[which.min(sapply(cluster.list[Q.max.idx], length))]
        }
        signal.region[[1]] = cluster.list[[Q.max.idx]]
        signal.cand = cluster.list[cand.idx] 
        Q.cand = Q[cand.idx]
        signal.cand.max.region =  cluster.list[[Q.max.idx]]
        nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) ) == 0) )
        while( length(nonoverlap.idx)>0  ){
          signal.cand = signal.cand[nonoverlap.idx] 
          Q.cand = Q.cand[nonoverlap.idx]
          Q.max.idx = which(Q.cand==max(Q.cand))
          if( length(Q.max.idx)>1 ){
            Q.max.idx = Q.max.idx[which.min(sapply(signal.cand[Q.max.idx], length))]
          }
          signal.region = c(signal.region , signal.cand[Q.max.idx])
          signal.cand.max.region = signal.cand[[Q.max.idx]]
          nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) ) == 0) )
        }
      }
    }
    
    Q.scan2.gz =  UC / sqrt(diag(VC)) 
    
  }

  
  return(list(Q.scan2.stat = Q.max, Q.scan2.pval=pval, signal.region = signal.region, C = C, Q.scan2.gz = Q.scan2.gz ))
  
}

## test for variance
scanTest3 <- function(U, V, MAC, weight=NULL, cords.ls, nperm=NULL, signal.alpha = NULL){
  
  
  MAC.LB = 10
  if(sum(MAC) < MAC.LB){  
  
    Qp.min = NA; pval= NA; signal.region=NA; C=NA; Qp=NA;
   
  }else{
    m = length(U)
    n.dom = length(cords.ls)
    m.ls = sapply(cords.ls, nrow)
    m.mapro = sum(m.ls)
    m.imapro = m - m.mapro
    
    flag.search = as.numeric(m.ls>1)
    
    
    cluster.list = NULL
    
    m.idx = 0
    for(d in 1:n.dom){
      
      if(flag.search[d]){
        
        cords.dist = as.matrix(dist(cords.ls[[d]]))
        dist.uniq = sort(unique(  as.numeric(cords.dist[lower.tri(cords.dist)])  ))
        n.poss = length(dist.uniq)
        cluster.list.d = as.list(c(1:m.ls[d])) 
        cluster.pre = cluster.list.d
        cluster.num.pre = m.ls[d]
        
        for(j in 1:n.poss){
          tmp = clusters(graph_from_adjacency_matrix(0+(cords.dist<=dist.uniq[j]), mode="undirected"))
          cluster.num.cur = tmp$no
          
          if(cluster.num.cur == cluster.num.pre){
            next
          }else{
            cluster.cur = split(1:m.ls[d], tmp$membership)
            cluster.add = cluster.cur[!(cluster.cur %in% cluster.pre)]
            cluster.list.d = c(cluster.list.d, cluster.add)
            
            cluster.pre = cluster.cur
            cluster.num.pre = cluster.num.cur
          }
          if(cluster.num.cur ==1){
            break
          }
        }
        
        for(j in 1:length(cluster.list.d)){
          cluster.list.d[[j]] = m.idx+cluster.list.d[[j]]
        }
        
        cluster.list = c(cluster.list, cluster.list.d)
        
      }else{
        
        cluster.list = c(cluster.list, list(m.idx + 1))
      }
      
      m.idx = m.idx + m.ls[d]
      
    }
    
    
    if( m.imapro ){
      
      cluster.list = c(cluster.list, list( (1:m)[-c(1:m.mapro)] ) )
      
    }
    
    n.scan = length(cluster.list)
    C = matrix(0, nrow=n.scan, ncol=m)
    for(j in 1:n.scan){
      C[j, cluster.list[[j]]] = 1
    }
    
    if( length( which(rowSums(C)==m) )==0 ){
      C = rbind(C, rep(1, m))
      n.scan = n.scan + 1
      cluster.list = c(cluster.list, list(rep(1:m)) )
    }
    
    
    C.MAC = C %*% MAC
    rm.idx = which(C.MAC < MAC.LB)
    if(length(rm.idx)>0){
      C = C[-rm.idx,,drop=FALSE]
      n.scan = n.scan - length(rm.idx)
      cluster.list = cluster.list[-rm.idx]
    }
    
    if(!is.null(weight)){
      
      C = t(t(C) * weight)
      
    }
    
    
    Qp = NULL
    for(j in 1:n.scan){
      idx = which(C[j,]!=0)
      Qp = c(Qp, skatTest(U[idx], V[idx,idx], weight=C[j,idx], r = 0)$Q.ker.optim.pval )
    }
    Qp[Qp<0] = 0 
    Qp.min = min(Qp)
    
    pval = ACAT( Qp )
    
   
    
    
    signal.region = NULL
    
    if(!is.null(signal.alpha)){
      
      set.seed(11)
      if(is.null(nperm)){
        nperm = 50/signal.alpha
      }
      Qp.min.perm = NULL
      for(l in 1:nperm){
        U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        
        Qp.perm = NULL
        for(j in 1:n.scan){
          idx = which(C[j,]!=0)
          Qp.perm = c(Qp.perm, skatTest(U.perm[idx], V[idx,idx], weight=C[j,idx], r = 0)$Q.ker.optim.pval )
        }
        Qp.perm[Qp.perm<0] = 0  
        Qp.min.perm = c(Qp.min.perm, min(Qp.perm) )
      }
      
      signal.thresh = quantile(Qp.min.perm, signal.alpha)
      cand.idx = which(Qp <= signal.thresh)
      
      if( length(cand.idx)>0 ){
        # Qp.min.idx is the idx for the cluster that gives minimum test pval. If there is tie, the cluster with smallest region will be selected
        Qp.min.idx = which(Qp==Qp.min)
        if( length(Qp.min.idx)>1 ){
          Qp.min.idx = Qp.min.idx[which.min(sapply(cluster.list[Qp.min.idx], length))]
        }
        signal.region[[1]] = cluster.list[[Qp.min.idx]]
        signal.cand = cluster.list[cand.idx] 
        Qp.cand = Qp[cand.idx]
        signal.cand.min.region =  cluster.list[[Qp.min.idx]]
        nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) ) == 0) )
        while( length(nonoverlap.idx)>0  ){
          signal.cand = signal.cand[nonoverlap.idx] 
          Qp.cand = Qp.cand[nonoverlap.idx]
          Qp.min.idx = which(Qp.cand==min(Qp.cand))        
          if( length(Qp.min.idx)>1 ){
            Qp.min.idx = Qp.min.idx[which.min(sapply(signal.cand[Qp.min.idx], length))]
          }
          signal.region = c(signal.region , signal.cand[Qp.min.idx])
          signal.cand.min.region = signal.cand[[Qp.min.idx]]
          nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) ) == 0) )
        }
      }
      
    }
    
  }

  
  return(list(Q.scan3.stat = Qp.min, Q.scan3.pval=pval, signal.region = signal.region, C=C, Q.scan3.gzpval = Qp))
  
}


######## Single-variant
scanTest.sv <- function(U, V, weight=NULL, cords.ls, nperm=1000, signal.alpha=0.05){
  
  m = length(U)

  cluster.list = as.list(1:m)
  
  n.scan = length(cluster.list)
  C = matrix(0, nrow=n.scan, ncol=m)
  for(j in 1:n.scan){
    C[j, cluster.list[[j]]] = 1
  }
  
  if(!is.null(weight)){
    
    C = t(t(C) * weight)
    
  }
  
  UC = C %*% U
  VC = C %*% V %*% t(C)
  Q = UC*UC / diag(VC)
  Q.max =  max(Q)

  
  ## permutation test
  
  Q.max.perm = NULL
  
  set.seed(11)
  for(l in 1:nperm){
    U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    tmp = C %*% U.perm
    Q.max.perm = c(Q.max.perm, max( tmp*tmp/diag(VC) ) )
  }
  
  
  pval =  sum(Q.max.perm >= Q.max)/nperm
  
  signal.region = NULL
  
  if(!is.null(signal.alpha)){
    
    signal.thresh = quantile(Q.max.perm, 1-signal.alpha)
    cand.idx = which(Q >= signal.thresh)
    
    if( length(cand.idx)>0 ){
      signal.region = which(Q >= signal.thresh)
    }
    
  }
  
  
  
  return(list(Q.scan2.stat = Q.max, Q.scan2.pval=pval, signal.region = signal.region, C = C, Q.scan2.gz = UC / sqrt(diag(VC)) ))
  
}



######## scan Test 1D
## test for mean
scanTest2.1D <- function(U, V, MAC, weight=NULL, POS, nperm=NULL, signal.alpha=NULL){
  
  MAC.LB = 10
  if(sum(MAC) < MAC.LB){  
    
    Q.max = NA; pval=NA; signal.region =NA; C = NA; Q.scan2.gz = NA
    
  }else{
    
    m = length(U)
    pos = POS - POS[1]
    pos.max = max(pos)
    W.size = round(pos.max/c(20, 16, 12, 8, 4, 2))
    skip.length = round(min(W.size)/2)
    
    
    C = NULL
    cluster.list = list()
    cluster.mac = NULL
    
    n.size = length(W.size)
    for(j in 1:n.size){
      
      pos.a = 0; pos.b = W.size[j]
      mac.prior = 0
      flag = 1
      
      while(flag>0){
        
        tmp.idx = which(pos>pos.a & pos<=pos.b)
        tmp.C = as.numeric(pos>pos.a & pos<=pos.b)
        tmp.mac = sum(MAC[tmp.idx])
        if( tmp.mac>=MAC.LB & abs(tmp.mac - mac.prior)>0 ){ # remove window with less than 10 MAC, and at least 1 MAC difference between the prior and current window
          C = rbind(C, tmp.C)
          cluster.mac = c(cluster.mac, tmp.mac)
          mac.prior = tmp.mac
          cluster.list = c(cluster.list, list(tmp.idx) )
        }
        
        if(pos.b == pos.max){
          flag = 0
        }
        
        pos.a = pos.a + skip.length
        pos.b = pos.b + skip.length
        
        if(pos.b>pos.max){
          pos.b = pos.max
        }
        
      }
      
    }
    
    
    n.scan = length(cluster.list)
    
    
    if(!is.null(weight)){
      
      C = t(t(C) * weight)
      
    }
    
    UC = C %*% U
    VC = C %*% V %*% t(C)
    Q = UC*UC / diag(VC)
    Q.max =  max(Q)
    
    Qp = 1-pchisq(Q, 1)
    pval = ACAT( Qp )
    
    
    signal.region = NULL
    
    if(!is.null(signal.alpha)){
      
      set.seed(11)
      if(is.null(nperm)){
        nperm = 50/signal.alpha
      }
      Q.max.perm = NULL
      for(l in 1:nperm){
        U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        tmp = C %*% U.perm
        Q.max.perm = c(Q.max.perm, max( tmp*tmp/diag(VC) ) )
      }
     
      signal.thresh = quantile(Q.max.perm, 1-signal.alpha)
      cand.idx = which(Q >= signal.thresh)
      
      if( length(cand.idx)>0 ){
        # Q.max.idx is the idx for the cluster that gives maximum test stat. If there is tie, the cluster with smallest region will be selected
        Q.max.idx = which(Q==Q.max)
        if( length(Q.max.idx)>1 ){
          Q.max.idx = Q.max.idx[which.min(sapply(cluster.list[Q.max.idx], length))]
        }
        signal.region[[1]] = cluster.list[[Q.max.idx]]
        signal.cand = cluster.list[cand.idx] 
        Q.cand = Q[cand.idx]
        signal.cand.max.region =  cluster.list[[Q.max.idx]]
        nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) ) == 0) )
        while( length(nonoverlap.idx)>0  ){
          signal.cand = signal.cand[nonoverlap.idx] 
          Q.cand = Q.cand[nonoverlap.idx]
          Q.max.idx = which(Q.cand==max(Q.cand))
          if( length(Q.max.idx)>1 ){
            Q.max.idx = Q.max.idx[which.min(sapply(signal.cand[Q.max.idx], length))]
          }
          signal.region = c(signal.region , signal.cand[Q.max.idx])
          signal.cand.max.region = signal.cand[[Q.max.idx]]
          nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) ) == 0) )
        }
      }
    }
    
    Q.scan2.gz =  UC / sqrt(diag(VC)) 
    
  }
  
  
  return(list(Q.scan2.stat = Q.max, Q.scan2.pval=pval, signal.region = signal.region, C = C, Q.scan2.gz = Q.scan2.gz ))
  
}

## test for variance
scanTest3.1D <- function(U, V, MAC, weight=NULL, POS, nperm=NULL, signal.alpha = NULL){
  
  
  MAC.LB = 10
  if(sum(MAC) < MAC.LB){  
    
    Qp.min = NA; pval= NA; signal.region=NA; C=NA; Qp=NA;
    
  }else{
    m = length(U)
    
    pos = POS - POS[1]
    pos.max = max(pos)
    W.size = round(pos.max/c(20, 16, 12, 8, 4, 2))
    skip.length = round(min(W.size)/2)
    
    
    C = NULL
    cluster.list = list()
    cluster.mac = NULL
    
    n.size = length(W.size)
    for(j in 1:n.size){
      
      pos.a = 0; pos.b = W.size[j]
      mac.prior = 0
      flag = 1
      
      while(flag>0){
        tmp.idx = which(pos>pos.a & pos<=pos.b)
        tmp.C = as.numeric(pos>pos.a & pos<=pos.b)
        tmp.mac = sum(MAC[tmp.idx])
        if( tmp.mac>=MAC.LB & abs(tmp.mac - mac.prior)>0 ){ # remove window with less than 10 MAC, and at least 1 MAC difference between the prior and current window
          C = rbind(C, tmp.C)
          cluster.mac = c(cluster.mac, tmp.mac)
          mac.prior = tmp.mac
          cluster.list = c(cluster.list, list(tmp.idx) )
        }
        
        if(pos.b == pos.max){
          flag = 0
        }
        
        pos.a = pos.a + skip.length
        pos.b = pos.b + skip.length
        
        if(pos.b>pos.max){
          pos.b = pos.max
        }
        
      }
      
    }
    
    
    n.scan = length(cluster.list)
    
    
    if(!is.null(weight)){
      
      C = t(t(C) * weight)
      
    }
    
    # observed stat
    Qp = NULL
    for(j in 1:n.scan){
      idx = which(C[j,]!=0)
      Qp = c(Qp, skatTest(U[idx], V[idx,idx], weight=C[j,idx], r = 0)$Q.ker.optim.pval )
    }
    Qp[Qp<0] = 0 
    Qp.min = min(Qp)
    
    pval = ACAT( Qp )
    
    
    
    
    signal.region = NULL
    
    if(!is.null(signal.alpha)){
      
      set.seed(11)
      if(is.null(nperm)){
        nperm = 50/signal.alpha
      }
      Qp.min.perm = NULL
      for(l in 1:nperm){
        U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        
        Qp.perm = NULL
        for(j in 1:n.scan){
          idx = which(C[j,]!=0)
          Qp.perm = c(Qp.perm, skatTest(U.perm[idx], V[idx,idx], weight=C[j,idx], r = 0)$Q.ker.optim.pval )
        }
        Qp.perm[Qp.perm<0] = 0  
        Qp.min.perm = c(Qp.min.perm, min(Qp.perm) )
      }
      
      signal.thresh = quantile(Qp.min.perm, signal.alpha)
      cand.idx = which(Qp <= signal.thresh)
      
      if( length(cand.idx)>0 ){
        # Qp.min.idx is the idx for the cluster that gives minimum test pval. If there is tie, the cluster with smallest region will be selected
        Qp.min.idx = which(Qp==Qp.min)
        if( length(Qp.min.idx)>1 ){
          Qp.min.idx = Qp.min.idx[which.min(sapply(cluster.list[Qp.min.idx], length))]
        }
        signal.region[[1]] = cluster.list[[Qp.min.idx]]
        signal.cand = cluster.list[cand.idx] 
        Qp.cand = Qp[cand.idx]
        signal.cand.min.region =  cluster.list[[Qp.min.idx]]
        nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) ) == 0) )
        while( length(nonoverlap.idx)>0  ){
          signal.cand = signal.cand[nonoverlap.idx] 
          Qp.cand = Qp.cand[nonoverlap.idx]
          Qp.min.idx = which(Qp.cand==min(Qp.cand))        
          if( length(Qp.min.idx)>1 ){
            Qp.min.idx = Qp.min.idx[which.min(sapply(signal.cand[Qp.min.idx], length))]
          }
          signal.region = c(signal.region , signal.cand[Qp.min.idx])
          signal.cand.min.region = signal.cand[[Qp.min.idx]]
          nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) ) == 0) )
        }
      }
      
    }
    
  }
  
  
  return(list(Q.scan3.stat = Qp.min, Q.scan3.pval=pval, signal.region = signal.region, C=C, Q.scan3.gzpval = Qp))
  
}

################################# 
Beta.Weights<-function(MAF,weights.beta, Cutoff=1, Is.MAF=TRUE){
  # no change
  n<-length(MAF)
  weights<-rep(0,n)
  Sign<-rep(1,n)
  #print(MAF)
  
  IDX1<-which(MAF > 0.5)
  if(length(IDX1) > 0){
    Sign[IDX1]<--1
    MAF[IDX1]<-1-MAF[IDX1]
  }
  
  IDX_0<-union(which(MAF == 0), which(MAF > Cutoff))
  if(length(IDX_0) == n){
    #stop("No polymorphic SNPs")
    weights<-rep(0,n)
  } else if( length(IDX_0) == 0){
    weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
  } else {
    weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
  }
  
  weights = weights * Sign	
  return(weights)
  
}
SKAT_META_Optimal_Get_Q<-function(Score, r.all){
  # no change
  n.r<-length(r.all)
  Q.r<-rep(0,n.r)
  
  for(i in 1:n.r){
    r.corr<-r.all[i]
    Q.r[i]<-(1-r.corr) * sum(Score^2) + r.corr * sum(Score)^2
  }
  Q.r = Q.r /2
  
  re<-list(Q.r=Q.r)
  return(re)
}
SKAT_META_Optimal_Param<-function(Phi,r.all){
  # no change
  p.m<-dim(Phi)[2]
  r.n<-length(r.all)
  
  # ZMZ
  Z.item1.1<- Phi %*% rep(1,p.m)
  ZZ<-Phi
  ZMZ<- Z.item1.1 %*% t(Z.item1.1) / sum(ZZ)
  
  # W3.2 Term : mixture chisq
  W3.2.t<-ZZ - ZMZ
  lambda<-SKAT:::Get_Lambda(W3.2.t)
  
  # W3.3 Term : variance of remaining ...
  W3.3.item<-sum(ZMZ *(ZZ-ZMZ)) * 4
  
  # tau term 
  z_mean_2<- sum(ZZ)/p.m^2
  tau1<- sum(ZZ %*% ZZ) / p.m^2 / z_mean_2
  
  # Mixture Parameters
  MuQ<-sum(lambda)
  VarQ<-sum(lambda^2) *2 + W3.3.item
  KerQ<-sum(lambda^4)/(sum(lambda^2))^2 * 12
  Df<-12/KerQ
  
  # W3.1 Term : tau1 * chisq_1
  tau<-rep(0,r.n)
  for(i in 1:r.n){
    r.corr<-r.all[i]
    term1<-p.m^2*r.corr * z_mean_2 + tau1 * (1-r.corr)
    tau[i]<-term1
  }
  
  out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau,
            z_mean_2=z_mean_2, p.m=p.m,
            tau.1 = tau1,
            tau.2= p.m*z_mean_2 )
  
  #param2<<-out
  return(out)
}
SKAT_META_Optimal_Get_Pvalue<-function(Q.all, Phi, r.all, method){
  # no change
  n.r<-length(r.all)
  n.q<-dim(Q.all)[1]
  p.m<-dim(Phi)[2]
  
  lambda.all<-list()
  for(i in 1:n.r){
    r.corr<-r.all[i]
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Phi_rho<- L %*% (Phi %*% t(L))
    lambda.all[[i]]<-SKAT:::Get_Lambda(Phi_rho)
  }
  
  # Get Mixture param 
  param.m <- SKAT_META_Optimal_Param(Phi,r.all)
  Each_Info <- SKAT:::SKAT_Optimal_Each_Q(param.m, Q.all, r.all, lambda.all, method=method)
  pmin.q<-Each_Info$pmin.q
  pval <- rep(0,n.q)
  
  # added
  pmin<-Each_Info$pmin
  
  for(i in 1:n.q){
    pval[i]<-SKAT:::SKAT_Optimal_PValue_Davies(pmin.q[i,],param.m,r.all, pmin[i])
  }
  
  # Check the pval 
  # Since SKAT-O is between burden and SKAT, SKAT-O p-value should be <= min(p-values) * 2
  # To correct conservatively, we use min(p-values) * 3
  
  multi<-3
  if(length(r.all) < 3){
    multi<-2
  }
  
  for(i in 1:n.q){
    pval.each<-Each_Info$pval[i,]
    IDX<-which(pval.each > 0)
    
    pval1<-min(pval.each) * multi
    if(pval[i] <= 0 || length(IDX) < length(r.all)){
      pval[i]<-pval1
    }
    # if pval==0, use nonzero min each.pval as p-value
    if(pval[i] == 0){
      if(length(IDX) > 0){
        pval[i] = min(pval.each[IDX])
      }
    }
    
  }
  
  return(list(p.value=pval,p.val.each=Each_Info$pval))
}
SKAT_META_Optimal <- function(Score, Phi, r.all, method="davies"){
  # no change
  # if r.all >=0.999 ,then r.all = 0.999
  IDX<-which(r.all >= 0.999)
  if(length(IDX) > 0){
    r.all[IDX]<-0.999	
  }
  
  p.m<-dim(Phi)[2]
  n.r<-length(r.all)
  
  ###########################################
  # Compute Q.r and Q.r.res
  ##########################################
  out.Q <- SKAT_META_Optimal_Get_Q(Score, r.all)
  Q.res=NULL
  Q.all<-rbind(out.Q$Q.r, Q.res) 
  
  ##################################################
  # Compute P-values 
  #################################################
  
  out<-SKAT_META_Optimal_Get_Pvalue(Q.all, Phi/2, r.all, method)
  
  param<-list(p.val.each=NULL, q.val.each=NULL)
  param$p.val.each<-out$p.val.each[1,]
  param$q.val.each<-Q.all[1,]
  param$rho<-r.all
  param$minp<-min(param$p.val.each)
  
  id_temp<-which(param$p.val.each == min(param$p.val.each))
  id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
  if(length(id_temp1) > 0){
    param$rho[id_temp1] = 1
  }
  
  param$rho_est<-param$rho[id_temp]
  p.value<-out$p.value[1]
  re<-list(p.value = p.value, param=param)
  return(re)	
}

Met_SKAT_Get_Pvalue<-function(Score, Phi, r.corr, method){
  # change SKAT 
  Q.res = NULL
  p.m<-nrow(Phi)
  # if Phi==0
  if(sum(abs(Phi)) == 0){
    warning("No polymorphic SNPs!",call.=FALSE)
    return(list(p.value=1, p.value.resampling= NULL, pval.zero.msg=NULL))
  }
  
  if(length(Phi) <=1){
    r.corr=0
  } else{
    if(ncol(Phi) <=10){
      if(qr(Phi)$rank <= 1){
        r.corr=0
      }
    }
  }
  
  if(length(r.corr) > 1){
    re = SKAT_META_Optimal(Score, Phi, r.corr, method=method)
    return(re)
  } 
  
  if (r.corr == 0){
    Q<-sum(Score^2)/2
  } else if (r.corr==1){
    Q <- SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r
    a<- as.matrix(sum(Phi))
    re <- SKAT:::Get_Liu_PVal(Q, a, Q.res)
    return(re)
  } else {
    # like r.corr = 0.1 or 0.2
    Q <- SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Phi<- L %*% (Phi %*% t(L))
  }
  lambda <- SKAT:::Get_Lambda(Phi/2)
  re1 <- seqMeta:::pchisqsum2(Q, lambda, method = "integration", acc=1e-20)
  re <- list(p.value = re1$p, errflag = re1$errflag)
  # re<-SKAT:::Get_Davies_PVal(Q, Phi)
  return(re)
}


skatTest <- function(U, V, weight=NULL, r = 0){

  U =  U * weight
  V =  t(t(V * weight) * weight)
  Q.ker.optim.pval = Met_SKAT_Get_Pvalue(U, V, r, method="davies")$p.value
  return(list(Q.ker.optim.pval=Q.ker.optim.pval))

}
