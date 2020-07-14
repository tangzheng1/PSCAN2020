library(data.table)
library(stringr)
library(caret)
source("source.R")


args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

# seed = 1000
# n.sim = 100
# shape = "circular"
# eff.dir = "unidir"
# snr = 1
# causal = 0  ## 0%, 10% or 50% causal
# eff.size = 0
#nperm = 10000

if(causal==0){
  nperm = 5000
}else{
  nperm = 1000
}



n = 5000 

############ simulation start ############
load("haplo_data.Rdata")
load("haplo_pos.Rdata")
pos.cand = haplo_pos$CHROM_POS[-which(haplo_pos$CHROM_POS > (999999-3000))]

fwer.grid = c(0.01, 0.05)
#fwer.grid = 0.05
n.fwer = length(fwer.grid)
#n.sim = 10
pscan2.sen = matrix(NA, n.sim, n.fwer)
pscan2.spe = matrix(NA, n.sim, n.fwer)
pscan2.kappa = matrix(NA, n.sim, n.fwer)
pscan3.sen = matrix(NA, n.sim, n.fwer)
pscan3.spe = matrix(NA, n.sim, n.fwer)
pscan3.kappa = matrix(NA, n.sim, n.fwer)
sv.sen = matrix(NA, n.sim, n.fwer)
sv.spe = matrix(NA, n.sim, n.fwer)
sv.kappa = matrix(NA, n.sim, n.fwer)

sv.flag = NULL
pscan2.flag = NULL
pscan3.flag = NULL

for(l in 1:n.sim){
  
  set.seed(seed + l)
  
  ####==== genotypes ====#### 
  start_pos <- sort(sample(pos.cand, 1))
  end_pos <- start_pos + 3000
  
  windows.subdat <- haplo_pos[haplo_pos$CHROM_POS >= start_pos & haplo_pos$CHROM_POS <= end_pos,]
  SNP_id <- windows.subdat[, 1]
  SNP.subset <- haplo_data[, which(colnames(haplo_data) %in% paste0("SNP", SNP_id))]
  SNP.subset <- ifelse(SNP.subset == 2, 0, 1)
  sub.ind <- sample(1:nrow(haplo_data), n) 
  SNP.subset1 <- SNP.subset[sub.ind, ]
  SNP.subset2 <- SNP.subset[-sub.ind, ]
  G <- SNP.subset1 + SNP.subset2  
  
  ## flip if allele if maf>0.5
  maf = colMeans(G)/2
  idx = which(maf>0.5)
  for(i in idx){
    G[which(G[,i]==2),i] = 3
    G[which(G[,i]==0),i] = 2
    G[which(G[,i]==3),i] = 0
  }
  
  ## trim SNPs that are highly corrlated
  hc = findCorrelation(cor(G), cutoff=0.9) 
  hc = sort(hc)
  if(length(hc)>0){
    G = G[,-hc]
  }
  
  ## combine SNPs that have mac<10
  mac = colSums(G)
  G.tmp = NULL
  G.burden = rep(0, n)
  flag = 0
  for(j in 1:ncol(G)){
    if(mac[j]<10){
      flag = 1
      G.burden = G.burden + G[,j]
    }else{
      G.tmp = cbind(G.tmp, G[,j])
    }
    
  }
  
  if(flag){
    G.tmp = cbind(G.burden, G.tmp)
  }
  
  G = G.tmp
  if(ncol(G)==1){
    next
  }
  maf = colMeans(G)/2
  
  
  ####==== set-up genetic effects====#### 
  m = ncol(G)  # total SNP
  m.c = as.integer(causal * m) # causal SNP
  if(causal>0 & m.c==0){
    m.c = 1  
  }
  m.nc = m - m.c # neutral SNP
  
  
  par.beta = rep(0,m)
  
  if(m.c>0){
    causal.idx = sample(1:m, m.c)
    
    par.beta.causal = eff.size * abs(log10(maf[causal.idx])) 
    
    if(eff.dir == "bidir"){
      #m.cp = rbinom(1,m.c,0.5) # positive SNP
      m.cp =  as.integer(m.c/2)
      neg.idx = sample(1:m.c, m.c-m.cp)
      par.beta.causal[neg.idx] = -par.beta.causal[neg.idx]
    }
    
    par.beta[causal.idx] = par.beta.causal
  }
  
  
  
  ####==== sample SNP coordinates ====#### 
  all.cords = matrix(rnorm(m*3), nrow=m, ncol=3)
  
  if(m.c>0){
    
    if(shape=="circular"){
      
      all.cords[causal.idx,1] = rnorm(m.c, 0, snr)
      all.cords[causal.idx,2] = rnorm(m.c, 0, snr)
      all.cords[causal.idx,3] = rnorm(m.c, 0, snr)
      
    }else if(shape=="banded"){
      
      all.cords[causal.idx,1] = rnorm(m.c, 0, snr)
      all.cords[causal.idx,2] = rnorm(m.c, 0, 0.1)
      all.cords[causal.idx,3] = rnorm(m.c, 0, 0.1)
      
    }else{
      
      stop("Unknown shape.")
    }    
    
  }
  
  cords = list()
  cords = c(cords, list( all.cords ) )
  
  
  ####==== phenotypes ====#### 
  Z = rnorm(n)
  Y = 0.3*Z + as.numeric( G %*% par.beta ) + rnorm(n)
  
  ####==== summary stats ====#### 
  weight.burden = 1/sqrt(maf*(1-maf))
  weight.skat = dbeta(maf, 1, 25)
  
  mac = colSums(G)
  
  
  
  # generate U and V  
  covariates = cbind(1, Z)
  a.reduce = summary(glm(Y ~ Z))
  gamma.hat = c(a.reduce$coefficient[1,1], a.reduce$coefficient[2,1])
  covariates_gamma.hat = as.numeric(covariates%*%matrix(gamma.hat,nrow=2))
  sigma.sq = sum((Y -covariates_gamma.hat)^2)/n
  
  # U
  res.y.reduce = Y - covariates_gamma.hat
  U = as.numeric( res.y.reduce%*%(G)/sigma.sq )
  
  # V
  GX = t(G) %*% covariates
  V = (t(G) %*% G - GX %*% solve(t(covariates)%*%covariates) %*% t(GX))/sigma.sq 
  
  ####==== signal region detection ====####
  if(m.c>0){
    neutral.idx = setdiff( c(1:m), causal.idx)
  }else{
    neutral.idx = c(1:m)
  }
  
  
  #### SV ####
  sv.flag.tmp = rep(0, n.fwer)
  cluster.list = as.list(1:m)
  
  n.scan = length(cluster.list)
  C = matrix(0, nrow=n.scan, ncol=m)
  for(j in 1:n.scan){
    C[j, cluster.list[[j]]] = 1
  }
  
  UC = C %*% U
  VC = C %*% V %*% t(C)
  Q = UC*UC / diag(VC)
  
  Q.max.perm = NULL
  set.seed(11)
  for(s in 1:nperm){
    U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    tmp = C %*% U.perm
    Q.max.perm = c(Q.max.perm, max( tmp*tmp/diag(VC) ) )
  }
  
  for(i in 1:n.fwer){
    
    signal.alpha = fwer.grid[i]
    signal.thresh = quantile(Q.max.perm, 1-signal.alpha)
    cand.idx = which(Q >= signal.thresh)
    
    signal.region = NULL
    if( length(cand.idx)>0 ){
      signal.region = which(Q >= signal.thresh)
    }
    
    if( length(signal.region)>0 ){
      sv.flag.tmp[i]=1
    }
    
  }
  
  sv.flag = rbind(sv.flag, sv.flag.tmp)

  
  #### P-SCAN: define windows ####
  
  cluster.list = NULL
  cords.ls = cords
  
  m.idx = 0
  cords.dist = as.matrix(dist(cords.ls[[1]]))
  dist.uniq = sort(unique(  as.numeric(cords.dist[lower.tri(cords.dist)])  ))
  n.poss = length(dist.uniq)
  cluster.list.d = as.list(c(1:m)) 
  cluster.pre = cluster.list.d
  cluster.num.pre = m
  
  for(j in 1:n.poss){
    tmp = clusters(graph_from_adjacency_matrix(0+(cords.dist<=dist.uniq[j]), mode="undirected"))
    cluster.num.cur = tmp$no
    
    if(cluster.num.cur == cluster.num.pre){
      next
    }else{
      cluster.cur = split(1:m, tmp$membership)
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
  
  n.scan = length(cluster.list)
  C = matrix(0, nrow=n.scan, ncol=m)
  for(j in 1:n.scan){
    C[j, cluster.list[[j]]] = 1
  }
  
  C.base = C
  
  #### P-SCAN-M ####
  pscan2.flag.tmp = rep(0, n.fwer)
  C = t(t(C.base) * weight.burden)
  
  UC = C %*% U
  VC = C %*% V %*% t(C)
  Q = UC*UC / diag(VC)
  Q.max =  max(Q)
  
  set.seed(11)
  Q.max.perm = NULL
  for(s in 1:nperm){
    U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    tmp = C %*% U.perm
    Q.max.perm = c(Q.max.perm, max( tmp*tmp/diag(VC) ) )
  }
  
  for(i in 1:n.fwer){
    
    signal.alpha = fwer.grid[i]
    signal.thresh = quantile(Q.max.perm, 1-signal.alpha)
    cand.idx = which(Q >= signal.thresh)
    
    signal.region = NULL
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
    
    if( length(signal.region)>0 ){
      pscan2.flag.tmp[i]=1
    }
    
  }
  
  pscan2.flag = rbind(pscan2.flag, pscan2.flag.tmp)
  
  #### P-SCAN-V ####
  pscan3.flag.tmp = rep(0, n.fwer)
  C = t(t(C.base) * weight.skat)

  Qp = NULL
  for(j in 1:n.scan){
    idx = which(C[j,]!=0)
    Qp = c(Qp, skatTest(U[idx], V[idx,idx], weight=C[j,idx], r = 0)$Q.ker.optim.pval )
  }
  Qp[Qp<0] = 0  ## occasionally skatTest give negative p-val
  Qp.min = min(Qp)

  set.seed(11)
  Qp.min.perm = NULL
  for(s in 1:nperm){
    U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

    Qp.perm = NULL
    for(j in 1:n.scan){
      idx = which(C[j,]!=0)
      Qp.perm = c(Qp.perm, skatTest(U.perm[idx], V[idx,idx], weight=C[j,idx], r = 0)$Q.ker.optim.pval )
    }
    Qp.perm[Qp.perm<0] = 0  ## occasionally skatTest give negative p-val
    Qp.min.perm = c(Qp.min.perm, min(Qp.perm) )
  }


  for(i in 1:n.fwer){

    signal.alpha = fwer.grid[i]
    signal.thresh = quantile(Qp.min.perm, signal.alpha)
    cand.idx = which(Qp <= signal.thresh)

    signal.region = NULL
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

    if( length(signal.region)>0 ){
      pscan3.flag.tmp[i]=1
    }

  }

  pscan3.flag = rbind(pscan3.flag, pscan3.flag.tmp)
  
  print(l)
  

  save(sv.flag, file=paste("sv.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep=""))
  save(pscan2.flag, file=paste("pscan2.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep=""))
  save(pscan3.flag, file=paste("pscan3.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep=""))
  
  
}



save(sv.flag, file=paste("sv.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep=""))
save(pscan2.flag, file=paste("pscan2.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep=""))
save(pscan3.flag, file=paste("pscan3.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep=""))

