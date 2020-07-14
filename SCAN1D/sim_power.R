library(data.table)
library(stringr)
library(caret)
source("source.R")



args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

# perl submit_power.pl 10 100 circular unidir 0.1 0.1 0.3
# job = 1 # 1:10  parallele job
# seed = 1000000 + job
# # seed = 1000
# n.sim = 100
# shape = "circular"
# eff.dir = "unidir"
# snr = 1
# causal = 0.1  ## 0%, 10% or 50% causal
# eff.size = 0.3


n = 5000



############ simulation start ############
load("haplo_data.Rdata")
load("haplo_pos.Rdata")

pos.cand = haplo_pos$CHROM_POS[-which(haplo_pos$CHROM_POS > (999999-3000))]

results = NULL

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

  POS <- windows.subdat$CHROM_POS
  
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
    POS = POS[-hc]
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
  

  ####==== perform test ====####
  
  ## P-SCAN-M
  pscan.burden =  scanTest2(U, V, mac, weight=weight.burden, cords, signal.alpha=NULL)$Q.scan2.pval

  ## P-SCAN-V
  pscan.skat = scanTest3(U, V, mac, weight=weight.skat, cords, signal.alpha=NULL)$Q.scan3.pval
  
  # SCAN-M: 1D scan test
  scan1D.burden = scanTest2.1D(U, V, mac, weight=weight.burden, POS, signal.alpha=NULL)$Q.scan2.pval
  
  # SCAN-V: 1D scan test
  scan1D.skat = scanTest3.1D(U, V, mac, weight=weight.skat, POS, signal.alpha=NULL)$Q.scan3.pval
  
  ## existing methods
  burden.U = weight.burden %*% U
  burden.V = weight.burden %*% V %*% weight.burden
  burden =  1-pchisq( (burden.U^2/burden.V) ,1)
 
  skat = skatTest(U, V, weight.skat, r=0)$Q.ker.optim.pval
  
  results = rbind(results, c(pscan.burden, scan1D.burden, burden, pscan.skat, scan1D.skat, skat))
  print(l)
  
  
}


colnames(results) = c("pscan_burden", "scan1D_burden", "burden", "pscan_skat", "scan1D_skat", "skat")
save(results, file=paste("results_", shape, snr, "_causal", causal, "_", eff.dir, "_", eff.size,".Rdata", sep=""))

