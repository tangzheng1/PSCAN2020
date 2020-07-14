shape = "circular"
#shape = "banded"

eff.dir = c("unidir", "bidir")
n.dir = length(eff.dir)
#eff.dir = "bidir"

snr = c("0.1", "0.25", "1")
#snr = 1
n.snr = length(snr)

causal = c("0.1", "0.5")
n.causal = length(causal)

# shape = "circular"
# eff.dir = "bidir"
# snr = 0.1
# causal = 0.1


snr.ls = NULL
causal.ls = NULL
dir.ls = NULL
sv.sen.avg = NULL; sv.spe.avg = NULL; sv.kappa.avg = NULL
pscan2.sen.avg = NULL; pscan2.spe.avg = NULL; pscan2.kappa.avg = NULL
pscan3.sen.avg = NULL;pscan3.spe.avg = NULL;pscan3.kappa.avg = NULL

n.job = 10

for(i in 1:n.snr){
  
  for(j in 1:n.causal){
    
    
    for(k in 1:n.dir){
      
      string = paste(shape, snr[i], "_", "causal", causal[j], "_", eff.dir[k], sep="")
      if(causal[j]=="0.1"){
        eff.size = 0.25
      }else{
        eff.size = 0.1
      }
      
      sv.sen.all = NULL; sv.spe.all = NULL; sv.kappa.all = NULL
      pscan2.sen.all = NULL; pscan2.spe.all = NULL; pscan2.kappa.all = NULL
      pscan3.sen.all = NULL; pscan3.spe.all = NULL; pscan3.kappa.all = NULL
      
      for(l in 1:n.job){
        
        tmp = try( load( paste(string, "/n=", l,  "/sv.sen_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          sv.sen.all = rbind(sv.sen.all, sv.sen)
        }else{
          #print(paste("===",l,"===="))
        }
        tmp = try( load( paste(string, "/n=", l,  "/sv.spe_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          sv.spe.all = rbind(sv.spe.all, sv.spe)
        }else{
          #print(paste("===",l,"===="))
        }
        tmp = try( load( paste(string, "/n=", l,  "/sv.kappa_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          sv.kappa.all = rbind(sv.kappa.all, sv.kappa)
        }else{
          #print(paste("===",l,"===="))
        }        
        
        
        tmp = try( load( paste(string, "/n=", l,  "/pscan2.sen_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          pscan2.sen.all = rbind(pscan2.sen.all, pscan2.sen)
        }else{
          #print(paste("===",l,"===="))
        }
        tmp = try( load( paste(string, "/n=", l,  "/pscan2.spe_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          pscan2.spe.all = rbind(pscan2.spe.all, pscan2.spe)
        }else{
          #print(paste("===",l,"===="))
        }
        tmp = try( load( paste(string, "/n=", l,  "/pscan2.kappa_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          pscan2.kappa.all = rbind(pscan2.kappa.all, pscan2.kappa)
        }else{
          #print(paste("===",l,"===="))
        }
        
        tmp = try( load( paste(string, "/n=", l,  "/pscan3.sen_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          pscan3.sen.all = rbind(pscan3.sen.all, pscan3.sen)
        }else{
          #print(paste("===",l,"===="))
        }
        tmp = try( load( paste(string, "/n=", l,  "/pscan3.spe_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          pscan3.spe.all = rbind(pscan3.spe.all, pscan3.spe)
        }else{
          #print(paste("===",l,"===="))
        }
        tmp = try( load( paste(string, "/n=", l,  "/pscan3.kappa_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          pscan3.kappa.all = rbind(pscan3.kappa.all, pscan3.kappa)
        }else{
          #print(paste("===",l,"===="))
        }
        

      }
      
      
      sv.sen.avg = rbind(sv.sen.avg, colMeans(sv.sen.all, na.rm=TRUE)); sv.spe.avg = rbind(sv.spe.avg, colMeans(sv.spe.all, na.rm=TRUE)); sv.kappa.avg = rbind(sv.kappa.avg, colMeans(sv.kappa.all, na.rm=TRUE))
      pscan2.sen.avg = rbind(pscan2.sen.avg, colMeans(pscan2.sen.all, na.rm=TRUE)); pscan2.spe.avg = rbind(pscan2.spe.avg, colMeans(pscan2.spe.all, na.rm=TRUE)); pscan2.kappa.avg = rbind(pscan2.kappa.avg, colMeans(pscan2.kappa.all, na.rm=TRUE))
      pscan3.sen.avg = rbind(pscan3.sen.avg, colMeans(pscan3.sen.all, na.rm=TRUE)); pscan3.spe.avg = rbind(pscan3.spe.avg, colMeans(pscan3.spe.all, na.rm=TRUE)); pscan3.kappa.avg = rbind(pscan3.kappa.avg, colMeans(pscan3.kappa.all, na.rm=TRUE))
      
      snr.ls = c(snr.ls, snr[i])
      causal.ls = c(causal.ls, causal[j])
      dir.ls = c(dir.ls, ifelse(eff.dir[k]=="unidir", "Unidirectional", "Bidirectional") )
    }
  }
}


colnames(sv.sen.avg) = c("FWER=0.01", "FWER=0.05")
sv.sen.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, sv.sen.avg)
colnames(pscan2.sen.avg) = c("FWER=0.01", "FWER=0.05")
pscan2.sen.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan2.sen.avg)
colnames(pscan3.sen.avg) = c("FWER=0.01", "FWER=0.05")
pscan3.sen.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan3.sen.avg)

colnames(sv.spe.avg) = c("FWER=0.01", "FWER=0.05")
sv.spe.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, sv.spe.avg)
colnames(pscan2.spe.avg) = c("FWER=0.01", "FWER=0.05")
pscan2.spe.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan2.spe.avg)
colnames(pscan3.spe.avg) = c("FWER=0.01", "FWER=0.05")
pscan3.spe.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan3.spe.avg)


colnames(sv.kappa.avg) = c("FWER=0.01", "FWER=0.05")
sv.kappa.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, sv.kappa.avg)
colnames(pscan2.kappa.avg) = c("FWER=0.01", "FWER=0.05")
pscan2.kappa.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan2.kappa.avg)
colnames(pscan3.kappa.avg) = c("FWER=0.01", "FWER=0.05")
pscan3.kappa.avg.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan3.kappa.avg)


save(sv.sen.avg.rst, file="sv.sen.avg.rst.Rdata")
save(sv.spe.avg.rst, file="sv.spe.avg.rst.Rdata")
save(sv.kappa.avg.rst, file="sv.kappa.avg.rst.Rdata")

save(pscan2.sen.avg.rst, file="pscan2.sen.avg.rst.Rdata")
save(pscan2.spe.avg.rst, file="pscan2.spe.avg.rst.Rdata")
save(pscan2.kappa.avg.rst, file="pscan2.kappa.avg.rst.Rdata")

save(pscan3.sen.avg.rst, file="pscan3.sen.avg.rst.Rdata")
save(pscan3.spe.avg.rst, file="pscan3.spe.avg.rst.Rdata")
save(pscan3.kappa.avg.rst, file="pscan3.kappa.avg.rst.Rdata")

####==== type I error ====####


shape = "circular"
eff.dir = "unidir"
snr = "1"
eff.size = 0


n.job = 100
sv.flag.all = NULL
pscan2.flag.all = NULL
pscan3.flag.all = NULL

for(l in 1:n.job){
  
  tmp = try( load( paste("./n=", l,  "/sv.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep="") ) )
  if(!inherits(tmp, 'try-error')){
    sv.flag.all = rbind(sv.flag.all, sv.flag)
  }else{
    print(paste("===",l,"===="))
  }
  
  tmp = try( load( paste("./n=", l,  "/pscan2.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep="") ) )
  if(!inherits(tmp, 'try-error')){
    pscan2.flag.all = rbind(pscan2.flag.all, pscan2.flag)
  }else{
    print(paste("===",l,"===="))
  }
  
  tmp = try( load( paste("./n=", l,  "/pscan3.flag_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep="") ) )
  if(!inherits(tmp, 'try-error')){
    pscan3.flag.all = rbind(pscan3.flag.all, pscan3.flag)
  }else{
    print(paste("===",l,"===="))
  }
  
  print(paste("===",l,"===="))

  
}

save(sv.flag.all, file="sv.flag.all.Rdata")
save(pscan2.flag.all, file="pscan2.flag.all.Rdata")
save(pscan3.flag.all, file="pscan3.flag.all.Rdata")




