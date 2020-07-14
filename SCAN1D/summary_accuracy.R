# summary accuracy simulation results
shape = "circular"
#shape = "banded"

eff.dir = c("unidir", "bidir")
n.dir = length(eff.dir)
#eff.dir = "bidir"

snr = c("0.1", "0.25", "1")
#snr = 1
n.snr = length(snr)

#causal="0.1"
causal = c("0.1", "0.5")
n.causal = length(causal)

sv.sen = NULL; sv.spe = NULL
scan1DM.sen = NULL; scan1DM.spe = NULL
scan1DV.sen = NULL; scan1DV.spe = NULL
pscan2.sen = NULL; pscan2.spe = NULL
pscan3.sen = NULL; pscan3.spe = NULL

snr.ls = NULL
causal.ls = NULL
dir.ls = NULL

for(i in 1:n.snr){
  
  for(j in 1:n.causal){
    
    
    for(k in 1:n.dir){

      string = paste("results_",shape, "_snr", snr[i], "_causal", causal[j], "_", eff.dir[k], sep="")
      if(causal[j]=="0.1"){
        #eff.size = 0.3
        #eff.size = 0.15
        eff.size = 0.25
      }else{
        eff.size = 0.1
      }
      string = paste(string, "_effsize", eff.size, ".Rdata", sep="")
      
      load(string)
      
      tmp = colMeans(combined.ret)
      
      sv.sen = rbind(sv.sen, tmp[c(1,2)] );  sv.spe = rbind(sv.spe, tmp[c(3,4)] )
      scan1DM.sen = rbind(scan1DM.sen, tmp[c(5,6)] );  scan1DM.spe = rbind(scan1DM.spe, tmp[c(7,8)] )
      scan1DV.sen = rbind(scan1DV.sen, tmp[c(9,10)] );  scan1DV.spe = rbind(scan1DV.spe, tmp[c(11,12)] )
      pscan2.sen = rbind(pscan2.sen, tmp[c(13,14)] );  pscan2.spe = rbind(pscan2.spe, tmp[c(15,16)] )
      pscan3.sen = rbind(pscan3.sen, tmp[c(17,18)] );  pscan3.spe = rbind(pscan3.spe, tmp[c(19,20)] )
      #colnames(results) = c("pscan.burden.pval", "burden.pval", "pscan.skat.pval", "skat.pval")
      #save(results, file=paste("results_", shape, "_", eff.dir, "_", snr[j], "_", eff.size, ".Rdata", sep=""))
      snr.ls = c(snr.ls, snr[i])
      causal.ls = c(causal.ls, causal[j])
      dir.ls = c(dir.ls, ifelse(eff.dir[k]=="unidir", "Unidirectional", "Bidirectional") )
      
    }
  }

}  

  
colnames(sv.sen) = c("FWER.0.01", "FWER.0.05")
sv.sen.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, sv.sen)
save(sv.sen.rst, file="sv.sen.rst.Rdata")

colnames(sv.spe) = c("FWER.0.01", "FWER.0.05")
sv.spe.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, sv.spe)
save(sv.spe.rst, file="sv.spe.rst.Rdata")


colnames(pscan2.sen) = c("FWER.0.01", "FWER.0.05")
pscan2.sen.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan2.sen)
save(pscan2.sen.rst, file="pscan2.sen.rst.Rdata")

colnames(pscan2.spe) = c("FWER.0.01", "FWER.0.05")
pscan2.spe.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan2.spe)
save(pscan2.spe.rst, file="pscan2.spe.rst.Rdata")


colnames(pscan3.sen) = c("FWER.0.01", "FWER.0.05")
pscan3.sen.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan3.sen)
save(pscan3.sen.rst, file="pscan3.sen.rst.Rdata")

colnames(pscan3.spe) = c("FWER.0.01", "FWER.0.05")
pscan3.spe.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, pscan3.spe)
save(pscan3.spe.rst, file="pscan3.spe.rst.Rdata")

colnames(scan1DM.sen) = c("FWER.0.01", "FWER.0.05")
scan1DM.sen.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, scan1DM.sen)
save(scan1DM.sen.rst, file="scan1DM.sen.rst.Rdata")

colnames(scan1DM.spe) = c("FWER.0.01", "FWER.0.05")
scan1DM.spe.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, scan1DM.spe)
save(scan1DM.spe.rst, file="scan1DM.spe.rst.Rdata")


colnames(scan1DV.sen) = c("FWER.0.01", "FWER.0.05")
scan1DV.sen.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, scan1DV.sen)
save(scan1DV.sen.rst, file="scan1DV.sen.rst.Rdata")

colnames(scan1DV.spe) = c("FWER.0.01", "FWER.0.05")
scan1DV.spe.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, scan1DV.spe)
save(scan1DV.spe.rst, file="scan1DV.spe.rst.Rdata")


