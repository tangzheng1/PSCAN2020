# summary accuracy simulation results
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

power = NULL
snr.ls = NULL
causal.ls = NULL
dir.ls = NULL

for(i in 1:n.snr){
  
  for(j in 1:n.causal){
    
    
    for(k in 1:n.dir){

      string = paste("results_",shape, "_snr", snr[i], "_causal", causal[j], "_", eff.dir[k], sep="")
      if(causal[j]=="0.1"){
        eff.size = 0.3
      }else{
        eff.size = 0.15
      }
      string = paste(string, "_effsize", eff.size, ".Rdata", sep="")
      
      load(string)
      
      power = rbind(power, colMeans(combined.ret<1/1000000) )
      #colnames(results) = c("pscan.burden.pval", "burden.pval", "pscan.skat.pval", "skat.pval")
      #save(results, file=paste("results_", shape, "_", eff.dir, "_", snr[j], "_", eff.size, ".Rdata", sep=""))
      snr.ls = c(snr.ls, snr[i])
      causal.ls = c(causal.ls, causal[j])
      dir.ls = c(dir.ls, ifelse(eff.dir[k]=="unidir", "Unidirectional", "Bidirectional") )
      
    }
  }

}  

  
colnames(power) = c("PSCANM", "SCAN1DM", "Burden", "PSCANV", "SCAN1DV", "SKAT")
power.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, power)
save(power.rst, file="power.rst.Rdata")

