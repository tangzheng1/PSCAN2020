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

n.job = 10

power = NULL
snr.ls = NULL
causal.ls = NULL
dir.ls = NULL

for(i in 1:n.snr){
  
  for(j in 1:n.causal){
    
    
    for(k in 1:n.dir){

      string = paste(shape, snr[i], "_", "causal", causal[j], "_", eff.dir[k], sep="")
      if(causal[j]=="0.1"){
        eff.size = 0.3
      }else{
        eff.size = 0.15
      }
      
      results.all = NULL
      
      for(l in 1:n.job){
        
        tmp = try( load( paste(string, "/n=", l,  "/results_", shape, "_", eff.dir[k], "_", snr[i], "_", eff.size, ".Rdata", sep="") ) )
        if(!inherits(tmp, 'try-error')){
          results.all = rbind(results.all, results)
        }else{
          print(paste("===",l,"===="))
        }
        
        # if(nrow(results)!=1000){
        #   print(paste("===",j,"====",l))
        # }
        
      }
      
      power = rbind(power, colMeans(results.all<1/1000000) )
      #colnames(results) = c("pscan.burden.pval", "burden.pval", "pscan.skat.pval", "skat.pval")
      #save(results, file=paste("results_", shape, "_", eff.dir, "_", snr[j], "_", eff.size, ".Rdata", sep=""))
      snr.ls = c(snr.ls, snr[i])
      causal.ls = c(causal.ls, causal[j])
      dir.ls = c(dir.ls, ifelse(eff.dir[k]=="unidir", "Unidirectional", "Bidirectional") )
      
    }
  }

}  

  
colnames(power) = c("PSCANM", "Burden", "PSCANV", "SKAT")
power.rst = data.frame(causal=causal.ls, dir=dir.ls, snr=snr.ls, power)
save(power.rst, file="power.rst.Rdata")

####==== type I error ====####


shape = "circular"
eff.dir = "unidir"
snr = "1"
eff.size = 0


n.job = 100
results.all = NULL
for(l in 1:n.job){
  
  tmp = try( load( paste("./n=", l,  "/results_", shape, "_", eff.dir, "_", snr, "_", eff.size, ".Rdata", sep="") ) )
  if(!inherits(tmp, 'try-error')){
    results.all = rbind(results.all, results)
  }else{
    print(paste("===",l,"===="))
  }
  
  print(paste("===",l,"===="))
  # if(nrow(results)!=1000){
  #   print(paste("===",j,"====",l))
  # }
  
}
 
save(results.all, file="results.all.Rdata")
