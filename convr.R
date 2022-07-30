library(survival)
library(Biobase)
library(BiocGenerics)
library(parallel)
library(samr)
library(siggenes)

source('cleanData.r')
source('findDEG.r')
source('findpval.r')
source('remove_sig.r')
source('find_pval_dat.r')
source('remove_meta_PCNA.r')

dataset_dir = '../TCGA/'
PositiveGroup <- c("ACC", "BLCA", "BRCA","GBMLGG", "HNSC", "KIPAN", "KIRC",
                   "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", 
                   "THYM", "UCEC", "UVM")

for (d in PositiveGroup){
   # this may take a while
  load(paste0(dataset_dir,d,'.rData'))
  dataObj = cleanData(dataObj)
  q = 0.5
  SizeS = 50
  status = dataObj$Surv[,2]
  surv = dataObj$Surv
  DEG = list()
  number = c()
  DEGFinal = list()
  PVAL = c()
  INT = c()
  V = c()
  Vsig = c()
  RAND = c()
  pval_dat_PCNA_rand = c()
  pval_rand = c()
  M = 1
  for (k in 1:1000){
    message(paste(k))
    randgene = sample(colnames(dataObj$merged.dat), size =SizeS, replace = FALSE )
    RAND = rbind(RAND,randgene)
    mat_rand = dataObj$merged.dat[,randgene]
    pval_rand[k] = findpval(mat_rand)
    if (pval_rand[k] < 0.05){
      DEG[[1]] = findDEG(randgene)
      i = 2
      int = c()
      int[1] = length(intersect(randgene,DEG[[1]][[1]]))
      while (int[i-1] < SizeS) {
        DEG[[i]] = findDEG(DEG[[i-1]][[1]])
        int[i] = length(intersect(DEG[[i]][[1]],DEG[[i-1]][[1]]))
        i = i+1  
        }
      message(paste(int))
      INT = rbind(INT,int)
      number[M] = i
      DEGFinal[[M]] = DEG[[i-1]][[1]]
      PVAL[M] = DEG[[i-1]][[2]]
      V = c(V,DEGFinal[[M]])
      if (DEG[[i-1]][[2]] < 0.05){
        Vsig = c(Vsig,DEGFinal[[M]])
        }
      M = M+1
      print(M)
      }}


write.table(sort(table(V), decreasing = TRUE), file = paste0('TABLEV', d,'.txt'), quote = FALSE, row.names = F)
write.table(sort(table(Vsig), decreasing = TRUE), file =paste0('TABLEVSIG', d,'.txt'), quote = FALSE, row.names = F)
write.table(number, file = paste0('int', d,'.txt'), quote = FALSE, row.names = F)
write.table(PVAL,file = paste0('PVAL', d,'.txt'), quote = FALSE, row.names = F)
write.table(intersect(randgene,V), file = paste0('intersect', d,'.txt'), quote = FALSE, row.names = F)
}
