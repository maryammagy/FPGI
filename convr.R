setwd('/Users/maryammaghsoudi/R/maryam-random/BRCA/convergent/TCGA')
load('/Users/maryammaghsoudi/R/yishai/datasets/BRCA.rData')

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

dataObj = cleanData(dataObj)
q = 0.5
SizeS = 50
status = dataObj$Surv[,2]
surv = dataObj$Surv
DEG = list()
number = -c()
DEGFinal = list()
PVAL = -c()
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

#table(), performs categorical tabulation of data with the variable and its frequency.
write.table(sort(table(V), decreasing = TRUE), file = "Result/TABLEV.txt", quote = FALSE, row.names = F)
write.table(sort(table(Vsig), decreasing = TRUE), file = "Result/TABLEVSIG.txt", quote = FALSE, row.names = F)
write.table(number, file = "Result/int.txt", quote = FALSE, row.names = F)
write.table(PVAL, file = "Result/PVAL.txt", quote = FALSE, row.names = F)
A=read.table(file = "Result/TABLEV.txt", header = TRUE)
B=read.table(file = "Result/TABLEVSIG.txt", header = TRUE)
#C=A$V[which(A$Freq>=max(A$Freq))]
#D=B$V[which(B$Freq>=max(B$Freq))]

write.table(intersect(randgene,V), file = "Result/intersect.txt", quote = FALSE, row.names = F)
#write.table(C, file = "Result/important.txt", quote = FALSE, row.names = F)
#write.table(D, file = "Result/important2.txt", quote = FALSE, row.names = F)
#write.table(c(length(V), length(table(V)),length(table(Vsig)),length(table(V))/ length(V),length(table(Vsig))/ length(V),length(C),length(D),max(A$Freq),max(B$Freq)),file = "Result/Unique.txt", quote = FALSE, row.names = F)

# remove meta_PCNA signature and effect of that in P-values
dat_PCNA=remove_meta_PCNA(dataObj)
for (n in 1:nrow(RAND)) {
  random_gene=RAND[n,]
  dat_PCNA_rand=dat_PCNA[,random_gene]
  pval_dat_PCNA_rand[n]=findpval(dat_PCNA_rand)
}

result = cbind(pval_rand,pval_dat_PCNA_rand) 

t=c()
all_Freq=B$V[which(B$Freq>=1)]
t[1]=length(all_Freq)
pval_all_Freq=find_pval_dat(all_Freq,RAND)
result=cbind(result,pval_all_Freq)            

for (i in 2:1000) {
  p=i*0.001
  quan=quantile(B$Freq,probs = p)
  D=B$V[which(B$Freq>=quan)]
  t[i]=length(D)
  if(t[i]!=t[i-1]){
    pval_dat_rand=find_pval_dat(D,RAND)
    result=cbind(result,pval_dat_rand)
  }else{
    i=i+1}
}
 small_p=c()
for (j in 1:ncol(result)) {
  small_p[j]=sum(result[,j]<0.05)
}
t=unique(t)
tab=c()
tab=rbind(length(status),ncol(dataObj$merged.dat),length(PVAL),(sum(PVAL<0.05)),
          nrow(B),max(B$Freq),sum(pval_dat_PCNA_rand<0.05), min(small_p)
          ,
          t[which(small_p==min(small_p))-2], B$Freq[t[which(small_p==min(small_p))-2]]
)
rownames(tab)=c('sample','gene','pval rand < 0.05','pval DEGfinal < 0.05',
                'candidated genes','max Freq','venet correction'
             ,'my correction'
                ,
                "significant gene","MAX Freq significant genes"
                )
tab=as.data.frame(tab)  
colnames(tab)=c('UCEC1000')
write.xlsx(tab,row.names=TRUE,file = '/Users/maryammaghsoudi/R/maryam-random/UCEC/UCEC1000.xlsx')

VsigINT=Reduce(intersect, list(ucec1000$Vsig,ucec2000$Vsig))
intersect(ucec1000$Vsig[1:50],ucec2000$Vsig[1:50])
lusc1000=B

lusc2000=B


brca2000=B
brca2000_2=B

save(VsigINT,file = '/Users/maryammaghsoudi/R/maryam-random/VsigINT/UCECvSIGint.rData')
