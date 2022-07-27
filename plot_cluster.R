write.csv(B$Vsig,file = '/Users/maryammaghsoudi/Desktop/LUSCVsigINT.txt',quote = FALSE, row.names = F,col.names = NULL)

PositiveGroup <- c("ACC", "BLCA", "BRCA","GBMLGG", "HNSC", "KIPAN", "KIRC",
              "KIRP",  "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", 
              "THYM", "UCEC", "UVM")

stripchart(rep(1, length(group)) ~ group, vertical = FALSE, pch=16, col="white", 
           las=1, xlab='Proportion of significant random sets', xlim=c(0, 1), 
           main=paste('PhenoClust', 'clusters'))
#ratio of our 17 positive cancer cluster
rs <- NULL
group <- NULL
for (group_name in PositiveGroup){
  temp <- unlist(subRatios[group_name])
  if(all(is.na(temp))){
    rs <- c(rs, 2)
    group <- c(group, group_name)
  } else {
    #subcluster ratio
    rs <- c(rs, temp)
    group <- c(group, rep(group_name, length(temp)))
  }
}

stripchart(rs ~ group, vertical = FALSE, pch=16, col=1:length(PositiveGroup), add=TRUE)
for (i in 1:length(PositiveGroup)){
  abline(h=i, col='lightgrey', lty=3)
  j=PositiveGroup[i]
  lines(y=c(i-0.4, i+0.4), x=c(ratios[j, 7], ratios[j, 7]), col=8, lwd=2)
}
  # lines(y=c(i-0.4, i+0.4), x=c(ratios[i, 7], ratios[i, 7]), col=i, lwd=2)

abline(v=0.05, col='grey', lwd=1, lty=2)



