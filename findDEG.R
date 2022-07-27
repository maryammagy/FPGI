findDEG=function(geness){     
  mat=dataObj$merged.dat[,geness]
  pc <- prcomp(mat)
  pred <- pc$x[ ,1]

  group <- pred>quantile(pred, probs = q)
  group2 <- pred>quantile(pred, probs = 1-q)
  if (wilcox.test(surv[,1]~group)$statistic > wilcox.test(surv[,1]~group2)$statistic)
    formula <- surv~group
  else
    formula <- surv~group2
  sf <- survfit(formula)
  
  pVal <- pchisq(q=survdiff(formula)$chisq,df=1,lower.tail=F) 

  no_gp=which(group==FALSE)
  no_gp2=which(group==TRUE)
  #print(no_gp)
  gene_gp=dataObj$merged.dat[no_gp,]
  gene_gp2=dataObj$merged.dat[no_gp2,]
  allgenes=rbind.data.frame(gene_gp,gene_gp2)
  clabel1=rep(1,length(no_gp))
  clabel2=numeric(length(no_gp2))
  cl=c(clabel1,clabel2)
  data=t(allgenes)
  d.stat(data,cl, var.equal = FALSE, B = 100,R.fold = 1,na.method = "mean")
  sam.out<-sam(data, cl, method = d.stat)
  Delta<-findDelta(sam.out, fdr = NULL, genes =50, prec = 4, initial = NULL, verbose = FALSE)
  
  sam.sum <- summary(sam.out, Delta[1], entrez=FALSE)
  d<-attributes(sam.sum)$mat.sig
  a<-row.names(d)
  
  deg=list(a,pVal)
  return(deg)
}
