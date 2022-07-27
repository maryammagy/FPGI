findpval=function(mat){     
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
}