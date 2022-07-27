find_pval_dat=function(D,RAND){
  pval_dat_rand=c()
  dat=remove_sig(D)
  for (j in 1:nrow(RAND)) {
    random_g=RAND[j,]
    dat_rand=dat[,random_g]
    pval_dat_rand[j]=findpval(dat_rand)
  }
return(pval_dat_rand)
  }