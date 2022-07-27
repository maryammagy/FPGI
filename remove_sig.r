remove_sig <- function(C){
  imp_data <- intersect(C, colnames(dataObj$merged.dat))
  imp_sig <- dataObj$merged.dat[ , imp_data]

  imp_sig <- apply(X = imp_sig, MARGIN = 1, FUN = function(x) median(x, na.rm = T))
  
  dat <- apply(X = dataObj$merged.dat, MARGIN = 2, FUN = function(x, imp_sig){
    lm(x~imp_sig)$residuals + mean(x)
  }, imp_sig)
  return(dat)
}
