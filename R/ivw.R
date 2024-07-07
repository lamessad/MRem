IVW= function(betaY, betaX, betaYse) {

  beta=0;se=1;p=1
  if (length(betaY) > 0) {
    wv = 1/betaYse^2
    mod1 = lm(betaY ~ -1 + betaX, weights=wv)
    sod1 = summary(mod1)
    df   = length(betaY) - 1
    beta = sod1$coefficients[1,1]
    #se   = sod1$coefficients[1,2]/sod1$sigma
    se   = sod1$coefficients[1,2]/min(1,sod1$sigma) #check
    p    = 2*(pt(abs(beta/se),df, lower.tail=F))
  }
  z    = c(beta,se,beta/se,p)
  names(z)= c("Estimate","SE","t value","P")
  return(z)
}
