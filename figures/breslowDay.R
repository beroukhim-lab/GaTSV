#code source:https://online.stat.psu.edu/onlinecourses/sites/stat504/files/lesson04/breslowday.test_.R
######################################################################
# Function to perform the Breslow and Day (1980) test including
# the corrected test by Tarone
# Uses the equations in Lachin (2000) p. 124-125.
#
# Programmed by Michael Hoehle <http://www-m4.ma.tum.de/pers/hoehle>
# Note that the results of the Tarone corrected test do
# not correspond to the numbers in the Lachin book...
#
# Params:
#  x - a 2x2xK contingency table
#
# Returns:
#  a vector with three values
#   1st value is the Breslow and Day test statistic
#   2nd value is the correct test by Tarone
#   3rd value - p value based on the Tarone test statistic
#               using a \chi^2(K-1) distribution
######################################################################

breslowday.test <- function(x) {
  #Find the common OR based on Mantel-Haenszel
  or.hat.mh <- my_mantelhaen.test(x)$estimate #modified this to refer to my modification that's receptivr to large entries
  #Number of strata
  K <- dim(x)[3]
  #Value of the Statistic
  X2.HBD <- 0
  #Value of aj, tildeaj and Var.aj
  a <- tildea <- Var.a <- numeric(K)
  
  for (j in 1:K) {
    #Find marginals of table j
    mj <- apply(x[,,j], MARGIN=1, sum)
    nj <- apply(x[,,j], MARGIN=2, sum)
    
    #Solve for tilde(a)_j;#modified these sections again to allow for large value computations
    coef <- c(-(as.numeric(mj[1])*as.numeric(nj[1])) * or.hat.mh, as.numeric(nj[2])-as.numeric(mj[1])+
                or.hat.mh*(as.numeric(nj[1])+as.numeric(mj[1])),
              1-or.hat.mh)
    sols <- Re(polyroot(coef))
    #Take the root, which fulfills 0 < tilde(a)_j <= min(n1_j, m1_j)
    tildeaj <- sols[(0 < sols) &  (sols <= min(nj[1],mj[1]))]
    #Observed value
    aj <- x[1,1,j]
    
    #Determine other expected cell entries
    tildebj <- mj[1] - tildeaj
    tildecj <- nj[1] - tildeaj
    tildedj <- mj[2] - tildecj
    
    #Compute \hat{\Var}(a_j | \widehat{\OR}_MH)
    Var.aj <- (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)^(-1)
    
    #Compute contribution
    X2.HBD <- X2.HBD + as.numeric((aj - tildeaj)^2 / Var.aj)
    
    #Assign found value for later computations
    a[j] <- aj ;  tildea[j] <- tildeaj ; Var.a[j] <- Var.aj
  }
  
  #Compute Tarone corrected test
  X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.aj) )
  
  #Compute p-value based on the Tarone corrected test
  p <- 1-pchisq(X2.HBDT, df=K-1)
  
  res <- list(X2.HBD=X2.HBD,X2.HBDT=X2.HBDT,p=p)
  class(res) <- "bdtest"
  return(res)
}

print.bdtest <- function(x) {
  cat("Breslow and Day test (with Tarone correction):\n")
  cat("Breslow-Day X-squared         =",x$X2.HBD,"\n")
  cat("Breslow-Day-Tarone X-squared  =",x$X2.HBDT,"\n\n")
  cat("Test for test of a common OR: p-value = ",x$p,"\n\n")
}
