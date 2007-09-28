# power calculation for studies using baseline measure 
#   - simulation: continuous response, baseline and genotype as covariate (ANCOVA)
#   - can specify various modes of inheritance ('additive', 'dominant','recessive')
#   - use compound symmetry for covariance matrix 
# Author: Michael Man
#   Date: June 22, 2004
#      N: total number of subjects
#      p: frequency of A allele
#  alpha: significance level                     (used only in 'power.genotype.conti')
#    Rep: number of iterations to generate power (used only in 'power.genotype.conti')
#     pi: correlation coefficient
#    me1: mean of control group
#  delta: treatment/genotype effect
#  sd1/2: standard deviation of the control and treatment groups
#   minh: mode of inheritance
# genotype.delta: the effect due to individual genotype effect or overall effect
# Factor: whether treat 'Trt' as a factor in 'x'
#   TEST: debug
# reference:
#   Frison and Pocock (1992) "Repeated measures in clinical trials: analysis using mean summary statistics
#     and its implications for design" Statistics in Medicine 11:1685-1704
#   Vickers (2001) "The use of percentage change from baseline as an outcome in a controlled trial is
#     statistically inefficient: a simulation study" BMC Med Res Methodol. 2001; 1 (1): 6
# requirement: need library 'mvtnorm'

### power calculation
power.genotype.conti <- function(N, Rep=2000, alpha=.05, ...){
  pval <- sapply(rep(N, Rep), FUN=simu.genotype.conti, ...)
  retval <- list()
  retval$power <- mean(pval<=alpha)
  retval$ci    <- ci.binom(pval<=alpha)
  retval$call  <- match.call()
  class(retval) <- "power.genoytpe.conti"
  retval
}

print.power.genoytpe.conti <- function(x, ...)
  {
    cat("\n")
    cat("Power calculation for a study including genotyping\n")
    cat("\n")
    cat("call:\n")
    print(x$call)
    cat("\n")
    cat("Estimated power:", x$power, "\n")
    cat("\n")
    print(x$power)
    print("Simulation confidence region for estimated power:\n")
    cat("\n")
    print(x$ci)
  }


### simulation
simu.genotype.conti <- function (N, p=0.15, pi=0, me1=50, me2=me1, delta=-5, sd1=10, sd2=10, verbose=FALSE,
                                 minh=c('additive', 'dominant','recessive'), genotype.delta=TRUE, Factor=FALSE) 
{
  minh <- match.arg(minh)
  if ( min(N, sd1, sd2)<0 ) stop('N, sd1, and sd2 must be greater than 0')
  if ( p<=0 | abs(pi)<0 | p>=1 | abs(pi)>1 ) stop('p and abs(pi) must be between 0 and 1.') 
  f.mod <- switch(minh,
         dominant       = c(      0,     1, 1),
         additive       = c(      0 ,  0.5, 1),
         recessive      = c(      0,     0, 1)  ) 
  q <- 1 - p
  fhw <- c(q^2, 2*p*q, p^2)  # major allele first
  nhw <- round(N*fhw)
  if (sum(nhw)!=N) nhw[3] <- N-sum(nhw[1:2])
  covm <- matrix(c(1,    pi,    pi,    1   ), nr=2)*
          matrix(c(sd1^2, sd1*sd2, sd1*sd2, sd2^2), nr=2)
  if (!genotype.delta) delta <- delta/sum(fhw*f.mod)  # convert to overall delta due to all genotypes

  ## explicitly create x1:x3, t1:t3 to avoid R CMD check warning for R > 2.6.0
  x1 <- x2 <- x3 <- NULL
  t1 <- t2 <- t3 <- NULL
  
  for (i in 1:3) {
    if (nhw[i]!=0) {
      assign(paste('x',i,sep=''), rmvnorm(n=nhw[i], mean=c(me1,me2+f.mod[i]*delta), sigma=covm))
      assign(paste('t',i,sep=''), rep(f.mod[i],nhw[i]))
    } else {
      assign(paste('x',i,sep=''), NULL)
      assign(paste('t',i,sep=''), NULL)
    }
  }
  x <- data.frame(rbind(x1,x2,x3), Trt=c(t1,t2,t3))
  if (Factor) x$Trt <- as.factor(x$Trt)
  colnames(x) <- c('X', 'Y', 'Trt')
  mod <- lm(Y~X+Trt, data=x)   # ANCOVA
  if (verbose) {
    print( summary(mod) )
    plot(Y~X+Trt, data=x)
  }
  ret <- anova(mod, test='F')['Trt','Pr(>F)']
  ret
}




