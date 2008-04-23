GeneticPower.Quantitative.Factor <- function(N=1000,
                                             delta=1,
                                             freq=0.15,
                                             minh=c("additive","dominant","recessive"),
                                             sigma=1,
                                             OtherParms=0,
                                             alpha=0.05,
                                             numtests=1,
                                             moi=NULL,
                                             rsquared=NULL) 
{
## N = total samples in the analysis
## delta = mean(bb) - mean (AA), where 'b' is the disease allele, 'A' is the reference allele
## freq = allele frequency of disease allele 'b'
## minh = mode of inheritance:  "recessive", "additive", "dominant" same as moi=0,0.5, and 1.0, respectively
##        defaults to "additive" if no moi specified eith
## sigma = standard deviation of the response phenotype
## OtherParms = the number of additional parameters (really, DOF) in the model that will reduce your overall DOF
## alpha = the desired significance level
## numtests = the number of tests to be corrected by Bonferroni adjustment beforee achieving 'alpha'
## moi = mode of inheritance: 0 for recessive, 0.5 for additive, 1.0 for dominant, or anywhere in between,
##       this OVER-RIDES minh...useful for modeling i-between moi's...
## rsquared = fraction of total sum-of-squares explained by fit. OVER-RIDES delta AND sigma.

    alphy <- 1-(1-alpha)^(1/numtests);

    if (is.null(moi)) {
      minh <- match.arg(minh) # can use abbreviated names
      moi <- switch(minh,
                    additive   = 0.5,
                    dominant   = 1.0,
                    recessive  = 0  ) # define the mode of inheritance
    }
  
    N1 <- N*(1-freq)^2;
    N2 <- 2*N*freq*(1-freq);
    N3 <- N*freq^2;
     
    mu <- c(0,moi,1)*delta; 
    mu.bar <- (mu[1]*N1+mu[2]*N2+mu[3]*N3)/N;

    lambda <- ifelse(is.null(rsquared),
                     (N1*(mu[1]-mu.bar)^2+N2*(mu[2]-mu.bar)^2+N3*(mu[3]-mu.bar)^2)/(sigma^2),
                     (N-3)*rsquared/(1-rsquared)
                     );

    power <- pf(qf(1-alphy,df1=2,df2=(N-3-OtherParms)),ncp=lambda,df1=2,df2=(N-3-OtherParms),lower.tail=F)

    return(power)
}
