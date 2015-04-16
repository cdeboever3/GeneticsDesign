GeneticPower.Quantitative.Numeric <- function(N=1000,
                                              delta=1,
                                              freq=0.15,
                                              minh=c('additive','dominant','recessive'),
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
## minh = mode of inheritance:  "recessive", "additive", "dominant" same as moi=0,0.5,and 1.0, respectively.
##        defaults to "additive" if no moi specified either
## sigma = standard deviation of the response phenotype
## OtherParms = the number of additional parameters (really, DOF) in the model that will reduce your overall DOF
## alpha = the desired significance level
## numtests = the number of tests to be corrected by Bonferroni adjustment beforee achieving 'alpha'
## moi = mode of inheritance: 0 for recessive, 0.5 for additive, 1.0 for dominant, or anywhere in between,
##       this OVER-RIDES minh...useful for modeling in-between moi's...
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

    if (is.null(rsquared)) {
      
      xbar <- (N2 + N3*2)/N;
      ybar <- (N2*mu[2] + N3*mu[3])/N;
      
      slope <- (N1*xbar*ybar + (1-xbar)*N2*(mu[2]-ybar) + (2-xbar)*N3*(mu[3]-ybar))/(N1*xbar^2 + N2*(1-xbar)^2 + N3*(2-xbar)^2);
      a <- ybar - slope*xbar;
      
      fits <- a + slope*c(0,1,2);
      
      eNs <- c(N1,N2,N3);
      
      Rsqtop <- sum(eNs*(fits-ybar)^2);
      
      Rsqbottom <- Rsqtop + sum((eNs-1)*sigma^2 + eNs*(mu - fits) + eNs*(mu-fits)^2);
      
      NewRsquared <- Rsqtop/Rsqbottom;
      
    } else {
      
      NewRsquared <- rsquared;
    } 
    
    power <- pf(qf(1-alphy,df1=1,df2=(N-2-OtherParms)),ncp=(NewRsquared*(N-2)/((1-NewRsquared))),df1=1,df2=(N-2-OtherParms),lower.tail=F);

    return(power);
}




