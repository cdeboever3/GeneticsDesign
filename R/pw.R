## contributed by Michael Man.
## TODO: Needs documentation before exporting.


### simple simulation for two group design
pw <- function(n1, n2=n1*(1-fc)/fc, fc=.5, pi=0, me1=50, me2=45, sd1=10, sd2=10, TEST=F){
  covm <- matrix(c(1,    pi,    pi,    1   ), nr=2)*
          matrix(c(sd1^2, sd1*sd2, sd1*sd2, sd2^2), nr=2)
  x1 <- rmvnorm(n=n1, mean=c(me1,me1), sigma=covm)
  x2 <- rmvnorm(n=n2, mean=c(me1,me2), sigma=covm)
  x <- data.frame(rbind(x1,x2), Trt=c(rep(0,n1), rep(1,n2)))
  colnames(x) <- c('X', 'Y', 'Trt')
  mod <- lm(Y~X+Trt, data=x)
  if (TEST) {
    print( summary(mod) )
    plot(Y~X+Trt, data=x)
  }
  ret <- anova(mod, test='F')['Trt','Pr(>F)']
  ret
}

