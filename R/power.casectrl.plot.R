power.casectrl.plot <- function (N, gamma=1.6, p=1:9/10, kp=0.1, alpha=0.05, fc=0.5,
                                 minh=c('multiplicative', 'dominant','recessive'),
                                 Nsnp=1, vary=c('prevalence','SNPs'), ylim=c(0,1), PLOT=T, ... )
{
  minh <- match.arg(minh)
  vary <- match.arg(vary)
  if (length(p)<2) stop('Must have more than 1 value in p.')
  if (length(kp) > 1 & length(Nsnp) > 1) stop("Nsnps and kp can't be all > 1.")
  if (vary=='prevalence') {
    cmd <- expression(tapply(p, p, function(x, ...) power.casectrl(p=x,...), N=N, gamma=gamma, kp=kp[j], alpha=alpha/Nsnp, fc=fc, minh=minh))
    Xvary <- kp 
  } else if (vary=='SNPs') {
    cmd <- expression(tapply(p, p, function(x, ...) power.casectrl(p=x,...), N=N, gamma=gamma, kp=kp, alpha=alpha/Nsnp[j], fc=fc, minh=minh))
    Xvary <- Nsnp 
  }
  J <- length(Xvary)
  ret <- matrix(NA, nc=J, nr=length(p))
  colnames(ret) <- paste(vary, '=', Xvary)
  for (j in 1:J) ret[,j] <- eval(cmd)

  if (PLOT) {
    nc <- 1:ncol(ret)
    subt <- paste("( RR", gamma, "; total subjects", N,"; SNPs", Nsnp[1], "; prevalence", kp[1],
                   "; mode of inheritance:", minh, "; overall sig.level", alpha, ")" )
    matplot(p, ret, type="l", ylim=ylim, lty=1, col=nc, xlab="Allele Frequency", ylab="Power", sub=subt, ...)
    abline(h=c(.8), lty=1)
    legend( locator(1), colnames(ret), lty=1, col=nc )
  }
  ret
}

