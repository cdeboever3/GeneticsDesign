"power.casectrl" <-
function (N, gamma = 4.5, p = 0.15, kp = 0.1, alpha = 0.05, fc = 0.5, 
    minh = c("multiplicative", "dominant", "recessive")) 
{
    minh <- match.arg(minh)
    if (!all(gamma > 0, N > 0)) 
        stop("N and gamma must be greater than 0")
    if (min(p, kp, alpha, fc) <= 0 | max(p, kp, alpha, fc) >= 
        1) 
        stop("p, kp, alpha, and fc must be between 0 and 1.")
    f.mod <- switch(minh, multiplicative = c(gamma^2, gamma, 
        1), dominant = c(0, 0, 1), recessive = c(0, 1, 1))
    q <- 1 - p
    fhw <- c(p^2, 2 * p * q, q^2)
    pi <- kp/sum(f.mod * fhw)
    if (pi <= 0 | pi >= 1) {
        warning("The combination of p, kp, and gamma produces an unrealistic value of pi.")
        ret <- NA
    }
    else {
        fe <- rbind(fhw, fhw)
        dimnames(fe) <- list(c("Case", "Control"), c("AA", "Aa", 
            "aa"))
        f <- fe * rbind(f.mod * pi, 1 - f.mod * pi)
        Pct <- apply(f, 1, sum)
        f2 <- f * c(fc, 1 - fc)/Pct
        fe2 <- fe * c(fc, 1 - fc)
        fe2
        apply(fe2, 1, sum)
        f2
        apply(f2, 1, sum)
        lambda <- sum((f2 - fe2)^2/fe2) * N
        ret <- 1 - pchisq(qchisq(1 - alpha, df = 1), df = 1, 
            ncp = lambda, lower.tail = TRUE)
    }
    ret
}
