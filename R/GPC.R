# Genetics power calculator for linear trend association studies
#
# Purcell's Power estimation method
# (http://pngu.mgh.harvard.edu/~purcell/gpc/cc2.html)

# Linear Tread Test
# http://linkage.rockefeller.edu/pawe3d/help/Linear-trend-test-ncp.html


# Given
# pA -- High risk allele frequency (A)
# pD -- disease prevalence
# RRAa -- Genotype relative risk Aa = RR(Aa|aa)=Pr(D|Aa)/Pr(D|aa)
# RRAA -- Genotype relative risk AA = RR(AA|aa)=Pr(D|AA)/Pr(D|aa)
# Dprime -- LD measure
# pB -- Marker allele frequency (B)
# nCase -- Number of cases
# ratio -- Control: case ratio = nControl/nCase
# alpha -- User-defined type I error rate
GPC.default<-function(pA, pD, RRAa, RRAA, Dprime, pB, 
                   nCase=500, ratio=1, alpha=0.05, quiet=FALSE)
{
  if(!(pA>0 && pA<1) || !(pB>0 && pB<1) || !(Dprime>=0 && Dprime<=1)
     || !(pD>0 && pD<1) || !(RRAa>1) || !(RRAA>1) ||
     !(nCase>1) || !(ratio>0) || !(alpha>0 && alpha<0.5))
  { cat("Some arguments are out of bounds!\n") 
    cat("0<pA<1\n") 
    cat("0<pD<1\n") 
    cat("1<RRAa\n") 
    cat("1<RRAA\n") 
    cat("0<=Dprime<=1\n") 
    cat("0<pB<1\n") 
    cat("1<nCase\n") 
    cat("0<ratio\n") 
    cat("0.05<alpha<0.5\n") 
    stop("Program ends due to input errors!\n")
  }
  # get penetrances Pr(D|aa), Pr(D|Aa), Pr(D|AA)
  pa<-1-pA
  pb<-1-pB
  denom<-RRAA*pA^2+RRAa*2*pA*pa+pa^2
  PrDgaa<-pD/denom

  PrDgAa<-RRAa*PrDgaa
  PrDgAA<-RRAA*PrDgaa

  pen<-c(PrDgaa, PrDgAa, PrDgAA)

  # estimate haplotype frequencies Pr(AB), Pr(Ab), Pr(aB), and Pr(ab)
  # based on pA, pB, and Dprime. We assume D>0
  myHapFreqs<-hapFreq(Dprime, pA, pB)

  # estimate the sampling probabilities Pr(BB|D), Pr(Bb|D), Pr(bb|D)
  tmp<-samplingProb(myHapFreqs, pen, pD, pA, pB)
  PrBBgD<-tmp[1]
  PrBbgD<-tmp[2]
  PrbbgD<-tmp[3]
  PrBBgDbar<-tmp[4]
  PrBbgDbar<-tmp[5]
  PrbbgDbar<-tmp[6]

  PrDgBB<-tmp[7]
  PrDgBb<-tmp[8]
  PrDgbb<-tmp[9]

  PrBgD<-tmp[10]
  PrbgD<-tmp[11]
  PrBgDbar<-tmp[12]
  PrbgDbar<-tmp[13]

  n0c<-(nCase*PrbbgD)
  n1c<-(nCase*PrBbgD)
  n2c<-nCase-n0c-n1c

  nControl<-(nCase*ratio)

  n0n<-(nControl*PrbbgDbar)
  n1n<-(nControl*PrBbgDbar)
  n2n<-nControl-n0n-n1n

  nc.vec<-c(n2c, n1c, n0c)
  nn.vec<-c(n2n, n1n, n0n)

  mat<-data.frame(case=nc.vec, control=nn.vec, code=c(2,1,0))

  x.vec<-c(2, 1, 0)

  R<-sum(nc.vec)
  S<-sum(nn.vec)

  n0<-n0c+n0n
  n1<-n1c+n1n
  n2<-n2c+n2n
  n.vec<-c(n2, n1, n0)

  N<-R+S

  p1.vec<-c(PrBBgD, PrBbgD, PrbbgD)
  p0.vec<-c(PrBBgDbar, PrBbgDbar, PrbbgDbar)

  numer<-sum(x.vec*(p1.vec-p0.vec))
  numer<-numer^2
  part1<-sum(x.vec^2*(R*p0.vec+S*p1.vec))
  part2<-sum(x.vec*(R*p0.vec+S*p1.vec))
  denom<-part1-part2^2/N

  # non-centrality parameter for the linear trend test
  ncp<-R*S*numer/denom

  # format outputs
  # Case-control parameters
  mat.para<-matrix(0, nrow=7, ncol=1)
  rownames(mat.para)<-c("Number of cases", "Number of controls", 
    "High risk allele frequency (A)", "Prevalence",
    "Genotypic relative risk Aa", "Genotypic relative risk AA",
    "Genotypic risk for aa (baseline)")
  mat.para[1,1]<-nCase
  mat.para[2,1]<-nControl
  mat.para[3,1]<-pA
  mat.para[4,1]<-pD
  mat.para[5,1]<-RRAa
  mat.para[6,1]<-RRAA
  mat.para[7,1]<-PrDgaa

  # Marker locus B
  mat.B<-matrix(0,nrow=7, ncol=1)
  rownames(mat.B)<-c("High risk allele frequency (B)",
    "Linkage disequilibrium (D')", "Penetrance at marker genotype bb",
    "Penetrance at marker genotype Bb", "Penetrance at marker genotype BB",
    "Genotypic odds ratio Bb", "Genotypic odds ratio BB") 
  mat.B[1,1]<-pB
  mat.B[2,1]<-Dprime
  mat.B[3,1]<-PrDgbb
  mat.B[4,1]<-PrDgBb
  mat.B[5,1]<-PrDgBB
  # OR(Bb|bb)
  mat.B[6,1]<-PrBbgD*PrbbgDbar/(PrBbgDbar*PrbbgD)
  # OR(BB|bb)
  mat.B[7,1]<-PrBBgD*PrbbgDbar/(PrBBgDbar*PrbbgD)

  # Expected allele frequencies Pr(B|D), Pr(b|D), Pr(B|\bar{D}), Pr(b|\bar{D})
  mat.aFreq<-matrix(0,nrow=2,ncol=2)
  rownames(mat.aFreq)<-c("B","b")
  colnames(mat.aFreq)<-c("Case","Control")
  mat.aFreq[1,1]<-PrBgD
  mat.aFreq[1,2]<-PrBgDbar
  mat.aFreq[2,1]<-PrbgD
  mat.aFreq[2,2]<-PrbgDbar

  # Expected genotype frequencies
  # Pr(BB|D), Pr(Bb|D), Pr(bb|D)
  # Pr(BB|\bar{D}), Pr(Bb|\bar{D}), Pr(bb|\bar{D})
  mat.gFreq<-matrix(0, nrow=3, ncol=2)
  rownames(mat.gFreq)<-c("BB","Bb", "bb")
  colnames(mat.gFreq)<-c("Case", "Control")
  mat.gFreq[1,1]<-PrBBgD
  mat.gFreq[1,2]<-PrBBgDbar
  mat.gFreq[2,1]<-PrBbgD
  mat.gFreq[2,2]<-PrBbgDbar
  mat.gFreq[3,1]<-PrbbgD
  mat.gFreq[3,2]<-PrbbgDbar

  alpha.vec<-c(0.1, 0.05, 0.01, 0.001, alpha)
  power.vec<-rep(0,5)
  for(i in 1:5)
  { a<-alpha.vec[i]
    const<-qchisq(1-a, df=1)
    power.vec[i]<-1-pchisq(const, df=1, ncp=ncp)
  }

  # Case-Control statistics
  mat.stat<-cbind(alpha.vec, power.vec)
  colnames(mat.stat)<-c("Alpha", "Power")
  rownames(mat.stat)<-rep("",5)
 
  res<-list(power=power.vec[5], ncp=ncp,
            mat.para=mat.para, mat.B=mat.B, mat.aFreq=mat.aFreq,
            mat.gFreq=mat.gFreq, mat.stat=mat.stat)

  if(quiet==FALSE)
  {
    cat("\n Case-control parameters>>\n");
    print(mat.para)
    cat("\n Marker locus B>>\n");
    print(mat.B)
    cat("\n Expected allele frequencies>>\n");
    print(mat.aFreq)
    cat("\n Expected genotype frequencies>>\n");
    print(mat.gFreq)
    cat("\n Case-control statistics>>\n");
    print(mat.stat)
    cat("\n power (alpha=",alpha, ")=", power.vec[5], " ncp=", ncp, "\n")
  }

  invisible(res)
}

GPC<-function(pA, pD, RRAa, RRAA, r2, pB, 
                   nCase=500, ratio=1, alpha=0.05, quiet=FALSE)
{
  if(!(pA>0 && pA<1) || !(pB>0 && pB<1) || !(r2>=0 && r2<=1)
     || !(pD>0 && pD<1) || !(RRAa>1) || !(RRAA>1) ||
     !(nCase>1) || !(ratio>0) || !(alpha>0 && alpha<0.5))
  { cat("Some arguments are out of bounds!\n") 
    cat("0<pA<1\n") 
    cat("0<pD<1\n") 
    cat("1<RRAa\n") 
    cat("1<RRAA\n") 
    cat("0<=r2<=1\n") 
    cat("0<pB<1\n") 
    cat("1<nCase\n") 
    cat("0<ratio\n") 
    cat("0.05<alpha<0.5\n") 
    stop("Program ends due to input errors!\n")
  }

  # estimate Dprime based on r2, pA, pB
  Dprime<-Dprime.fun2(r2, pA, pB)

  res<-GPC.default(pA, pD, RRAa, RRAA, Dprime, pB, 
                   nCase, ratio, alpha, quiet)

  invisible(res)
}

# estimate haplotype frequencies Pr(AB), Pr(Ab), Pr(aB), and Pr(ab) based
# on pA, pB, and Dprime. We assume D>0
hapFreq<-function(Dprime, pA, pB)
{
  pa<-1-pA
  pb<-1-pB
  dmax<-min(c(pA*pb, pa*pB))
  D<-Dprime * dmax

  PrAB<-pA*pB+D
  PraB<-pa*pB-D
  PrAb<-pA*pb-D
  Prab<-pa*pb+D

  return(c(PrAB, PraB, PrAb, Prab))
}


# estimate the sampling probabilities Pr(BB|D), Pr(Bb|D), Pr(bb|D)
samplingProb<-function(myHapFreqs, pen, pD, pA, pB)
{
  pa<-1-pA
  pb<-1-pB
  # penetrances Pr(D|aa), Pr(D|Aa), Pr(D|AA)
  PrDgaa<-pen[1]
  PrDgAa<-pen[2]
  PrDgAA<-pen[3]

  # haplotype frequencies Pr(AB), Pr(aB), Pr(Ab), Pr(ab)
  PrAB<-myHapFreqs[1]
  PraB<-myHapFreqs[2]
  PrAb<-myHapFreqs[3]
  Prab<-myHapFreqs[4]
  
  # sampling probabilities for cases
  # Pr(BB|D), Pr(Bb|D), Pr(bb|D)
  numer<-PrDgAA*PrAB^2+PrDgAa*2*PrAB*PraB+PrDgaa*PraB^2
  PrBBgD<-numer/pD

  numer<-PrDgAA*2*PrAB*PrAb+PrDgAa*2*(PrAB*Prab+PrAb*PraB)+PrDgaa*2*PraB*Prab
  PrBbgD<-numer/pD

  numer<-PrDgAA*PrAb^2+PrDgAa*2*PrAb*Prab+PrDgaa*Prab^2
  PrbbgD<-numer/pD

  PrDgBB<-PrBBgD*pD/(pB^2)
  PrDgBb<-PrBbgD*pD/(2*pB*pb)
  PrDgbb<-PrbbgD*pD/(pb^2)

  # sampling probabilities for controls
  # Pr(BB|\bar{D}), Pr(Bb|\bar{D}), Pr(bb|\bar{D})
  PrBBgDbar<-(1-PrDgBB)*pB^2/(1-pD)
  PrBbgDbar<-(1-PrDgBb)*2*pB*pb/(1-pD)
  PrbbgDbar<-(1-PrDgbb)*pb^2/(1-pD)

  # Expected allele frequencies Pr(B|D), Pr(b|D), Pr(B|Dbar), Pr(b|Dbar)
  PrBgD<-PrBBgD+PrBbgD/2
  PrbgD<-PrbbgD+PrBbgD/2
  PrBgDbar<-PrBBgDbar+PrBbgDbar/2
  PrbgDbar<-PrbbgDbar+PrBbgDbar/2

  return(c(PrBBgD, PrBbgD, PrbbgD, PrBBgDbar, PrBbgDbar, PrbbgDbar,
           PrDgBB, PrDgBb, PrDgbb, PrBgD, PrbgD, PrBgDbar, PrbgDbar))
}


# r2 -- LD measure r^2 between SNP 1 and SNP 2
# pA -- frequency of minor allele (A) for SNP 1
# pB -- frequency of minor allele (B) for SNP 2
#
# D = pA.pB - pAB
# 
# dmax = min(pA.(1-pB),(1-pA).pB)
# 
# dmin = max(-pA.pB, -(1-pA).(1-pB))
# 
# if D < 0, D' = D/dmin else D' = D/dmax
# 
# r2 = D.D/(pA.pB.(1-pA).(1-pB) 
#
# D' = |r|*sqrt(pA*pB*(1-pA)*(1-pB)) / dmax if D > 0
# D' = -|r|*sqrt(pA*pB*(1-pA)*(1-pB)) / dmin if D < 0
# 
#
# suppose that D < 0
Dprime.fun1<-function(r2, pA, pB)
{ tmpr2<-r2.upp1(pA, pB)
  if(r2>tmpr2)
  { msg<-paste("r2 = ", r2, " > upper bound of r2 = ", tmpr2, ". r2 is changed to floor(tmpr2*100)/100!\n");
    warning(msg);
    r2<-floor(tmpr2*100)/100
  }
  numer<- - sqrt(r2*pA*(1-pA)*pB*(1-pB));
  dmin<- max(-pA*pB, -(1-pA)*(1-pB));
  res<-numer/dmin;
  return(res)
}

# suppose that D > 0
Dprime.fun2<-function(r2, pA, pB)
{ 
  tmpr2<-r2.upp2(pA, pB)
  if(r2>tmpr2)
  { msg<-paste("r2 = ", r2, " > upper bound of r2 = ", tmpr2, ". r2 is changed to floor(tmpr2*100)/100!\n");
    warning(msg);
    r2<-floor(tmpr2*100)/100
  }
  numer<- sqrt(r2*pA*(1-pA)*pB*(1-pB));
  dmax<- min(pA*(1-pB), pB*(1-pA));
  res<-numer/dmax;
  return(res)
}

# upper bound of r2 given pA, pB for D<0
r2.upp1<-function(pA, pB)
{
  numer<-pA*pB*(1-pA)*(1-pB)
  dmin<- max(-pA*pB, -(1-pA)*(1-pB));
  res<-numer/(dmin^2)
  return(1/res)
}


# upper bound of r2 given pA, pB for D>0
r2.upp2<-function(pA, pB)
{
  numer<-pA*pB*(1-pA)*(1-pB)
  dmax<- min(pA*(1-pB), pB*(1-pA));
  res<-numer/(dmax^2)
  res<-dmax^2/numer
  return(res)
}

