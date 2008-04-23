\name{GPC}
\alias{GPC}
\alias{GPC.default}
\title{Genetics power calculator for linear trend association studies}
\description{
  Genetics power calculator for linear trend association studies.
}
\usage{
  GPC(pA, pD, RRAa, RRAA, r2, pB, 
                   nCase=500, ratio=1, alpha=0.05, quiet=FALSE)
  GPC.default(pA, pD, RRAa, RRAA, Dprime, pB, 
                   nCase=500, ratio=1, alpha=0.05, quiet=FALSE)
}
\arguments{
  \item{pA}{High risk allele frequency (\code{A}).}

  \item{pD}{Disease prevalence.}

  \item{RRAa}{Genotype relative risk (\code{Aa}) = \code{RR(Aa|aa)=Pr(D|Aa)/Pr(D|aa)}.}

  \item{RRAA}{Genotype relative risk (\code{AA}) = \code{RR(AA|aa)=Pr(D|AA)/Pr(D|aa)}.}

  \item{r2}{LD measure. Assume that \code{D>0}.}

  \item{Dprime}{LD measure.}

  \item{pB}{Marker allele frequency (\code{B}).}

  \item{nCase}{Number of cases.}

  \item{ratio}{Control:case ratio \code{= nControl/nCase}.}

  \item{alpha}{User-defined type I error rate.}

  \item{quiet}{Print some intermediate results if \code{quiet=FALSE}.}
}
\details{
  The power is for the test that disease is associated with a marker, given high risk allele frequency (\code{A}), disease prevalence, genotype relative risk (\code{Aa}), genotype relative risk (\code{AA}), LD measure (\code{D'} or \code{r^2}), marker allele frequency (\code{B}), number of cases, control:case ratio, and probability of the Type I error. The linear trend test (Cochran 1954; Armitage 1955) is used.
}
\value{
  \item{power}{The estimated power for the association test.}

  \item{ncp}{Non-centrality parameter.}

  \item{mat.para}{A matrix of case-control parameters, including number of cases, number of controls, high risk allele frequency, prevalence, genotypic relative risk (\code{Aa}), genotypic relative risk (\code{AA}), genotypic risk for \code{aa} (baseline).}

  \item{mat.B}{A matrix of marker locus \code{B} parameters, including marker allele frequency, linkage disequilibrium (\code{D'}), penetrance at marker genotype \code{bb}, penetrance at marker genotype \code{Bb}, penetrance at marker genotype \code{BB}, genotypic odds ratio \code{Bb}, genotypic odds ratio \code{BB}.}

  \item{mat.aFreq}{A 2 by 2 matrix of expected allele frequencies \code{Pr(B|D), Pr(b|D), Pr(B|non D), Pr(b|non D)}.}

  \item{mat.gFreq}{A 3 by 2 matrix of expected genotype frequencies \code{Pr(BB|D), Pr(Bb|D), Pr(bb|D), Pr(BB|non D), Pr(Bb|non D), Pr(bb|non D)}.}

  \item{mat.stat}{Power estimates for a sequence of Type I errors.}
}
\references{
  Armitage, P. (1955) 
  Tests for linear trends in proportions and frequencies. 
  \emph{Biometrics}, 11, 375-386.

  Cochran, W.G. (1954) 
  Some methods for strengthening the common chi-squared tests. 
  \emph{Biometrics}, 10, 417-451.

  Gordon D, Finch SJ, Nothnagel M, Ott J (2002) 
  Power and sample size calculations for case-control genetic association tests when errors are present: application to single nucleotide polymorphisms. 
  \emph{Hum. Hered.}, 54:22-33.

  Gordon D, Haynes C, Blumenfeld J, Finch SJ (2005) 
  PAWE-3D: visualizing Power for Association With Error in case/control genetic studies of complex traits. 
  \emph{Bioinformatics}, 21:3935-3937.

  Purcell S, Cherny SS, Sham PC. (2003).
  Genetic Power Calculator: design of linkage and association genetic mapping 
  studies of complex traits. 
  \emph{Bioinformatics}, 19(1):149-150.

  Sham P. (1998).
  \emph{Statistics in Human Genetics.}
  Arnold Applications of Statistics.
}
\author{
  Weiliang Qiu \email{stwxq@channing.harvard.edu},
  Ross Lazarus \email{ross.lazarus@channing.harvard.edu}
}
\examples{
  res1<-GPC(pA=0.05, pD=0.1, RRAa=1.414, RRAA=2, r2=0.9, pB=0.06, 
                   nCase=500, ratio=1, alpha=0.05, quiet=FALSE)

  res2<-GPC.default(pA=0.05, pD=0.1, RRAa=1.414, RRAA=2, Dprime=0.9, pB=0.06, 
                   nCase=500, ratio=1, alpha=0.05, quiet=FALSE)
}
\keyword{htest}

