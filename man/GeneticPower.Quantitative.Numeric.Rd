\name{GeneticPower.Quantitative.Numeric}
\alias{GeneticPower.Quantitative.Numeric}
\alias{GeneticPower.Quantitative.Factor}
\title{Power of Genetics Study}
\description{
  Compute power of quantitative genetics studies, when the genotype is
  handled as a numeric value (0,1,2)
  \code{GeneticPower.Quantitative.Numeric} or as a factor
  \code{GeneticPower.Quantitative.Factor}.
  
}
\usage{
GeneticPower.Quantitative.Numeric(
                    N=1000,
                    delta=1,
                    freq=0.15,
                    minh=c("additive", "dominant", "recessive"),
                    sigma=1,
                    OtherParms=0,
                    alpha=0.05,
                    numtests=1,
                    moi=NULL,
                    rsquared=NULL)
GeneticPower.Quantitative.Factor( 
                    N=1000,
                    delta=1,
                    freq=0.15,
                    minh=c("additive","dominant","recessive"),
                    sigma=1,
                    OtherParms=0,
                    alpha=0.05,
                    numtests=1,
                    moi=NULL, 
                    rsquared=NULL) 
}
\arguments{
  \item{N}{total samples in the analysis}
  \item{delta}{Treatment effect for an individual homozygote for the
    disease allele ('b') relative to an individual homozygote for the
    reference allele ('A')}
  \item{freq}{allele frequency of disease allele 'b'}
  \item{minh}{mode of inheritance:  "additive","dominant","recessive",
    Default is "additive". \emph{This parameter is OVER-RIDDEN by \code{moi}.}}
  \item{sigma}{standard deviation of the response phenotype}
  \item{OtherParms}{number of additional parameters (really, DOF) in the
    model that will reduce your overall DOF} 
  \item{alpha}{desired significance level}
  \item{numtests}{number of tests to be corrected by Bonferroni adjustment beforee achieving 'alpha'}
  \item{moi}{continuous value between 0 and 1 (inclusive) specifying the
    mode of inheritance: 0 for recessive, 0.5 for additive, 1.0
    for dominant. \emph{This parameter OVER-RIDES \code{minh}.}
    }
  \item{rsquared}{fraction of total sum-of-squares explained by
    fit. \emph{This parameter OVER-RIDES \code{delta} AND
      \code{sigma}.}}
  }
\details{
  The value of \code{moi} overrides any value specified for
  \code{minh}. Specifying a \code{minh="recessive"} is equivalent to
  specifying \code{moi=0}, \code{minh="additive"} is equivalent to
  \code{moi=0.5}, and \code{minh="dominant"} is equivalent to
  \code{moi=1.0}.
}
% \value{
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
% }
%\references{ ~put references to the literature/web site here ~ }
\author{
  Craig L.Hyde \email{Craig.L.Hyde@pfizer.com}
  and Feng Gao \email{feng.gao1@pfizer.com}
}
%\note{ ~~further notes~~  }

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

GeneticPower.Quantitative.Numeric(
                                  N=50,
                                  freq=0.1,
                                  minh="recessive",
                                  alpha=0.05
                                  )

GeneticPower.Quantitative.Factor(
                                  N=50,
                                  freq=0.1,
                                  minh="recessive",
                                  alpha=0.05
                                 )


## 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}
%\keyword{genetics}
