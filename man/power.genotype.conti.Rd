\name{power.genotype.conti}
\alias{power.genotype.conti}
\alias{simu.genotype.conti}
\title{power for genetic studies using baseline measure}
\description{
  Estimate power for genetice studies using baseline measurements via
  simulation. 
}
\usage{
power.genotype.conti(N, Rep = 2000, alpha = 0.05, ...)
simu.genotype.conti(N, p=0.15, pi=0, me1=50, me2=me1, delta=-5,
                    sd1=10, sd2=10, verbose=FALSE,
                    minh=c('additive', 'dominant', 'recessive'),
                    genotype.delta=TRUE, Factor=FALSE) 
}
\arguments{
  \item{N}{total number of subjects}
  \item{p}{frequency of A (affected) allele}
  \item{Rep}{number of simulatin runs used to estimate power}
  \item{alpha}{significance level}
  \item{pi}{ correlation coefficient}
  \item{me1, me2}{mean of control and treatment groups}
  \item{delta}{treatment/genotype effect}
  \item{sd1,sd2}{standard deviation of the control and treatment groups}
  \item{minh}{mode of inheritance, one of 'additive', 'dominant', or 'recessive'}
  \item{genotype.delta}{logical indicating whether the treatment effect
    occurs only for an individual genotype (\code{genotype.delta=TRUE})
    or for all genotypes (\code{genotype.delta=FALSE}) }
  \item{Factor}{Should the simulated treatment variable 'Trt' be 
    be treated as a factor variable (\code{Factor=TRUE}) or as a
    numeric variable (\code{Factor=FALSE}).}
  \item{verbose}{Should information about each simulated data set and
    model fit be displayed.}
  \item{\dots}{Arguments to be passed to \code{simu.genotype.conti}}
}
\details{

}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{
   Frison and Pocock (1992) "Repeated measures in clinical trials: analysis using mean summary statistics
   and its implications for design" Statistics in Medicine 11:1685-1704
   
   Vickers (2001) "The use of percentage change from baseline as an outcome in a controlled trial is
   statistically inefficient: a simulation study" BMC Med Res
   Methodol. 2001; 1 (1): 6
  }
\author{Michael Man, minor changes by Gregory R. Warnes
  \email{greg@random-technologies-llc.com} }
\seealso{ \code{\link{power.casectrl}} }
\examples{

\dontrun{
  # use defaults, 100 subjects
  power.genotype.conti(N=100)

  # same calculation, specifying all values
  power.genotype.conti(N=100, Rep=2000, p=0.15, pi=0, me1=50, me2=50, delta=-5,
                       sd1=10, sd2=10, verbose=FALSE, minh='additive',
                       genotype.delta=TRUE, Factor=FALSE) 

  # Show details for small simulation study
  power.genotype.conti(N=10, verbose=TRUE)
}
\dontshow{
  # test code
  set.seed(100)
  power.genotype.conti(N=100, Rep=10)

  power.genotype.conti(N=100, Rep=10, p=0.15, pi=0, me1=50, me2=50, delta=-5,
  sd1=10, sd2=10, verbose=FALSE, minh='additive',
  genotype.delta=TRUE, Factor=FALSE) 
}


}
\keyword{design}

