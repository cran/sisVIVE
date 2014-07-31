\name{plot.sisVIVE}
\alias{plot.sisVIVE}
\title{Plots the entire solution path for coefficients estimated for the sisVIVE object}
\description{ Produces a plot of estimates for alpha and beta across lambda(s)}
\usage{\method{plot}{sisVIVE}(x, ...)}
\arguments{
  \item{x}{a sisVIVE object}
  \item{...}{Additional arguments for generic plot}
}
\details{The default plot shows the estimated coefficients for beta. Beta, the causal effect of exposure on outcome, is denoted by a solid black line.}
\value{NULL}
\author{Hyunseung Kang}
\seealso{predict, summary, and cv.sisVIVE}
\examples{
library(MASS)
library(lars)

n = 1000; L = 10; s= 3;
Z <- matrix(rnorm(n*L),n,L)
error <- mvrnorm(n,rep(0,2),matrix(c(1,0.8,0.8,1),2,2))
intD = rnorm(1); ZtoD =   rnorm(L,0,1); ZtoY = c(rnorm(s),rep(0,L-s)); DtoY = 1; intY = rnorm(1)
D = intD + Z \%*\% ZtoD + error[,1]
Y = intY + Z \%*\% ZtoY + D * DtoY + error[,2]

result = sisVIVE(Y,D,Z)
plot(result)
}
\keyword{sisVIVE}
