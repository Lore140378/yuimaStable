\name{setNumerical}
\alias{setNumerical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set numerical details for the Quasi Stable Likelihood Estimation procedure.
}
\description{
\code{setNumerical} is a constructor for \code{\link{yuimaStable.Info}}.
}
\usage{
setNumerical(NumericalHess = FALSE, info = NULL, parallel = T, joint = T, Todorov = F, N = 175, 
   num_of_cores = max(detectCores() - 2, 1), pos = 7, posInt = 1, r0 = 0.9, maxcount = 100, scale = 2, 
   onlypar = FALSE, W = function(x) {1/(1 + x[1:(length(x) - 1)]^2)}, Residual = T, 
   aa_alt = gauss.quad(10000, "laguerre", 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NumericalHess}{ A logical variable, if \code{NumericalHess = FALSE} the asymptotic variance-covariance matrix is analitycally computed using the explicit formula. 
                        By default \code{FALSE}}
  \item{info}{ A list containing additional info. By default an empty \code{list}.}
  \item{parallel}{ A logical variable, if \code{parallel = TRUE}, parallel computing is used to evaluate the density in the qsle procedure. By default \code{TRUE}.}
  \item{joint}{A logical variable, if \code{TRUE} joint quasi likelihood approach is used.}
  \item{Todorov}{A logical variable, if \code{TRUE} \code{beta0} is estimated using Todorov's approach.}
  \item{N}{ Order of polynomials used in both quadratures (Laguerre-Legendre).}
  \item{num_of_cores}{ Number of cores used.}
  \item{pos}{ A positive \code{numeric} object that identifies the subset of the domain where 
      the Gauss-Legendre is applied (\code{-pos <=x < -posInt U posInt < x <= pos}).}
  \item{posInt}{ A positive \code{numeric} object that identifies the subset of the domain where
      the Gauss-Laguerre is applied (\code{-posInt <= x <= posInt}).
}
  \item{r0}{ Parameter for the Todorov's approach.
}
  \item{maxcount}{ Maximum iteration number for the Todorov's approach.}
  \item{scale}{ A numeric object. By default 1. }
  \item{onlypar}{A \code{logical} variable. By default \code{False}, \code{\link{qsle}} also returns the unit-time filtered noise.}
  \item{W}{Weithing function. By default \code{W = 1/(1+x^2)}  }
  \item{Residual}{ A \code{logical} variable.  By default \code{TRUE}.}
  \item{aa_alt}{ A \code{list}.}
}
\value{
A function return an object of \code{\link{yuimaStable.Info}}-class that contains details for all quadratures and procedures applied internally in the \code{\link{qsle}} function.
}
\references{
Insert here url working paper.
}
\author{
Hiroki Masuda and Lorenzo Mercuri.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
