\name{setLocStableModel}
\alias{setLocStableModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Mathematical description of a stochastic differential  equations driven by a locally Stable Noise
}
\description{
The function returns an object of \code{\link{yuimaLocStableSDE}} that describes an univariate stochastic differential equation 
where the driving noise is locally Stable.
}
\usage{
setLocStableModel(drift, jump.coeff, measure, solve.variable = "x", 
   jump.variable = "z", time.variable = "t", xinit = NULL, 
   Deriv_drift = NULL, Log_Deriv_Jump = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drift}{
a string that defines the drift coefficient.
}
  \item{jump.coeff}{
a string for the jump component.
}
  \item{measure}{
	
Levy measure for jump variables.
}
  \item{solve.variable}{
The name of the left hand side of the SDE.
}
  \item{jump.variable}{
The name of the jump variable.
}
  \item{time.variable}{
The name of the time variable.
}
  \item{xinit}{
The initial condition of SDE. By default 0.
}
  \item{Deriv_drift}{
A function that returns the gradient vector of the drift. By default \code{NULL} 
that implies that the function for the computation of the gradient vector is constructed internally
using symbolic calculus.
}
  \item{Log_Deriv_Jump}{
A function that returns the log derivatives of the jump function. By default \code{NULL} 
that implies that the function for the computation of the log derivatives is constructed internally
using symbolic calculus.
}
}
\details{
Please refer to the vignettes and examples for more details.
}
\value{
The function returns an object of \code{yuimaLocStableSDE}-class that extends the \code{\link{yuima.model}} by introducing four additional slots:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{drift_func }{ that contains a \code{R} function describing the drift of the  SDE.}
\item{jump_funct }{ that contains a \code{R} function describing the jump of the  SDE.}
\item{Deriv_drift }{ that contains a \code{R} function that computes the gradient vector of the drift function.}
\item{Deriv_logjump }{ that contains a \code{R} function that computes the gradient vector of the jump coefficient. }
%% ...
}
\references{
it will be filled as soon as possible.  
}
\author{
Hiroki Masuda and Lorenzo Mercuri
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
# Add examples here.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
