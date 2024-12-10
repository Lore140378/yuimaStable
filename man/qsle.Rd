\name{qsle}
\alias{qsle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quasi Stable Likelihood estimation for a SDE driven by a local Stable Levy Noise
}
\description{
The function returns the estimators of drift and jump parameters.}
\usage{
qsle(model, start, beta0 = 1.5, data, method = "Nelder-Mead", 
    upper = NULL, lower = NULL, numericalinfo = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ An object of \code{\link{yuimaLocStableSDE}}.}
  \item{start}{ Initial values to be passed to the optimizer.}
  \item{beta0}{ Initial value of \code{beta0} to be passed to the optimizer.}
  \item{data}{  \code{
yuima.data} or \code{\link{zoo}} object that contains the data.}
  \item{method}{Optimization algorithm. See \code{\link{optim}} for details.}
  \item{upper}{Upper bound. See \code{\link{optim}} for details. }
  \item{lower}{Lower bound. See \code{\link{optim}} for details.}
  \item{numericalinfo}{An object of \code{\link{yuimaStable.Info}}. If \code{numericalinfo=NULL} (default value), the inputs for all quadratures and procedures are setted internally. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
## Model Definition##
# Set law
rIG <- function(x,delta,gamma){
  if( delta <= 0 )
    stop("delta must be positive value.")
  if( gamma <= 0 )
    stop("gamma must be positive value.")
  V <- rchisq(x,df=1);
  
  z1 <- ( delta/gamma + V/(2*gamma^2) ) - sqrt( V*delta/(gamma^3) + ( V/(2*gamma^2) )^2 )  
  U <- runif(x,min=0,max=1)
  idx <- which( U < (delta/(delta+gamma*z1)) )
  z2 <- (delta/gamma)^2 /z1[-idx]
  ret <- numeric(x)
  ret[idx] <- z1[idx]
  ret[-idx] <- z2
  
  return(ret)
}

rNIG <- function(n,alpha,delta,t){
  x <-n
  mu<-0
  beta <- 0
  delta <-delta*t
  gamma <- sqrt(alpha^2 - beta^2)
  if (gamma == 0) {
    V = rnorm(x)^2
    Z = delta * delta/V
    X = sqrt(Z) * rnorm(x)
  }else{ 
    Z <- rIG(x,delta,gamma)
    N <- rnorm(x,0,1)
    X <- mu + beta*Z + sqrt(Z)*N
  }
  return(X)
}

mylaw <- setLaw(rng = rNIG, dim = 1)

# SDE definition
model<- setLocStableModel(drift="alpha1*x + alpha2/(1+x^2)",
                  jump.coeff="exp(gam1*cos(x)+gam2*sin(x))",
                  measure = list(df=mylaw))

# Simulation trajectory
h <- 1/500
Deltat<-h/50
FinalT <- 10
truepar <- list(alpha1=-1,alpha2=1,gam1=1.5,
                gam2 =0.5,alpha=5, delta=1)
Samp1 <- setSampling(Initial=0, Terminal=FinalT,n=FinalT/Deltat) # Set Sampling for the approximation a continuous time.
Samp2 <- setSampling(Initial=0, Terminal=FinalT,n=FinalT/h) # Set sampling for the observation.
set.seed(1)
Sim1 <- simulate(object=model, true.parameter = unlist(truepar),
                 sampling = Samp1, subsampling = Samp2)
plot(Sim1)

# Estimation
res <- qsle(model=model, start=truepar, beta0=1.5, data=data)

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
