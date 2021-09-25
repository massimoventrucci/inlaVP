\name{pc.gamma}
\alias{pc.gamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the PC prior density for the mixing parameter 'gamma'
}
\description{
This function computes the density of the PC prior for 'gamma', for a vector of values of 'gamma'; as inputs the user can pass either the scale parameter or the parameters (u,alpha)
}
\usage{
pc.gamma(gamma, lambda, u, alpha, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gamma}{
vector with values of the mixing parameter 'gamma'; values must be inside the (0,1) interval
}
  \item{lambda}{
  the scale parameter of the PC prior
}
  \item{u}{
the upper bound for the mixing parameter 'gamma', from the statement Prob(gamma < U) = alpha
}
  \item{alpha}{
the probability level in Prob(gamma < U) = alpha
}
  \item{log}{
logical; if TRUE the log density is returned (default FALSE)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns the density, or log-density of the PC prior
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
pc.gamma(seq(0, 1, length.out = 10), u=0.95, alpha=0.99, log=TRUE)
}
