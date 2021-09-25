\name{vp.postsample}
\alias{vp.postsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{bla
}
\description{
Returns a sample of length 'nsim' from the posterior distribution of the mixing parameters of the VP model: gamma and phi for model (1); gamma, phi, psi1 and psi2 for model (2)}
\usage{
vp.postsample(inla.res, vp.model = "striid", nsim = 10000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inla.res}{
inla output
}
  \item{vp.model}{
character indicating the type of VP model; either 'str' or 'striid'; the former is model (1) in the VP paper, the latter is model (2)
}
  \item{nsim}{
number of samples drawn from the posterior, to pass on to inla.posterior.sample()
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data.frame with 'nsim' rows and 3 columns (tau, gamma and phi, for model (1)), 5 columns (tau, gamma and phi, psi1 and psi2 for model (2))}
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
