\name{pcprior.interaction.lambda}
\alias{pcprior.interaction.lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the scale parameter of the PC prior for the mixing parameter 'gamma'}
\description{Given'U' and 'alpha', this function returns the scale parameter of the PC prior for the mixing parameter 'gamma' balancing main versus interaction; the scale parameter specifies the decay rate from the base model 'gamma=0', this can be deifined by the user through a statement such as Prob(gamma < U) = alpha}
\usage{
pcprior.interaction.lambda(u, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{
the upper bound for the mixing parameter 'gamma', from the statement Prob(gamma < U) = alpha
}
  \item{alpha}{
the probability level in Prob(gamma < U) = alpha
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
It returns 'lambda', the scael parameter of the PC prior for tghe mixing balancing main versus interaction}
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
# Prob(gamma < u) = alpha
pcprior.interaction.lambda(u=0.95, alpha=0.99)
}
