\name{vp.table}
\alias{vp.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variance Partitioning table}
\description{The Variance Partitioning (VP) table summarizes the contribution of the different sources of variation - main, interaction, spatial and temporal effects, etc - in terms of proportion of explained (generalized) variance}
\usage{
vp.table(inla.res, vp.model = "striid", table.type = "mixing", nsim = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inla.res}{
inla output
}
  \item{vp.model}{
character indicating the type of VP model; either 'str' or 'striid'; the former is model (1) in the VP paper, the latter is model (2)
}
  \item{table.type}{
character indicating the type of table, either 'mixing' or 'explained_sd_v1'; mixing is the preferred VP table because it reports the values of the mixing parameters gamma (1-gamma); phi (1-phi); psi1 (1-psi1); psi2 (1-psi2)}
  \item{nsim}{
number of samples drawn from the posterior, to pass on to inla.posterior.sample()
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
It returns a table summarizing the contribution of the different sources of variation - main, interaction, spatial and temporal effects, etc - in terms of proportion of explained (generalized) variance
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
