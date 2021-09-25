\name{vp.plot}
\alias{vp.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variance Partitioning plot}
\description{The Variance Partitioning (VP) plot is a graphical version of the Variance Partitioning table. It summarizes the contribution of the different sources of variation - main, interaction, spatial and temporal effects, etc - in terms of proportion of explained (generalized) variance. This plot is inspired to Gelman (2005) - Analysis of variance: why it is more important than ever. The Annals of Statistics - where anova results are described by a graph which shows the estimated standard deviation for each bunch of random effects in the model. Our VP plot follows the same idea but represents the contribution of each source in a scale (0, 1)}
\usage{
vp.plot(inla.res, vp.model = "striid", table.type = "mixing", cex.sources = 0.95, title.plot = "", nsim = 10000, ...)
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
  \item{cex.sources}{
size of the labels on the left side with the names of the sources (main, int, space, time, etc...)
}
  \item{title.plot}{
character, title of the plot
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
It returns a plot summarizing the contribution of the different sources of variation - main, interaction, spatial and temporal effects, etc - in terms of proportion of explained (generalized) variance
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
