\name{m}
\alias{m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Define the IGMRF model on the main effects (it builds the structure matrix)}
\description{
This function is intended for internal use, its output is taken in input by control.interaction().}
\usage{
m(x, igmrf.type = "rw1", R = NULL, g, p.spline = FALSE, n.bsplines = 20, deg.bsplines = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector of increasing integer values indicating time or spatial locations (for a temporal main effect, it can be a vector of not equally spaced time locations)
}
  \item{igmrf.type}{
the main effect model, e.g. "rw2", "besag"
}
  \item{R}{
optional; an obj of class 'Cmatrix' with the structure matrix
}
  \item{g}{
only for gmrf.type='besag'; the structure matrix defined through obj 'inla.graph'
}
  \item{p.spline}{
experimental; TRUE of FALSE, whether we want or not to define a P-spline (default FALSE)
}
  \item{n.bsplines}{
experimental; only for p.spline=TRUE, number of b-spline basis function (default is 20)
}
  \item{deg.bsplines}{
experimental; only for p.spline=TRUE, degree of the b-splines (default is 3, cubic b-splines)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
TO DO}
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
x <- 1:50
m1 <- m(x, igmrf.type='rw1')
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
