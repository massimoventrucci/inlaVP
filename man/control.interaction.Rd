\name{control.interaction}
\alias{control.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set up the interaction model}
\description{
This function takes in input the main effects - defined through m() - and the type of interaction and defines the kronecker product IGMRF for the interaction and the associated constraints. Note; this function works under the assumption the data are sorted in a convenient manner, where both id.space and id.time contain increasing integers and id.time runs faster than id.space.}
\usage{
control.interaction(m1, m2, interaction.type = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{
temporal main effect; it is the output of function m()
}
  \item{m2}{
spatial main effect; it is the output of function m()
}
  \item{interaction.type}{
a number from 1 to 4 (the value 'none', which specifies the additive model, is experimental)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of several objects, the most useful are 'Rkron' the structure matrix of the Kronecker product IGMRF for the interaction model, and 'constr' the constraints for the interaction random effects to pass on to inla.
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
control.interaction(m(1:10, igmrf.type='rw1'),
m(1:20, igmrf.type='rw1'), interaction.type = 4)
}
%control.interaction(m(1:10, igmrf.type='rw1'), m(1:20, igmrf.type='rw1'), interaction.type = 4)
