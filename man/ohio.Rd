\name{ohio}
\alias{ohio}
\docType{data}
\title{Lung Cancer Ohio data}
\description{Yearly lung cancer death counts for white males observed from 1968 to 1988 on the 88 counties of Ohio; downloaded from http://www.biostat.umn.edu/~brad/data2}
\usage{data("ohio")}
\format{
  data: A data frame with 1848 rows and the following 4 variables:
\describe{
  \item{County}{an index from 1 to 88 indicating the county}
  \item{Year}{an index from 1 to 21 indicating the province}
  \item{Counts}{number of deaths}
  \item{Population}{population at risk}
}
ohio_graph: an object of class 'inla.graph' with the neighbouring relationships between Ohio counties
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ohio)
}
\keyword{datasets}
