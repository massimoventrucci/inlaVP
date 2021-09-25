\name{covid_italy}
\alias{covid_italy}
\docType{data}
\title{Covid-19 incidence weekly data in Italian provinces.}
\description{
Data cover all of the 107 Italian provinces and span a period of time that goes from the onset of the pandemic on 24-th February 2020 to late July 2021 for a total of 70 weeks; the full dataset is made available by the Italian National Institute of Health through the website https://github.com/ pcm-dpc/COVID-19}
\usage{data("covid_italy")}
\format{
  data: A data frame with 7490 rows and the following 6 variables:
\describe{
  \item{id.week}{an index from 1 to 70 indicating the week which the count refers to}
  \item{id.provinces}{an index from 1 to 107 indicating the province}
  \item{region}{name of the region}
  \item{province}{name of the province}
  \item{new_cases}{number of new cases of covid-19}
  \item{pop_province}{population at risk}
}
italy_graph: an object of class 'inla.graph' with the neighbouring relationships between provinces (3 connected components)
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
data(covid_italy)
}
\keyword{datasets}
