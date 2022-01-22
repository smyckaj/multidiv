\name{restrict_classe2geosse}
\alias{restrict_classe2geosse}
\title{Convert a set of biogeographic instructions in GeoSSE-like logic into constraints for ClaSSE procedure in diversitree}

\description{
}
\usage{
restrict_classe2geosse(table, forbidden=NULL, single.sympatry=F, no.sympatry=F, single.vicariance=F, no.vicariance=F, pooled.founder=F, single.founder=F,no.founder=F, single.extinction=F, no.extinction=F, pooled.migration=F, single.migration=F, no.migration=F)}
\arguments{
  \item{table}{a data frame coding all present combinations of regions forming geographic ranges}

  \item{forbidden}{a vector of parameter names in ClaSSE convention forcedly set to 0}

  \item{single.sympatry}{set sympatric speciation rate equal for all regions}

  \item{no.sympatry}{set sympatric speciation rates for all regions to 0}

  \item{single.vicariance}{set vicariance speciation rate equal for all region combinations}

  \item{no.sympatry}{set vicariance speciation rates for all region combinations to 0}

  \item{pooled.founder}{set founder speciation rates equal for each group of events with the same target region}

  \item{single.founder}{set founder speciation rate equal for all region combinations}

  \item{no.founder}{set founder speciation rates for all region combinations to 0}

  \item{single.extinction}{set extinction rate equal for all regions}

  \item{no.extinction}{set extinction rates for all region combinations to 0}

  \item{pooled.migration}{set migration rates equal for each group of events with the same target region}

  \item{single.migration}{set migration rate equal for all region combinations}

  \item{no.migration}{set migration rates for all region combinations to 0}
}

\value{
A list of constraints to be used as argument formulae in constrain function of diversitree.
}

\author{
J Smyčka}

\seealso{
\code{\link{diversitree::constrain}}}

\examples{
table=data.frame(A=c(1,0,0,1,1,0,1),B=c(0,1,0,1,0,1,1),C=c(0,0,1,0,1,1,1))
restrict_classe2geosse(table,single.sympatry = T)

}
