\name{sim_bdmulti}
\alias{sim_bdmulti}
\title{
Simulating a list of trees with birth-death model}
\description{
The function simulates a list of birth-death trees using rbdsim procedure from ape. In case parallel package is attached (UNIX only), the simulation can run on multiple CPU cores in parallel.
}
\usage{
sim_bdmulti(n, lambda, mu, time, n.cores = 1)
}
\arguments{
  \item{n}{
an integer specifying number of trees
}
  \item{lambda}{
a numeric value specifying speciation rate}
  \item{mu}{
a numeric value specifying extinction rate}
  \item{time}{
a numeric value giving the length of simulation}
  \item{n.cores}{
an integer specifying number of cores for parallel run}
}

\value{
a list of objects of type 'phylo'
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
J Smyčka
}

\examples{
#simulate 10 pure birth trees with speciation rate 0.1 for 30 time units
tlist=sim_bdmulti(10,0.1,0,30)
tlist

#the same ran across 10 CPU cores
library(parallel)
tlistpar=sim_bdmulti(10,0.1,0,30, n.cores = 10)
tlistpar
}
