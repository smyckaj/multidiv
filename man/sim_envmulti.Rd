\name{sim_envmulti}
\alias{sim_envmulti}
\title{
Simulating a list of trees with environment-dependent diversification model}
\description{
The function simulates a list of trees with speciation rate exponentially dependent on past environment using sim_env_bd procedure from RPANDA. The null trees are replaced by non-null topologies via inerative calculation. In case parallel package is attached (UNIX only), the simulation can run on multiple CPU cores in parallel.}
\usage{
sim_envmulti(n, env_data, lambda, mu, alpha, time, n.cores = 1)
}
\arguments{
  \item{n}{
an integer specifying number of trees
}
  \item{env_data}{
environmental data, given as a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).}
  \item{lambda}{
a numeric value specifying baseline speciation rate}
  \item{mu}{
a numeric value specifying extinction rate}
  \item{alpha}{
a mumeric value specifying exponential dependence of speciation rate on past environment}
  \item{time}{
a numeric value giving the length of simulation}
  \item{n.cores}{
an integer specifying number of cores for parallel run}
}

\value{
a list of objects of type 'phylo'
}

\author{
J Smyčka}


\examples{
#load dataset of past temperature reconstructions
data(InfTemp)

#simulate 10 pure birth trees with speciation dependent on past temperature with alpha=0.01 and baseline speciation rate lambda=0.1 over 30 time units
tlist=sim_envmulti(10,InfTemp,0.1,0,0.01,30)
tlist

#the same ran across 10 CPU cores
library(parallel)
tlistpar=sim_envmulti(10,InfTemp,0.1,0,0.01,30, n.cores=10)
tlistpar
}
