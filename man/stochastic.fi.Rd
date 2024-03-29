% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/front.R
\name{stochastic.fi}
\alias{stochastic.fi}
\title{Stochastic generalized fragility indices}
\usage{
stochastic.fi(
  X,
  Y,
  get.replacements,
  get.p.val,
  r = 0.5,
  nsim = 10,
  qfi.init = 10L,
  alpha = 0.05,
  verbose = FALSE,
  cl = NULL,
  D = 20,
  gamma = 0.2,
  init.step = TRUE
)
}
\arguments{
\item{X}{a data frame of covariates which are not subject to modification.}

\item{Y}{a data frame of responses which are subject to modification.}

\item{get.replacements}{a function which outputs a data frame containing all possible row replacements
for Y which are to be considered. The functions inputs the row of Y under consideration,
the row of X under consideration, the row name, and the full original data frames X and Y.}

\item{get.p.val}{a function that inputs X and Y and returns a p value}

\item{r}{the index of the stochastic fragility index, by default 0.5. Having r=0 is equivalent to the generalized fragility
index and having r=1 means that all patient combinations of the output size can reverse significance.}

\item{nsim}{The number of simulations in the root finding algorithm, by default 10}

\item{alpha}{a numeric for the significance cutoff}

\item{verbose}{a logical value for whether to print status updates while running}

\item{cl}{a cluster from the \code{parallel} package, used to compute fragility index over
each modified observation at each stage of the greedy algorithm}

\item{D}{a parameter for Polyak-Ruppert averaging, by default 20}

\item{gamma}{a parameter for Polyak-Ruppert averaging, by default 0.2}

\item{init.step}{a parameter for Polyak-Ruppert averaging, by default TRUE}

\item{gfi.init}{An initialization of the output size, by default 10}
}
\value{
a length 2 list, with the first entry giving the stochastic generalized fragility index and the
last entry giving the history of the root finding algorithm.
}
\description{
This function calculates the stochastic generalized fragility indices, which are a stochastic version of the
generalized fragility indices. They ensure that a random collection of patients can reverse significance with
sufficiently high probability, thus ensuring that typical patients can reverse significance.
}
\examples{
n <- 100
X <- data.frame("tr_group" = sample(c("treated", "not treated"), n, TRUE))
Y <- data.frame("outcome" = sample(c("sick", "healthy"), n, TRUE))
get.p.val <- function(X, Y) fisher.test(table(X[[1]], Y[[1]]))$p.value
get.replacements <- function(y, x, rn, Y, X) data.frame(Y=setdiff(unique(Y[[1]]), y))

stochastic.fi(X, Y, get.replacements, get.p.val)

}
