devtools::load_all()
library(dplyr)
library(magrittr)
library(ggplot2)
library(tictoc)
library(xtable)
library(latex2exp)
#### set up data and get FI ####
dat <- matrix(c(40, 176-40, 23, 174-23), nrow=2, byrow=TRUE)
rownames(dat) <- c("FFR", "Angiography")

# get fragility index
bin.fi(dat,alg='exact',test='chisq.prop')

##### try to reproduce their sample size calculation ####
# They said 322
TrialSize::TwoSampleProportion.Equality(.05, .1, .15, .3, 1)*2

#### set up functions for calculation ####
get.p.val <- function(mat) {
  pearson.test.stata <- function(mat) {
    n1 <- sum(mat[1,])
    n2 <- sum(mat[2,])
    p1 <- mat[1,1]/n1
    p2 <- mat[2,1]/n2

    #pbar <- (n1*p1+n2*p2)/(n1+n2)
    #ts <- (p1-p2)/sqrt(pbar*(1-pbar)*(1/n1+1/n2))
    ts <- (p1-p2)/sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)

    p_value <- 2*pnorm(abs(ts), lower.tail=FALSE)
    return(ifelse(is.nan(p_value), 1, p_value))
  }
  suppressWarnings(p <- pearson.test.stata(mat))

  return(ifelse(is.na(p), 1, p))
}

get.replacements <- NULL#list(function(y,x,rn,Y,X) setdiff(c(0,1), y))

my.theta1 <- .15 # angiography
my.delta <- .15
get.sample.null <- function(ss) draw.binom(ss, theta1=my.theta1, theta2=my.theta1, row.prop=1/2, matrix=TRUE)
get.sample.alt <- function(ss) draw.binom(ss, theta1=my.theta1, theta2=my.theta1+my.delta, row.prop=1/2, matrix=TRUE)
get.sample.alt.f <- function(ss, delta) draw.binom(ss, theta1=my.theta1, theta2=my.theta1+delta, row.prop=1/2, matrix=TRUE)

#### set up simulation parameters ####
k <- 6
eps <- .01

phi.grid <- seq(from=0, to=30, by=5)
tau.grid.short <- c(.5, .2)
tau.grid.long <- c(.5, .35, .2)

nsim_low <- 30
nsim_high <- 80

chosen.phi <- 25
chosen.tau <- .5
chosen.pi <- .9

#### run calculation script ####
cl <- parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl, varlist=c("k", "phi.grid", "tau.grid.short", "tau.grid.long", "my.theta1", "my.delta",
                                      "chosen.phi", "chosen.tau", "chosen.pi", "eps", "get.sample.alt.f", "nsim_low",
                                      "nsim_high"), envir=environment())
source('exec/example_simulation_script.R', print.eval=TRUE)
parallel::stopCluster(cl)

#### print output ####
xtable(my_rates_phivary,digits=k-1)
xtable(my_rates_tauvary,digits=k-1)
xtable(trad.out, digits=k-1)
