#' Calculate a fragility index for data with binary responses
#'
#' Estimate a fragility index using for the fisher exact test using a table
#'
#' This is a function which is a wrapper around internal functions which calculate
#' the fragility index for different data types and test specifications.
#'
#' @param data.table a 2x2 contingency table, stored as a matrix or table
#' @param X a dataframe representing the covariates
#' @param Y a dataframe representing the responses
#' @param alg a string specifying the FI algorithm, 'greedy' (default), 'exact', or 'walsh'
#' @param test a string specifying the test, defaulting to 'fisher' for the Fisher exact test.
#' An alternative is 'fisher.midp' and 'pearson.chisq'.
#' @param max.f the maximum fragility index considered
#' @param warm.start If using the exact algorithm without specifying a `max.f` argument,
#' a greedy algorithm will be ran first to upper bound the search for the exact algorithm.
#' Set to FALSE to only run the exact algorithm.
#' @param alpha a number for the size of test
#' @param verbose a logical value for whether to print status updates while running
#' @param delta the noninferiority margin for when test='ni.normal'
#' @return the output of greedy.fi for the given test
#'
#' @examples
#' x <- matrix(nrow=2,byrow=TRUE,rgeom(4, 1/50))
#' colnames(x) <- c('event', 'nonevent')
#' rownames(x) <- c('control', 'treatment')
#' out <- bin.fi(data.table=x, alg='exact')
#'
#' @export
bin.fi <- function(data.table=NULL, X=NULL, Y=NULL, alg='greedy', test='fisher', max.f=Inf, alpha=.05,
                   verbose=FALSE, warm.start=TRUE, delta=NULL) {
  if (alg=='greedy' | (max.f==Inf & warm.start)) {
    # convert data.table to X,Y
    if (is.null(Y)) {
      dat <- mat_to_xy(data.table, event.col=2)
      X <- dat[[1]]
      colnames(X) <- 'group'
      Y <- dat[[2]]
      colnames(Y) <- 'resp'
    }
    # get p value function
    if (test=='fisher') {
      get.p.val <- function(X,Y) fisher.test(table(cbind(X, Y)))$p.value
    } else if (test=='fisher.midp') {
      get.p.val <- function(X,Y) {
        m <- sum(Y[[1]]==Y[[1]][1])
        n <- nrow(Y)-m
        k <- sum(X[[1]]==X[[1]][1])
        x <- sum(X[[1]]==X[[1]][1] & Y[[1]]==Y[[1]][1])

        lo <- max(0L, k-n)
        hi <- min(k,m)
        support <- lo:hi
        logdc <- dhyper(support,m,n,k,log=TRUE)
        dnhyper<-function(ncp){
          d<-logdc+log(ncp)*support
          d<-exp(d-max(d))
          d/sum(d)
        }
        d<-dnhyper(1)

        sum(d[d<d[x-lo+1]])+.5*sum(d[d==d[x-lo+1]])
      }
    } else if (test=='pearson.chisq') {
      get.p.val <- function(X,Y) stats::chisq.test(table(cbind(X, Y)))$p.value
    } else if (test=='chisq.prop') {
      get.p.val <- function(X,Y) {
        n1 <- sum(X[[1]]==0)
        n2 <- sum(X[[1]]==1)
        p1 <- sum(X[[1]]==0 & Y[[1]]==1)/n1
        p2 <- sum(X[[1]]==1 & Y[[1]]==1)/n2

        pbar <- (n1*p1+n2*p2)/(n1+n2)
        ts <- (p1-p2)/sqrt(pbar*(1-pbar)*(1/n1+1/n2))
        p_value <- 2*pnorm(abs(ts), lower.tail=FALSE)
        return(ifelse(is.nan(p_value), 1, p_value))
      }
    } else if (test=='ni.normal') {
      get.p.val <- function(X,Y) {
        n1 <- sum(X[[1]]==0)
        n2 <- sum(X[[1]]==1)
        p1 <- sum(X[[1]]==0 & Y[[1]]==1)/n1
        p2 <- sum(X[[1]]==1 & Y[[1]]==1)/n2

        ts <- (p1-p2-delta)/sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
        if (is.nan(ts)) ts <- Inf # if p1=p2=0 or p1=p2=1, then I should reject!

        return(pnorm(ts, lower.tail=FALSE))
      }
    } else {
      stop('Please select an available test option')
    }

    # run alg
    get.replacements <- list(function(y,x,rn,Y,X) setdiff(unique(Y[[1]]), y)) # recall Y is left unchanged in the algorithm
    out <- greedy.fi(X=X, Y=Y, get.replacements=get.replacements, get.p.val=get.p.val, alpha=alpha, verbose=verbose)

    # post process "patients" if started as table?
  }

  if (alg=='exact' || alg=='walsh') {
    # convert X,Y to data.table
    if (is.null(data.table)) {
      data.table <- table(cbind(X, Y))
    }
    # get p value function
    if (test=='fisher') {
      get.p <- function(tab) fisher.test(tab)$p.value
    } else if (test=='fisher.midp') {
      get.p <- function(x) {
        m <- sum(x[,1L])
        n <- sum(x[,2L])
        k <- sum(x[1L,])
        x <- x[1L,1L]

        lo <- max(0L, k-n)
        hi <- min(k,m)
        support <- lo:hi
        logdc <- dhyper(support,m,n,k,log=TRUE)
        dnhyper<-function(ncp){
          d<-logdc+log(ncp)*support
          d<-exp(d-max(d))
          d/sum(d)
        }
        d<-dnhyper(1)

        sum(d[d<d[x-lo+1]])+.5*sum(d[d==d[x-lo+1]])
      }
    } else if (test=='pearson.chisq') {
      get.p <- function(tab) stats::chisq.test(tab)$p.value
    } else if (test=='chisq.prop') {
      get.p <- function(tab) {
        n1 <- sum(tab[1,])
        n2 <- sum(tab[2,])
        p1 <- tab[1,1]/n1
        p2 <- tab[2,1]/n2

        pbar <- (n1*p1+n2*p2)/(n1+n2)
        ts <- (p1-p2)/sqrt(pbar*(1-pbar)*(1/n1+1/n2))
        p_value <- 2*pnorm(abs(ts), lower.tail=FALSE)
        return(ifelse(is.nan(p_value), 1, p_value))
      }
    } else if (test=='ni.normal') {
      get.p <- function(tab) {
        n1 <- sum(tab[1,])
        n2 <- sum(tab[2,])
        p1 <- tab[1,2]/n1
        p2 <- tab[2,2]/n2

        ts <- (p1-p2-delta)/sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
        if (is.nan(ts)) ts <- Inf # if p1=p2=0 or p1=p2=1, then I should reject!
        ts <- unname(ts)

        return(pnorm(ts, lower.tail=FALSE))
      }
    } else {
      stop('Please select an available test option')
    }

    # run alg if exact
    if (alg=='exact') {
      if (warm.start) { # use upper bound as greedy alg to make faster
        if (max.f==Inf) {
          out <- bin.fi.exact(x=data.table, alpha=alpha, get.p=get.p, max.f=abs(out$FI), verbose=verbose)
        } else {
          out <- bin.fi.exact(x=data.table, alpha=alpha, get.p=get.p, max.f=max.f, verbose=verbose)
        }
      } else {
        out <- bin.fi.exact(x=data.table, alpha=alpha, get.p=get.p, max.f=max.f, verbose=verbose)
      }
    }

    # run alg if walsh
    if (alg=='walsh') {
      out <- bin.fi.walsh(mat=data.table, get.p=get.p, alpha=alpha)
    }
  }
  return(out)
}

#' Sample size calculator, taking into account power and fragility index
#'
#' This function takes in a function to compute a p value and a function
#' to simulate data. Using these, it finds the smallest sample size which
#' produces a power and fragility index larger than the input
#' thresholds.
#'
#' @param get.p.val a function that inputs X and Y and returns a p value
#' @param get.replacements a list of of function that inputs a response, covariate, and row name and
#' outputs the replacements to be considered in the algorithm (one function for each response column).
#' If there is only one column, a function can be passed in instead of a list containing the one function.
#' @param get.sample a function which inputs a sample size and outputs a sample, with
#' the covariates in a dataframe and the response in a dataframe
#' @param min.fi the smallest acceptable QUANTILE fragility index. When NULL, the FI calculation is skipped
#' and sample_size_init_fi is taken to produce the desired FI.
#' @param tau the quantile of FI to bound, default 1/2
#' @param min.power the smallest acceptable power. When NULL, the power calculation is skipped and
#' sample_size_init_power is taken to produce the desired power.
#' @param verbose boolean to indicate whether to print out algorithmic information for
#' the sample size calculation
#' @param sample_size_init_power a sample size to initialize the algorithm
#' (not necessary, defaults to 10) to find the sample size for power
#' @param sample_size_init_fi a sample size to initialize the algorithm
#' (not necessary, defaults to the sample size for power) to find the sample size for fi
#' @param alpha a number for the size of test
#' @param nsim the number of simulated draws to consider when estimating the power
#' and expected fragility index
#' @param eps a parameter to control the error. The smaller is it, the more precise the
#' output but the longer the function will take to run.
#' @param algorithm A string specifying the algorithm to use to calculate fragility indices.
#' The default is "greedy"
#' @param gamma the power of n^{-1} in the gradient descent in the Polyak-Ruppert averaging.
#' The default is 0.60.
#' @param niters The number of times to repeat the sample size calculations. The median of the
#' repitions is then reported.
#' @param cl a cluster from the `parallel` package, used for the loop over 1:niters
#' @return the calculated sample sizes for the desired power and fragility index
#'
#' @examples
#' get.p.val <- function(X, Y) {
#'   x <- X[[1]]
#'   y <- Y[[1]]
#'   mat <- matrix(nrow=2,ncol=2)
#'   mat[1,1] <- sum((x==0)&(y==0))
#'   mat[1,2] <- sum((x==0)&(y==1))
#'   mat[2,1] <- sum((x==1)&(y==0))
#'   mat[2,2] <- sum((x==1)&(y==1))
#'   fisher.test(mat)$p.value
#' }
#' get.replacements <- list(function(y,x,rn,Y,X) setdiff(unique(Y[[1]]), y))
#'
#' ss <- general.fi.samplesize(min.fi=10, min.power=.8, get.p.val=get.p.val, niters=1,
#' get.replacements=get.replacements, get.sample=draw.binom, nsim=25, verbose=TRUE)
#'
#' @export
general.fi.samplesize <- function(min.fi=10, min.power=.8,
                                  sample_size_init_power=100L, sample_size_init_fi=NULL,
                                  get.p.val, get.replacements, get.sample, gamma=0.6, niters=50,
                                  cl=NULL, verbose=FALSE, alpha=.05, tau=1/2, nsim=30, eps=.1, algorithm='greedy') {
  # helper functions
  get_power_or_findex <- function(sample_size,getp=TRUE,getfi=FALSE,get.p.val,get.replacements,alpha,verbose,cl=NULL) { #pval and fi of simulated data
    # function which does simulation
    sim_vals <- function(empty) {
      p_vals_i <- NA
      fi_vals_i <- NA

      # get init data
      dat <- get.sample(sample_size)
      if (algorithm=='greedy') {
        X <- dat[[1]]
        Y <- dat[[2]]
      } else if (algorithm=='walsh') {
        mat <- dat
      }

      # get p value
      if (getp) {
        if (algorithm=='greedy') {
          p_vals_i <- get.p.val(X,Y)
        } else if (algorithm=='walsh') {
          p_vals_i <- get.p.val(mat)
        }
      }

      # get fi
      if (getfi) {
        if (algorithm=='greedy') {
          fi_vals_i <- suppressWarnings(greedy.fi(X, Y, get.p.val=get.p.val,
                        get.replacements=get.replacements, alpha=alpha, verbose=FALSE))[[1]]
        } else if (algorithm=='walsh') {
          fi_vals_i <- bin.fi.walsh(mat, get.p.val, alpha)$FI
        } else {
          stop('Please input an appropriate algorithm choice.')
        }# end algorithm conditional
      } # end getfi conditional

      return(c(p_vals_i, fi_vals_i))
    }

    # run the simulation
    if (is.null(cl)) {
      repped <- sapply(1:nsim, sim_vals)
    } else {
      parallel::clusterExport(cl, varlist=c("get.sample", "get.p.val"), envir=environment())
      parallel::clusterExport(cl, varlist=c("draw.binom", "mat_to_xy"), envir=environment())
      repped <- parallel::parSapply(X=1:nsim, FUN=sim_vals, cl=cl)
    }
    p_vals <- repped[1,]
    fi_vals <- repped[2,]

    # stop if bad
    if (any(is.nan(p_vals))) stop('get.p.val returned NaN! Please check edge cases')

    # return output
    if (getp & !getfi) return(mean(p_vals < alpha)) # proportion of times rejecting
    if (!getp & getfi) return(fi_vals) # actual f values
    if (getp & getfi) return(list('p vals'=p_vals, 'fragility indices'=fi_vals)) # all info
  }

  # a secondary helpful function
  get_quantile_fi <- function(sample_size, get.p.val,get.replacements,alpha,verbose,cl=NULL) {
    out <- get_power_or_findex(sample_size=sample_size, getp=FALSE, getfi=TRUE,
                               get.p.val=get.p.val,get.replacements=get.replacements,
                               alpha=alpha,verbose=verbose,cl=NULL)
    counter<-0
    while(min(abs(out))==Inf) {
      out <- get_power_or_findex(sample_size=sample_size, getp=FALSE, getfi=TRUE,
                                 get.p.val=get.p.val,get.replacements=get.replacements,
                                 alpha=alpha,verbose=verbose,cl=NULL)
      counter <- counter+1
      if (counter>10) return((sample_size+1)*sign(out[1]))
    }
    return(quantile(out[abs(out)<Inf], tau))
  }

  # do power based sample size calculation
  if(!is.null(min.power)) {
    if (verbose) print('starting power calculation')
    f_pow <- function(sample_size) get_power_or_findex(sample_size, TRUE, FALSE, get.p.val=get.p.val,get.replacements=get.replacements,alpha=alpha,verbose=verbose,cl=NULL)-min.power

    get_power_zero <- function(empty) find_zero(f_pow, x.init=sample_size_init_power, fz.verbose=FALSE, D=2000, eps=eps, proj=round, gamma=gamma, limits=c(1,9999999))[[1]]
    if (is.null(cl)) {
      ss.s <- sapply(1:niters, FUN=get_power_zero)
    } else {
      parallel::clusterExport(cl, varlist=c("get.sample", "get.p.val"), envir=environment())
      parallel::clusterExport(cl, varlist=c("draw.binom", "find_zero", "bin.fi.walsh"), envir=environment())
      ss.s <- parallel::parSapply(X=1:niters, FUN=get_power_zero, cl=cl)
    }

    sample_size1 <- ceiling(quantile(ss.s, .5, names=FALSE))
  } else {
    if (verbose) print('skipped power calculation')
    sample_size1 <- sample_size_init_power
  }

  # increase for desired average FI among rejections
  if(!is.null(min.fi)) {
    if (is.null(sample_size_init_fi)) sample_size_init_fi <- sample_size1

    if (verbose) print('starting fi calculation')
    f_fi <- function(sample_size) get_quantile_fi(sample_size, get.p.val=get.p.val,get.replacements=get.replacements,alpha=alpha,verbose=verbose,cl=NULL)-min.fi

    get_fi_zero <- function(empty) find_zero(f_fi, x.init=sample_size_init_fi, fz.verbose=FALSE,D=25,eps=eps, proj=round, gamma=gamma, limits=c(1,9999999))[[1]]
    if (is.null(cl)) {
      ss.s <- sapply(1:niters, FUN=get_fi_zero)
    } else {
      parallel::clusterExport(cl, varlist=c("get.sample", "get.p.val"), envir=environment())
      parallel::clusterExport(cl, varlist=c("draw.binom", "find_zero", "bin.fi.walsh"), envir=environment())
      ss.s <- parallel::parSapply(X=1:niters, FUN=get_fi_zero, cl=cl)
    }

    sample_size2 <- ceiling(quantile(ss.s, .5, names=FALSE))
  } else {
    if(verbose) print('skipped fi calculations')
    sample_size2 <- sample_size_init_fi
  }

  # take largest sample size
  sample_size <- c('power_ss'=sample_size1, 'fi_ss'=sample_size2)
  return(sample_size)
}

#' Loops over general.fi.samplesize for a grid of sample sizes.
#'
#' @param min.fi a vector of the smallest acceptable quantile fragility index
#' @param tau the FI quantile, by default 1/2
#' @param get.p.val a function that inputs X and Y and returns a p value
#' @param get.replacements a function that inputs a response, covariate, and row name and outputs the
#' replacements to be considered in the algorithm
#' @param get.sample a function which inputs a sample size and outputs a sample, with
#' the covariates in a dataframe and the response in a dataframe
#' @param cl a cluster from the `parallel` package, used in `greedy.fi`
#' @param min.power the smallest acceptable power
#' @param verbose boolean to indicate whether to print out algorithmic information for the sample size calculation
#' @param sample_size_init_power a sample size to initialize the algorithm
#' (not necessary, defaults to 10) to find the sample size for power
#' @param sample_size_init_fi a sample size to initialize the algorithm
#' (not necessary, defaults to the sample size for power) to find the sample size for fi
#' @param alpha a number for the size of test
#' @param nsim the number of simulated draws to consider when estimating the power
#' and expected fragility index
#' @param eps a parameter to control the error. The smaller is it, the more precise the
#' output but the longer the function will take to run.
#' @param algorithm whether to short circuit the call to greedy.fi and instead use
#' a binary clasic alg. Default is 'greedy'
#' @param gamma the power of n^{-1} in the gradient descent in the Polyak-Ruppert averaging.
#' The default is 0.60.
#' @param niters The number of times to repeat the sample size calculations. The median of the
#' repitions is then reported.
#' @return the calculated sample size for the desired power and fragility index
#'
#' @examples
#' get.p.val <- function(X, Y) {
#'   x <- X[[1]]
#'   y <- Y[[1]]
#'   mat <- matrix(nrow=2,ncol=2)
#'   mat[1,1] <- sum((x==0)&(y==0))
#'   mat[1,2] <- sum((x==0)&(y==1))
#'   mat[2,1] <- sum((x==1)&(y==0))
#'   mat[2,2] <- sum((x==1)&(y==1))
#'   fisher.test(mat)$p.value
#' }
#' get.replacements <- list(function(y,x,rn,Y,X) setdiff(unique(Y[[1]]), y))
#' get.sample <- draw.binom
#'
#' out <- min.fi.curve(min.fi=seq(0, 10, by=5), get.p.val=get.p.val, get.replacements=get.replacements,
#' get.sample=get.sample, niters=1, algorithm='greedy')
#'
#' @export
min.fi.curve <- function(min.fi, get.p.val, get.replacements, get.sample, cl=NULL,
                         min.power=.8, alpha=.05, verbose=FALSE, niters=5,
                         sample_size_init_power=10L, sample_size_init_fi=NULL,
                         nsim=30,eps=.1,tau=1/2,algorithm='greedy',gamma=.6) {
  sample_sizes <- c()
  last_sample_size_fi <- sample_size_init_fi
  last_sample_size_power <- sample_size_init_power
  for (min.fi.val in min.fi) {
    sample_size <- general.fi.samplesize(get.p.val=get.p.val, get.replacements=get.replacements, get.sample=get.sample,
                                         cl=cl, verbose=verbose,
                                         min.fi=min.fi.val, min.power=min.power, alpha=alpha,
                                         sample_size_init_power=last_sample_size_power, sample_size_init_fi=last_sample_size_fi,
                                         nsim=nsim,eps=eps,algorithm=algorithm,tau=tau, gamma=gamma, niters=niters)

    last_sample_size_power <- unname(sample_size['power_ss'])
    last_sample_size_fi <- unname(sample_size['fi_ss'])
    min.power=NULL

    sample_sizes <- rbind(sample_sizes, sample_size)
  }
  rownames(sample_sizes) <- NULL
  return(cbind('min.fi'=min.fi, 'n'=sample_sizes))
}

#' Get rejection rates of the two main tests
#'
#' @param get.p.val a
#' @param get.replacements a
#' @param get.sample.null a set to NULL to skip null calculations
#' @param get.sample.alt a set to NULL to skip alternative calculations
#' @param phi a
#' @param n a
#' @param alpha a
#' @param algorithm a
#' @return a length 4 vector with the size and power of both tests.
#'
#' @examples
#' get.p.val <- function(mat) {
#'   pearson.test.stata <- function(mat) {
#'     n1 <- sum(mat[1,])
#'     n2 <- sum(mat[2,])
#'     p1 <- mat[1,1]/n1
#'     p2 <- mat[2,1]/n2
#'
#'     pbar <- (n1*p1+n2*p2)/(n1+n2)
#'     ts <- (p1-p2)/sqrt(pbar*(1-pbar)*(1/n1+1/n2))
#'     p_value <- 2*pnorm(abs(ts), lower.tail=FALSE)
#'     return(ifelse(is.nan(p_value), 1, p_value))
#'   }
#'   suppressWarnings(p <- pearson.test.stata(mat))
#'
#'   return(ifelse(is.na(p), 1, p))
#' }
#'
#' get.sample.null <- function(ss) draw.binom(ss, theta1=.14, theta2=.14, row.prop=1/2, matrix=TRUE)
#' get.sample.alt <- function(ss) draw.binom(ss, theta1=.14, theta2=.08, row.prop=1/2, matrix=TRUE)
#'
#' get.rejection.rates(get.p.val, NULL, get.sample.null, get.sample.alt, 5, 100, algorithm='walsh')
#'
#' @export
get.rejection.rates <- function(get.p.val, get.replacements=NULL, get.sample.null=NULL, get.sample.alt=NULL,
                                phi, n, alpha=0.05, algorithm='greedy', cl=NULL, nsim=1000) {
  # get function to do simulation
  sim_null_fis <- function(empty) {
    dat <- get.sample.null(n)

    if (algorithm=='walsh') {
      mat <- dat
      bin.fi.walsh(mat, get.p=get.p.val, alpha=alpha)$FI
    } else if (algorithm=='greedy') {
      X <- dat[[1]]
      Y <- dat[[2]]

      greedy.fi(X,Y,get.replacements=get.replacements,get.p.val=get.p.val, alpha=alpha)$FI
    }
  }
  sim_alt_fis <- function(empty) {
    dat <- get.sample.alt(n)

    if (algorithm=='walsh') {
      mat <- dat
      bin.fi.walsh(mat, get.p=get.p.val, alpha=alpha)$FI
    } else if (algorithm=='greedy') {
      X <- dat[[1]]
      Y <- dat[[2]]

      greedy.fi(X,Y,get.replacements=get.replacements,get.p.val=get.p.val, alpha=alpha)$FI
    }
  }

  # init, so that returning works generally
  null.FIs <- NULL
  alt.FIs <- NULL

  # run simulation...depending on parallelizing
  if (is.null(cl)) {
    if (!is.null(get.sample.null)) null.FIs <- sapply(1:nsim, sim_null_fis)
    if (!is.null(get.sample.alt)) alt.FIs <- sapply(1:nsim, sim_alt_fis)
  } else {
    parallel::clusterExport(cl, varlist=c("get.sample.null", "get.sample.alt", "get.replacements", "get.p.val", "n", "alpha", "algorithm"), envir=environment())
    parallel::clusterExport(cl, varlist=c("draw.binom", "greedy.fi", "bin.fi.walsh"), envir=environment()) # setdiff(ls("package:frglty"), c('system.file', 'library.dynam.unload', 'get.rejection.rates'))

    if (!is.null(get.sample.null)) null.FIs <- parallel::parSapply(X=1:nsim, FUN=sim_null_fis, cl=cl)
    if (!is.null(get.sample.alt)) alt.FIs <- parallel::parSapply(X=1:nsim, FUN=sim_alt_fis, cl=cl)
  }

  # return rejection proportion... could be size or power depending on model.size
  return(c('size1'=mean(null.FIs>0), 'power1'=mean(alt.FIs>0),
           'size2'=mean(null.FIs>phi), 'power2'=mean(alt.FIs>phi) #### playing with > vs >=
  ))
}

#' Get equivalent parameters in traditional power calculation
#'
#' @examples
#' get.p.val <- function (mat)
#' {
#'  pearson.test.stata <- function(mat) {
#'    n1 <- sum(mat[1, ])
#'    n2 <- sum(mat[2, ])
#'    p1 <- mat[1, 1]/n1
#'    p2 <- mat[2, 1]/n2
#'    pbar <- (n1 * p1 + n2 * p2)/(n1 + n2)
#'    ts <- (p1 - p2)/sqrt(pbar * (1 - pbar) * (1/n1 + 1/n2))
#'    p_value <- 2 * pnorm(abs(ts), lower.tail = FALSE)
#'    return(ifelse(is.nan(p_value), 1, p_value))
#'  }
#'  suppressWarnings(p <- pearson.test.stata(mat))
#'  return(ifelse(is.na(p), 1, p))
#' }
#' get.sample.alt.f <- function (ss, delta) draw.binom(ss, theta1 = 0.08, theta2 = 0.08+delta,
#' row.prop = 1/2, matrix = TRUE)
#' equivalent_usual_parameters(0.05, .8, .06, 850, get.sample.alt.f, get.p.val=get.p.val)
#'
#' @export
equivalent_usual_parameters <- function(alpha, pi, delta, n,
                                        get.sample.alt.f, get.p.val,
                                        verbose=FALSE, cl=NULL, limits=c(0,1),
                                        nsim=10000, mc.iters=100, delta.iters=100) {
  ## fix n. select 2 of alpha, pi, delta. find the value of the other that gives n
  ## TO DO: put in parallel support...
  p_vals <- replicate(nsim, get.p.val(get.sample.alt.f(n, delta)))

  # get delta
  delta_to_beta <- function(delta.) {
    p.vals <- replicate(mc.iters, get.p.val(get.sample.alt.f(n, delta.)))
    mean(p.vals<alpha)
  }
  delta_vec <- replicate(delta.iters, find_zero(function(delta.) delta_to_beta(delta.)-pi, x.init=delta,
                     limits=limits, D=1/5, eps=0.001, fz.verbose=verbose)$x)

  # return
  c(alpha=quantile(p_vals, pi, names=FALSE), pi=mean(p_vals<alpha), delta=median(delta_vec))
}

#
# get_power_function <- function(delta.grid, phi, alpha=.05, alpha.fi, n,
#                                get.p.val, get.replacements=NULL, get.sample.alt.f,
#                                algorithm='greedy', cl=NULL, nsim=1000) {
#   # assumes that n and alpha.fi are given.
#   rates.p <- rates.fi <- vector(length=length(delta.grid))
#   i <- 1
#   for (delta in delta.grid) {
#     get.sample.alt <- function(ss) get.sample.alt.f(ss, delta)
#     rates.p[i] <- get.rejection.rates(get.p.val=get.p.val, get.replacements=get.replacements, get.sample.null=NULL,
#                                    get.sample.alt=get.sample.alt, phi=0, n=n, alpha=alpha.fi, algorithm=algorithm,
#                                    cl=cl, nsim=nsim)[2]
#     rates.fi[i] <- get.rejection.rates(get.p.val=get.p.val, get.replacements=get.replacements, get.sample.null=NULL,
#                                    get.sample.alt=get.sample.alt, phi=phi, n=n, alpha=alpha, algorithm=algorithm,
#                                    cl=cl, nsim=nsim)[4]
#     i <- i+1
#   }
#   return(cbind(delta=delta.grid, rates.p, rates.fi))
# }


#' Calculate a fragility index using a greedy algorithm
#'
#' Estimate a fragility index using covariates, response(s), and a p value function.
#'
#' This is a general function which uses a greedy algorithm to compute a fragility
#' index.The function arguments accept data from two sources: covariates in X and
#' responses in Y. Covariates are unchanged by the algorithm, while responses are
#' iteratively changed. The type of each response is specified to determine which
#' substitute outcomes to consider.
#'
#' @param X a data frame of covariates
#' @param Y a data frame of responses
#' @param get.replacements a list of of function that inputs a response, covariate, and row name and
#' outputs the replacements to be considered in the algorithm (one function for each response column).
#' If there is only one column, a function can be passed in instead of a list containing the one function.
#' @param get.p.val a function that inputs X and Y and returns a p value
#' @param alpha a number for the size of test
#' @param verbose a logical value for whether to print status updates while running
#' @param cl a cluster from the `parallel` package, used to compute fragility index over
#' each modified observation at each stage of the greedy algorithm
#' @param only.consider a vector of row names to only consider
#' @param dont.consider a vector of row names to not consider
#' @return a list containing a fragility index and related information
#' \describe{
#'   \item{FI}{A fragility index, positive if the initial $p$-value was < alpha and negative otherwise}
#'   \item{p_val_sequence}{An atomic vectors of the p-values at each step of the iteratative algorithm,
#'   where the first is the starting p value and the last is the p value which crossed the alpha
#'   threshold}
#'   \item{reverse}{A boolean indicating whether the original p value is larger than alpha}
#'   \item{num_patients}{The number of patients whose responses were subject to change}
#'   \item{patients}{An atomic vector with the rownames of the rows which had their responses changed}
#'   \item{old_responses}{A list of dataframes of the original responses of each changed patient}
#'   \item{new_responses}{A list of dataframes of the new responses of each changed patient}
#' }
#'
#' @examples
#' n <- 100
#' X <- data.frame('tr_group'=sample(c('treated', 'not treated'), n, TRUE))
#' Y <- data.frame('outcome'=sample(c('sick', 'healthy'), n, TRUE))
#' get.p.val <- function(X,Y) fisher.test(table(X[[1]],Y[[1]]))$p.value
#' get.replacements <- list(function(y,x,rn,Y,X) setdiff(unique(Y[[1]]), y))
#' greedy.fi(X, Y, get.p.val = get.p.val, get.replacements=get.replacements)
#'
#' @export
greedy.fi <- function(X, Y,
                      get.replacements, get.p.val, alpha=0.05,
                      verbose=FALSE, cl=NULL, only.consider=c(), dont.consider=c()) {
  # alpha=.05; verbose=TRUE; only.consider=c(); dont.consider=c();

  # get num_patients
  if (length(only.consider)>0) {
    num_patients <- length(only.consider)
  } else {
    num_patients <- nrow(Y) - length(dont.consider)
  }

  # edit new copies of X and Y so that init versions are stored get_replacements
  XX <- X
  YY <- Y

  # turn get.replacements into list if only a function
  if (!is.list(get.replacements)) get.replacements <- list(get.replacements)

  # loop while same_significance
  starting.p.val <- get.p.val(XX, YY)
  all.p.vals <- c(starting.p.val)
  old.p.val <- 999 # not possible value
  same.significance <- TRUE
  frag.ind.counter <- 0
  changed.patients <- c()
  old.resp.df <- YY[c(),] # empty data frame with same column types as YY
  new.resp.df <- YY[c(),]
  while (same.significance) {
    # get reduced data values to search through
    rows.to.take <- !duplicated(cbind(XX,YY))
    if (length(only.consider)>0) {
      rows.to.take <- rows.to.take & rownames(XX)%in%only.consider
    }
    if (length(dont.consider)>0) {
      rows.to.take <- rows.to.take & !(rownames(XX)%in%dont.consider)
    }
    XX.red <- subset(XX, rows.to.take)
    YY.red <- subset(YY, rows.to.take)

    if (nrow(XX.red)==0) {
      warning("did not converge since dont.consider got too big")
      frag.ind.counter <- Inf
      break
    }

    # get possible replacements list (same length as num of rows of reduced)
    repl <- vector('list', length=nrow(YY.red))
    names(repl) <- rownames(YY.red)
    for (i in 1:length(repl)) {
      # init
      holder <- vector('list', length=ncol(YY.red)) # temporary storage
      names(holder) <- colnames(YY.red)

      ## get the possible responses for each column, stored as list
      for (k in 1:ncol(YY)) {
        holder[[k]] <- get.replacements[[k]](YY.red[[k]][i], XX.red[i,], row.names(YY.red)[i], Y, X)

        ### check if there's no possible responses
        if (length(holder[[k]])==0) stop(paste0('Outcome ', k, ' had no possible alternative values.'))
      }

      ## reorganize to a dataframe with # of columns = # of cols of YY, and number of rows = number of combinations
      holder <- expand.grid(holder, stringsAsFactors=FALSE)

      ## make each row (ie each combination) a separate entry of a list
      repl[[i]] <- split(holder, seq(nrow(holder)))
      for (j in 1:length(repl[[i]])) { # fix the row names to the same as original
        rownames(repl[[i]][[j]]) <- rownames(YY.red)[i]
      }
    } # end for loop over i(repl)

    # init matrix to store output p values
    max.num.resp.options <- max(unlist(lapply(repl, FUN=length)))
    out.p <- matrix(NA, nrow=length(repl), ncol=max.num.resp.options)

    # define function to get p vals
    ## import: repl, YY.red, XX, YY, get.p.val
    get.out.p.row <- function(i) {
      repl.i <- repl[[i]]
      out.p.row <- vector(length=length(repl.i))
      for (j in 1:length(repl.i)) {
        YYY <- YY # just changed this to YYY while copy and pasting through
        YYY[rownames(YY.red)[i],] <- repl.i[[j]]
        out.p.row[j] <- get.p.val(XX, YYY)
      }
      return(out.p.row)
    }

    # get p value for all possible changes
    if (!is.null(cl)) {
      parallel::clusterExport(cl, varlist=c("repl", "YY.red", "XX", "YY", "get.p.val"), envir=environment())
      out.p <- parallel::parSapply(cl=cl, X=1:nrow(YY.red), FUN=get.out.p.row, simplify=FALSE)
    } else {
      out.p <- sapply(X=1:nrow(YY.red), FUN=get.out.p.row, simplify=FALSE)
    }
    out.p <- matrix(unlist(out.p), nrow=length(repl), byrow = TRUE)

    # find best p value
    if (starting.p.val < alpha) {
      best.inds <- arrayInd(which.max(out.p), dim(out.p))[1,]
    } else {
      best.inds <- arrayInd(which.min(out.p), dim(out.p))[1,]
    }
    best.i <- best.inds[1]
    best.j <- best.inds[2]

    # update variables
    my_row <- (rownames(YY)==rownames(YY.red)[best.i])
    current.p.val <- out.p[best.i, best.j]
    all.p.vals <- c(all.p.vals, current.p.val)
    dont.consider <- c(dont.consider, rownames(YY.red)[best.i]) # keeps fragility index definition honest--cant change same person twice, for restricted case
    changed.patients <- c(changed.patients, rownames(YY.red)[best.i])

    frag.ind.counter <- frag.ind.counter+1
    old.resp <- subset(YY, my_row)# YY[my_row,]
    YY[my_row,] <- repl[[best.i]][[best.j]]

    old.resp.df <- rbind(old.resp.df, old.resp)
    new.resp.df <- rbind(new.resp.df, repl[[best.i]][[best.j]])

    if (verbose) print(paste0('FI:', frag.ind.counter,
                              ', patient id:', rownames(XX.red)[best.i],
                              ', p val:', round(get.p.val(XX,YY),3),
                              ', replace:', toString(old.resp),
                              ' with ', toString(YY[rownames(XX.red)[best.i],])))

    # check for stopping
    if (starting.p.val < alpha) {
      if (current.p.val >= alpha) same.significance <- FALSE
    } else {
      if (current.p.val < alpha) same.significance <- FALSE
    }

    # check if spinning wheels (due to only.consider)
    if (old.p.val==current.p.val) {
      warning("did not converge since only.consider was too small: ie, all patients had their response change without altering significance of test")
      frag.ind.counter <- Inf
      break
    }
    old.p.val <- current.p.val
  } ## end while loop

  rev <- (starting.p.val > alpha)
  list('FI'=frag.ind.counter*(-1)^rev, 'p_value_sequence' = unname(all.p.vals),
       'reverse'= unname(rev), 'num_patients'=num_patients,
       'patients'=changed.patients,
       'old_responses'=old.resp.df, 'new_responses'=new.resp.df)
}


#' Simulate right-censored survival responses in two groups
#'
#' @param n the number of observations
#' @param num.times the number of time points at which events or right-censors can occur
#' @return a dataframe with three columns: group, time, and status
#'
#' @examples
#' dat <- get.survival.data(n=100, num.times=6)
#'
#' @export
get.survival.data <- function(n, num.times) {
  # roughly 50/50 in control and treatment
  dat <- data.frame(person_id = 1:n, group=sample(c("control", "treatment"), n, TRUE), stringsAsFactors=FALSE)

  # generate event times depending on group
  all.times <- sort(sample(1:(5*num.times), num.times)) #pick 10 possible event times times
  all.times <- all.times - min(all.times) #makes 0 the first time, for plotting purposes
  p <- list("control" = runif(num.times,0,min(1, 1/(num.times/4))), "treatment" = runif(num.times, 0, min(1, 1/(num.times/2))))
  event.time <- function(group) {
    p <- p[[group]]

    time.of.death <- Inf #for survival censoring at end of study
    for (i in 1:num.times) {
      if (runif(1) < p[i]) {
        time.of.death <- all.times[i]
        break
      }
    }

    censored <- sample(c("yes", "no"), 1, TRUE, c(.1, .9))
    if (time.of.death==Inf) {
      time.of.death <- max(all.times)+1 #edits to remove Inf
      censored = "yes"
    }
    list(censored, time.of.death)
  }

  info <- sapply(dat$group, event.time)
  data.frame('group'=dat[,2],
             "time" = unlist(info[2,]),
             "status" =  (unlist(info[1,])=='no')+0,
             stringsAsFactors=FALSE
  )
}

#' Calculate a fragility index for the coefficient tests for GLMs
#'
#' This is a function which outputs a GLM coefficient table with an extra column
#' which shows the fragility index corresponding to each p-value.
#'
#' @param formula a formula, see the documention of 'glm' for more details
#' @param family a description of the error distribution, see the documention of 'glm' for more details
#' @param data a data frame with measurements on the columns and cases on the rows
#' @param max.step an atomic vector for the max step of each response, for when the response type
#' is restricted
#' @param alpha a number for the size of test
#' @param cl a cluster from the `parallel` package, used in `greedy.fi`
#' @param verbose a logical value for whether to print status updates while running
#' @return the coefficient table from 'glm' with an extra column containing the fragility index of each coefficient
#'
#' @examples
#' data(PimaIndiansDiabetes2, package='mlbench')
#' primdat <- PimaIndiansDiabetes2[complete.cases(PimaIndiansDiabetes2),][1:100,]
#' glm.fi(diabetes ~ ., binomial(), primdat, verbose=FALSE)
#'
#' @export
glm.fi <- function(formula, family, data, max.step=1, alpha=.05, cl=NULL, verbose=FALSE) {
  if (family$family != gaussian()$family & family$family != binomial()$family) stop('Family must be binomial or gaussian.')
  # could possibly generalize to others by making max.step work on eta scale

  if (family$family==binomial()$family) {
    get.replacements <- list(function(y,x,rn,Y,X) setdiff(unique(Y[[1]]), y)) # recall Y is left unchanged in the algorithm
  } else if (family$family==gaussian()$family) {
    get.replacements <- list(function(y,x,rn,Y,X) y + c(-1,1)*max.step) # recall Y is left unchanged in the algorithm
  }

  out <- glm(formula, family, data)
  coef.table <- summary(out)$coefficients

  Z <- model.frame(formula, data=data)
  X <- Z[-1]
  Y <- Z[1]
  fi.coef <- c()
  for (i in 1:nrow(coef.table)) {
    get.p <- function(X, Y) summary(glm(formula, family, cbind(X,Y)))$coefficients[i, 4]

    this.fi <- greedy.fi(X=X, Y=Y, get.replacements=get.replacements, get.p.val=get.p, cl=cl, verbose=verbose, alpha=alpha)
    fi.coef <- c(fi.coef, this.fi$FI) # multiply by sign
  }
  return(cbind(coef.table, 'Fragility index'=fi.coef))
}

#' Calculate a fragility index for data with right-censored survival response
#'
#' This is a function which is a wrapper around internal functions which calculate
#' the fragility index for different data types and test specifications. It is used
#' for survival data.
#'
#' @param time the time of either an event or right censoring
#' @param status a 0,1 variable, with 1 if event and 0 if right-censored
#' @param group a factor with levels representing group membership
#' @param test a string specifying the test type, default is 'logrank'
#' @param max.step a numeric for the maximum allowable change of the time outcome for each
#' patient. Default 0.
#' @param cl a cluster from the `parallel` package, used to compute fragility index over
#' each modified observation at each stage of the greedy algorithm
#' @param alpha a number for the size of test
#' @param min.time the minimum time for replacement times. Default 0.
#' @param max.time the maximum time for replacement times. Default Inf.
#' @param tau an optional parameter for the rmst difference test
#' @param verbose a logical value for whether to print status updates while running
#' @return the output of greedy.fi for the given test
#'
#' @examples
#' dat <- get.survival.data(100, 6)
#' cl <- parallel::makeCluster(parallel::detectCores()-2)
#' out <- surv.fi(dat$time, dat$status, dat$group, max.step=0, cl=cl)
#' parallel::stopCluster(cl)
#'
#' @export
surv.fi <- function(time, status, group, test='logrank', max.step=0, cl=NULL, alpha=.05,
                    min.time=0, max.time=Inf, tau=NULL, verbose=FALSE) {

  Y <- data.frame('time'=time, 'status'=status)
  X <- data.frame('group'=group)

  if (test=='rmst.diff') {
    X$group <- (X$group==X$group[1])*1 # requires group to be 0,1
    get.p.val <- function(X, Y) survRM2::rmst2(Y$time, Y$status, X$group, tau=tau)$unadjusted.result[1,4]
  }
  if (test=='logrank') {
    get.p.val <- function(X,Y) {
      # this commented out way is faster but sometimes gives a wrong answer
      #groups <- levels(as.factor(X$group))
      #p <- fastlogranktest::logrank_test(Y$time[X$group==groups[1]], Y$time[X$group==groups[2]],
      #                  Y$status[X$group==groups[1]], Y$status[X$group==groups[2]])[3]
      #if (p < 10^(-8)) { # logrank_test is faster, but sometimes it rounds too hard
      p <- pchisq(survival::survdiff(survival::Surv(Y$time, Y$status)~X$group, rho=0)$chisq, lower.tail=FALSE, df=1)
      #}
      return(p)
    }
  }

  get.replacements <- list(function(y,x,rn,Y,X) pmin(pmax(min.time, y+c(-1,1)*max.step), max.time), # doesn't let times go negative
                           function(y,x,rn,Y,X) setdiff(unique(Y[[2]]), y)
  )

  out <- greedy.fi(X, Y, get.p.val=get.p.val, get.replacements=get.replacements, cl=cl, verbose=verbose)
  return(out)
}

#' Calculate a fragility index for a meta analysis with binary response
#'
#' This is a function which is a wrapper around internal functions which calculate
#' the fragility index for different data types and test specifications. It is used
#' for meta analyses on binary data. This function makes use of p-value calculators
#' from the `meta` R package.
#'
#' @param data a dataframe with the number of events and observations for the control and treatment
#' groups as columns. The names of the studies must be the rownames, and the names cannot contain two
#' consuctive dots.
#' @param restr.study a character vector giving the names of the studies whose patients are allowed
#' to have response changes. The default is all the studies.
#' @param method A character string indicating which method is to be used for pooling of studies.
#' One of "Inverse", "MH", "Peto", "GLMM", or "SSW", can be abbreviated. (From `meta`.`)
#' @param sm A character string indicating which summary measure ("RR", "OR", "RD", or "ASD").
#' See the `meta`` package documention for more details.
#' @param cl a cluster from the `parallel` package, used to compute fragility index over
#' each modified observation at each stage of the greedy algorithm
#' @param alpha a number for the size of test
#' @param verbose a logical value for whether to print status updates while running
#' @return the output of greedy.fi for the given test
#'
#' @examples
#' data(Olkin95, package='meta')
#' out <- binmeta.fi(data=Olkin95[1:5,-2], restr.study=c('Dewar', 'Lippschutz'))
#'
#' @export
binmeta.fi <- function(data, restr.study=NULL, method='MH', sm='RR', cl=NULL, alpha=.05, verbose=FALSE) {
  #restr.study=NULL; method='MH'; sm='RR'; cl=NULL; alpha=.05; verbose=FALSE
  #data=olkin95

  # convert to X, Y format
  ## remove all unnecessary columns
  ### assume first column is names, 2:5 columns are relative info

  ## get long format df
  row.to.df <- function(dat.row) {
    study.name <- dat.row[[1]]
    dat.row <- dat.row[-1] # remove study name

    dat.row <- unlist(dat.row)

    n <- sum(dat.row[c(2,4)])

    Z <- matrix(ncol=2, nrow=n)
    colnames(Z) <- c("X", "Y")

    row.ind <- 1
    for (i in 1:dat.row[1]) {
      Z[row.ind,] <- c(1,1)
      row.ind <- row.ind+1
    }
    for (i in 1:(dat.row[2]-dat.row[1])) {
      Z[row.ind,] <- c(1,0)
      row.ind <- row.ind+1
    }
    for (i in 1:dat.row[3]) {
      Z[row.ind,] <- c(0,1)
      row.ind <- row.ind+1
    }
    for (i in 1:(dat.row[4]-dat.row[3])) {
      Z[row.ind,] <- c(0,0)
      row.ind <- row.ind+1
    }

    Z <- as.data.frame(Z)
    return(cbind('study'=study.name, Z))
  }
  all.dfs <- plyr::alply(data, 1, row.to.df)
  df <- dplyr::bind_rows(all.dfs)

  ## extract X and Y from df
  X <- subset(df, select=c(1,2))
  Y <- subset(df, select=3)

  # get p value function
  get.p.val <- function(X,Y) {
    df.to.row <- function(df) {
      a <- sum(df$X==1 & df$Y==1)
      b <- sum(df$X==1 & df$Y==0)
      c <- sum(df$X==0 & df$Y==1)
      d <- sum(df$X==0 & df$Y==0)
      return(list('event.e'=a, 'n.e'=a+b, 'event.c'=c, 'n.c'=c+d))
    }
    Z <- cbind(X,Y)

    Z <- data.table::as.data.table(Z)
    dat <- dplyr::summarize(dplyr::group_by(Z, study), 'event.e'=sum(X==1 & Y==1), 'n.e'=sum(X==1), 'event.c'=sum(X==0 & Y==1), 'n.c'=sum(X==0))
    dat <- as.data.frame(dat) #remove data.table class

    p <- meta::metabin(event.e, n.e, event.c, n.c, data=dat, method=method,
                       sm=sm, comb.random=TRUE, comb.fixed=FALSE,
                       overall.hetstat=FALSE)$pval.random
    return(p)
  }

  # setup get.replacements
  get.replacements <- list(function(y,x,rn,Y,X) setdiff(unique(Y[[1]]), y))

  # setup only.consider through restr.study
  if (is.null(restr.study)) {
    only.consider <- c()
  } else {
    only.consider <- rownames(X)[X$study %in% restr.study]
  }

  # get ouput
  out <- greedy.fi(X, Y, get.p.val=get.p.val, get.replacements=get.replacements, only.consider=only.consider, cl=cl, verbose=verbose)
  # what's a good way to handle our not having good patient IDs? I should put in study name at this point.. or into rownames before (in addition to having column)?
  return(out)
}

#' Simulate independent draws from two independent Binomial models
#'
#' @param ss a positive integer representing the sample size
#' @param theta1 the probability of success of the first group
#' @param theta2 the probability of success of the second group
#' @param row.prop the proportion of observations in the first group
#' @param matrix a logical for whether to return a matrix or data frames. Default is FALSE.
#' @return a list containing a (0,1) group indicator data frame and a (0,1) response data frame,
#' each with ss rows and one column
#'
#' @examples
#' out <- draw.binom(100)
#'
#' @export
draw.binom <- function(ss, theta1=.3, theta2=.7, row.prop=1/2, matrix=FALSE) {
  # define model parameters
  m <- floor(ss*row.prop)
  n <- ss-m

  # get count of events in groups
  xx <- rbinom(1, m, theta1)
  yy <- rbinom(1, n, theta2)

  # get dfs if requested
  mat <- matrix(c(xx,yy,m-xx,n-yy),nrow=2)
  rownames(mat) <- 1:2

  if (matrix) {
    return(mat)
  } else {
    return(mat_to_xy(mat))
  }
}
