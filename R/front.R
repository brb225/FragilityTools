#' Simulate independent draws from two independent Binomial models
#'
#' This function can be conveniently used in simulations such as a monte carlo
#' sample size calculations.
#'
#' @param ss a positive integer representing the sample size
#' @param theta1 the probability of success of the first group
#' @param theta2 the probability of success of the second group
#' @param row.prop the proportion of observations in the first group. Default is 1/2.
#' @param matrix a logical for whether to return a matrix or data frames. Default is FALSE.
#'
#' @return If the argument matrix is FALSE, a list containing a (0,1) group indicator data frame and a (0,1) response data frame,
#' each with ss rows and one column. If the argument matrix is TRUE, a 2x2 matrix storing the data in a contingency table.
#'
#' @examples
#' out <- draw.binom(100)
#'
#' @export
draw.binom <- function(ss, theta1 = .3, theta2 = .7,
                       row.prop = 1 / 2, matrix = FALSE) {
  # define model parameters
  m <- floor(ss * row.prop)
  n <- ss - m

  # get count of events in groups
  xx <- rbinom(1, m, theta1)
  yy <- rbinom(1, n, theta2)

  # get dfs if requested
  mat <- matrix(c(xx, yy, m - xx, n - yy), nrow = 2)
  rownames(mat) <- 1:2
  colnames(mat) <- c('event', 'nonevent')

  if (matrix) {
    return(mat)
  } else {
    return(mat_to_xy(mat))
  }
}


#' Calculates a fragility index for 2x2 data
#'
#' This is a function which is a wrapper around several general purpose functions to calculate fragility
#' indices. Several algorithms and several methodologies can be employed to find various fragility indices
#' for 2x2 tables. The function supports allowing only sufficiently likely modifications, as described in
#' the article [Incidence fragility index].
#'
#' @param crosstab a 2x2 contingency table, stored as a matrix or table. Either input a data table or both X and Y.
#' @param X a dataframe representing the covariates. Either input a data table or both X and Y.
#' @param Y a dataframe representing the responses. Either input a data table or both X and Y.
#' @param alg a string specifying the FI algorithm, 'exact', 'original', 'original.bothdir', 'greedy', or 'hybrid'.
#' The exact algorithm is described in the [Incidence fragility index] article. It will return the exact fragility
#' index possibly at the cost of longer run time. The original algorithm runs the algorithm proposed by Walsh et al. (2014).
#' The 'original.bothdir' algorithm runs a natural variant of the original algorithm which allows for searching in both directions
#' of the chosen group. The greedy algorithm is described in [Generalized fragility index] and efficiently upper bounds the fragility index.
#' The hybrid approach is described in [Generalized fragility index] and first runs the greedy algorithm and then finetunes the output
#' similar to the exact algorithm.
#' @param test a string specifying the test, defaulting to 'fisher' for the Fisher exact test.
#' Some alternatives are 'fisher.midp', 'pearson.chisq', 'ni.normal'.
#' @param alpha a number for the size of test, default 0.05.
#' @param q the minimum probability of allowable outcome modifications, defaults to 0 (allowing all modifications).
#' Note that alg='original' or alg='original.bothdir' do not support this option.
#' @param verbose a logical value for whether to print status updates while running
#' @param fi.start the starting fragility index considered in the exact algorithm.
#' @param delta the noninferiority margin for when test='ni.normal'
#'
#' @return a list containing the fragility index and other values, depending on which algorithm is specified.
#' The name of the fragility index in the list is 'FI'
#'
#' @examples
#' x <- matrix(nrow = 2, byrow = TRUE, c(5, 100, 30, 70))
#' colnames(x) <- c("event", "nonevent")
#' rownames(x) <- c("treatment", "control")
#' bin.fi(crosstab = x, alg = "greedy", q = .5, verbose = FALSE)$FI
#' bin.fi(crosstab = x, alg = "exact", q = .01, verbose = FALSE)$FI
#'
#' @export
bin.fi <- function(crosstab = NULL,
                   X = NULL, Y = NULL,
                   alg = "greedy", test = "fisher",
                   alpha = .05, q = 0, verbose = FALSE,
                   fi.start = NULL, delta = NULL) {
  data.table <- crosstab # renamed
  min.p <- q # renamed

  # convert X,Y to data.table
  if (is.null(data.table)) {
    data.table <- table(cbind(X, Y))
  }

  # get p value function
  if (test == "fisher") {
    get.p <- function(tab) fisher.test(tab)$p.value
  } else if (test == 'greater.fisher') {
    get.p <- function(tab) fisher.test(tab, alternative='g')$p.value
  } else if (test == 'less.fisher') {
    get.p <- function(tab) fisher.test(tab, alternative='l')$p.value
  } else if (test == 'cambell') {
    get.p <- function(x) {
      r.marg <- rowSums(mat)
      c.marg <- colSums(mat)
      expected.mat <- outer(r.marg, c.marg, '*')/sum(mat)

      if (min(expected.mat)>=1) {
        #N-1 pearson chi square
        #print('n-1 chisq') ##test
        a <- x[1,1]
        b <- x[1,2]
        c <- x[2,1]
        d <- x[2,2]
        N <- sum(x)
        r.marg <- rowSums(x)
        c.marg <- colSums(x)

        STATISTIC <- (a*d-b*c)/prod(r.marg)*(a*d-b*c)/prod(c.marg)*(N-1)
        return(pchisq(STATISTIC, 1, lower.tail = FALSE))
      } else {
        #Irwins 2 sided Fisher exact test
        #print('fisher exact') ##test
        m <- sum(x[, 1L])
        n <- sum(x[, 2L])
        k <- sum(x[1L, ])
        x <- x[1L, 1L]

        lo <- max(0L, k - n)
        hi <- min(k, m)
        support <- lo:hi
        logdc <- dhyper(support, m, n, k, log = TRUE)
        dnhyper <- function(ncp) {
          d <- logdc + log(ncp) * support
          d <- exp(d - max(d))
          d / sum(d)
        }

        relErr <- 1 + 10^(-7)
        d <- dnhyper(1)
        return(sum(d[d <= d[x - lo + 1] * relErr]))
      }
    }
  } else if (test == "fisher.midp") {
    get.p <- function(x) {
      m <- sum(x[, 1L])
      n <- sum(x[, 2L])
      k <- sum(x[1L, ])
      x <- x[1L, 1L]

      lo <- max(0L, k - n)
      hi <- min(k, m)
      support <- lo:hi
      logdc <- dhyper(support, m, n, k, log = TRUE)
      dnhyper <- function(ncp) {
        d <- logdc + log(ncp) * support
        d <- exp(d - max(d))
        d / sum(d)
      }

      d <- dnhyper(1)
      sum(d[d < d[x - lo + 1]]) + .5 * sum(d[d == d[x - lo + 1]])
    }
  } else if (test == "pearson.chisq") {
    get.p <- function(tab) {
      p <- stats::chisq.test(tab)$p.value
      if (is.nan(p)) return(1) # if a marginal is 0, then table is rank one
      return(p)
    }
  } else if (test == "chisq.prop") {
    get.p <- function(tab) {
      n1 <- sum(tab[1, ])
      n2 <- sum(tab[2, ])
      p1 <- tab[1, 1] / n1
      p2 <- tab[2, 1] / n2

      pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
      ts <- (p1 - p2) / sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
      p_value <- 2 * pnorm(abs(ts), lower.tail = FALSE)
      return(ifelse(is.nan(p_value), 1, p_value))
    }
  } else if (test == "ni.normal") {
    get.p <- function(tab) {
      n1 <- sum(tab[1, ])
      n2 <- sum(tab[2, ])
      p1 <- tab[1, 2] / n1
      p2 <- tab[2, 2] / n2

      ts <- (p1 - p2 - delta) / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      if (is.nan(ts)) ts <- Inf # if p1=p2=0 or p1=p2=1, then I should reject!
      ts <- unname(ts)

      return(pnorm(ts, lower.tail = FALSE))
    }
  } else {
    stop("Please select an available test option")
  }

  if (alg=='original') {
    if (min.p != 0) warning('The original algorithm does not support min.p')
    dir <- ifelse(get.p(data.table)<alpha, 'left', 'right')
    out <- bin.fi.walsh(crosstab = data.table, get.p = get.p, alpha = alpha, dir = dir)
  }
  if (alg=='original.bothdir') {
    if (min.p != 0) warning('The original algorithm does not support min.p')
    out <- bin.fi.walsh(crosstab = data.table, get.p = get.p, alpha = alpha, dir = 'both')
  }

  mat.of.probs <- data.table
  mat.of.probs <- mat.of.probs[, c(2, 1)] / apply(mat.of.probs, 1, sum)
  can.vary <- mat.of.probs >= min.p

  if (alg=='greedy') {
    out <- bin.fi.greedy(data.table, get.p = get.p, alpha = alpha, can.vary = can.vary)
  }
  if (alg=='hybrid') {
    out.warmstart <- bin.fi.greedy(data.table, get.p = get.p, alpha = alpha, can.vary = can.vary)
    if (is.finite(out.warmstart$FI)) {
      out <- bin.fi.exact2(crosstab = data.table, get.p = get.p, alpha = alpha, can.vary = can.vary,
                           fi.start = out.warmstart$FI,
                           start.reversed = TRUE,
                           start.p = out.warmstart$p_value_sequence[length(out.warmstart$p_value_sequence)])
    } else { # start at 1 like normal
      out <- bin.fi.exact2(crosstab = data.table, get.p = get.p, alpha = alpha, can.vary = can.vary)
    }
  }
  if(alg == 'exact') {
    if (is.null(fi.start)) fi.start <- 1
    out <- bin.fi.exact2(crosstab = data.table, alpha = alpha, get.p = get.p,
                         fi.start = fi.start, can.vary = can.vary)
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
#' @param min.fi the smallest acceptable QUANTILE fragility index. When NULL, the FI calculation is skipped
#' and sample_size_init_fi is taken to produce the desired FI.
#' @param min.power the smallest acceptable power. When NULL, the power calculation is skipped and
#' sample_size_init_power is taken to produce the desired power.
#' @param sample_size_init_power a sample size to initialize the algorithm
#' (not necessary, defaults to 10) to find the sample size for power
#' @param sample_size_init_fi a sample size to initialize the algorithm
#' (not necessary, defaults to the sample size for power) to find the sample size for fi
#' @param get.p.val a function that inputs X and Y and returns a p value
#' @param get.sample a function which inputs a sample size and outputs a sample, with
#' the covariates in a dataframe X and the response in a dataframe Y
#' @param gamma the power of n^{-1} in the gradient descent in the Polyak-Ruppert averaging.
#' The default is 0.60.
#' @param niters The number of times to repeat the sample size calculations. The median of the
#' repitions is then reported.
#' @param cl a cluster from the `parallel` package, used for the loop over 1:niters
#' @param verbose boolean to indicate whether to print out algorithmic information for
#' the sample size calculation
#' @param alpha a number for the size of test
#' @param tau the quantile of FI to bound, default 1/2
#' @param nsim the number of simulated draws to consider when estimating the power
#' and expected fragility index quantile
#' @param eps a parameter to control the error. The smaller is it, the more precise the
#' output but the longer the function will take to run.
#' @param algorithm A string specifying the algorithm to use to calculate fragility indices.
#' The default is "greedy". Alternatives include "walsh" for the algorithm originally used by Walsh et al. (2014).
#'
#' @return the length two numeric vector of the calculated sample sizes for the desired power and fragility index
#'
#' @examples
#' ss <- general.fi.sscurve(min.fi = 10, min.power = .8, alpha = .05, tau=.5,
#' get.p.val = function(m)fisher.test(m)$p.value, get.sample = function(ss) draw.binom(ss, matrix=TRUE),
#' nsim = 10, niters = 2, verbose = FALSE)
#'
#' @export
general.fi.samplesize <- function(min.fi = 10, min.power = .8,
                                  sample_size_init_power = 100L, sample_size_init_fi = NULL,
                                  get.p.val, get.replacements, get.sample, gamma = 0.6, niters = 50,
                                  cl = NULL, verbose = FALSE, alpha = .05, tau = 1 / 2, nsim = 30,
                                  eps = .1, algorithm = "walsh") {
  # helper function
  get_power_or_findex <- function(sample_size, getp = TRUE, getfi = FALSE, get.p.val, get.replacements, alpha, verbose, cl = NULL) { # pval and fi of simulated data
    # function which does simulation
    sim_vals <- function(empty) {
      p_vals_i <- NA
      fi_vals_i <- NA

      # get init data
      dat <- get.sample(sample_size)
      if (algorithm == "greedy") {
        X <- dat[[1]]
        Y <- dat[[2]]
      } else if (algorithm == "walsh") {
        mat <- dat
      }

      # get p value
      if (getp) {
        if (algorithm == "greedy") {
          p_vals_i <- get.p.val(X, Y)
        } else if (algorithm == "walsh") {
          p_vals_i <- get.p.val(mat)
        }
      }

      # get fi
      if (getfi) {
        if (algorithm == "greedy") {
          fi_vals_i <- suppressWarnings(greedy.fi(X, Y,
                                                  get.p.val = get.p.val,
                                                  get.replacements = get.replacements, alpha = alpha, verbose = FALSE
          ))[[1]]
        } else if (algorithm == "walsh") {
          fi_vals_i <- bin.fi.walsh(mat, get.p.val, alpha)$FI
        } else {
          stop("Please input an appropriate algorithm choice.")
        } # end algorithm conditional
      } # end getfi conditional

      return(c(p_vals_i, fi_vals_i))
    }

    # run the simulation
    if (is.null(cl)) {
      repped <- sapply(1:nsim, sim_vals)
    } else {
      parallel::clusterExport(cl, varlist = c("get.sample", "get.p.val"), envir = environment())
      parallel::clusterExport(cl, varlist = c("draw.binom", "mat_to_xy"), envir = environment())
      repped <- parallel::parSapply(X = 1:nsim, FUN = sim_vals, cl = cl)
    }
    p_vals <- repped[1, ]
    fi_vals <- repped[2, ]

    # stop if bad
    if (any(is.nan(p_vals))) stop("get.p.val returned NaN! Please check edge cases")

    # return output
    if (getp & !getfi) {
      return(mean(p_vals < alpha))
    } # proportion of times rejecting
    if (!getp & getfi) {
      return(fi_vals)
    } # actual f values
    if (getp & getfi) {
      return(list("p vals" = p_vals, "fragility indices" = fi_vals))
    } # all info
  }

  # a secondary helpful function
  get_quantile_fi <- function(sample_size, get.p.val, get.replacements, alpha, verbose, cl = NULL) {
    out <- get_power_or_findex(
      sample_size = sample_size, getp = FALSE, getfi = TRUE,
      get.p.val = get.p.val, get.replacements = get.replacements,
      alpha = alpha, verbose = verbose, cl = NULL
    )
    counter <- 0
    while (min(abs(out)) == Inf) {
      out <- get_power_or_findex(
        sample_size = sample_size, getp = FALSE, getfi = TRUE,
        get.p.val = get.p.val, get.replacements = get.replacements,
        alpha = alpha, verbose = verbose, cl = NULL
      )
      counter <- counter + 1
      if (counter > 10) {
        return((sample_size + 1) * sign(out[1]))
      }
    }
    return(quantile(out[abs(out) < Inf], tau))
  }

  # do power based sample size calculation
  if (!is.null(min.power)) {
    if (verbose) print("starting power calculation")
    f_pow <- function(sample_size) get_power_or_findex(sample_size, TRUE, FALSE, get.p.val = get.p.val, get.replacements = get.replacements, alpha = alpha, verbose = verbose, cl = NULL) - min.power

    get_power_zero <- function(empty) find_zero(f_pow, x.init = sample_size_init_power, fz.verbose = FALSE, D = 2000, eps = eps, proj = round, gamma = gamma, limits = c(1, 9999999))[[1]]
    if (is.null(cl)) {
      ss.s <- sapply(1:niters, FUN = get_power_zero)
    } else {
      parallel::clusterExport(cl, varlist = c("get.sample", "get.p.val"), envir = environment())
      parallel::clusterExport(cl, varlist = c("draw.binom", "find_zero", "bin.fi.walsh"), envir = environment())
      ss.s <- parallel::parSapply(X = 1:niters, FUN = get_power_zero, cl = cl)
    }

    sample_size1 <- ceiling(quantile(ss.s, .5, names = FALSE)) # the points of the niters is to loop over the whole thing to stabilize.. document this better
  } else {
    if (verbose) print("skipped power calculation")
    sample_size1 <- sample_size_init_power
  }

  # increase for desired average FI among rejections
  if (!is.null(min.fi)) {
    if (is.null(sample_size_init_fi)) sample_size_init_fi <- sample_size1

    if (verbose) print("starting fi calculation")
    f_fi <- function(sample_size) get_quantile_fi(sample_size, get.p.val = get.p.val, get.replacements = get.replacements, alpha = alpha, verbose = verbose, cl = NULL) - min.fi

    get_fi_zero <- function(empty) find_zero(f_fi, x.init = sample_size_init_fi, fz.verbose = FALSE, D = 25, eps = eps, proj = round, gamma = gamma, limits = c(1, 9999999))[[1]]
    if (is.null(cl)) {
      ss.s <- sapply(1:niters, FUN = get_fi_zero)
    } else {
      parallel::clusterExport(cl, varlist = c("get.sample", "get.p.val"), envir = environment())
      parallel::clusterExport(cl, varlist = c("draw.binom", "find_zero", "bin.fi.walsh"), envir = environment())
      ss.s <- parallel::parSapply(X = 1:niters, FUN = get_fi_zero, cl = cl)
    }

    sample_size2 <- ceiling(quantile(ss.s, .5, names = FALSE))
  } else {
    if (verbose) print("skipped fi calculations")
    sample_size2 <- sample_size_init_fi
  }

  # take largest sample size
  sample_size <- c("power_ss" = sample_size1, "fi_ss" = sample_size2)
  return(sample_size)
}


#' Fragility index aware sample size calculator for a grid of sample sizes
#'
#' Loops over general.fi.samplesize to calculate the fragility index based sample size.
#' The p value based sample size will only be calculated once.
#'
#' @param min.fi a vector of the smallest acceptable quantile fragility index. Each element will be passed into
#' @param tau the FI quantile, by default 1/2
#' @param algorithm A string, equal to 'walsh' to perform fragility index calculation using the original algorithm
#' due to Walsh et al. (2014) or equal to 'greedy' to perform fragility index calculation using the greedy algorithm
#' described in [Generalized fragility index].
#' @param get.p.val a function that inputs X and Y and returns a p value when algorithm is 'greedy' or that
#' inputs a matrix when algorithm is 'walsh'.
#' @param get.replacements a function which outputs a data frame containing all possible row replacements
#' for Y which are to be considered. The functions inputs the row of Y under consideration,
#' the row of X under consideration, the row name, and the full original data frames X and Y.
#' This gets used in the greedy algorithm.
#' @param get.sample a function which inputs a sample size and outputs a sample, with
#' the covariates in a dataframe and the response in a dataframe
#' @param cl a cluster from the `parallel` package, used in `greedy.fi`
#' @param min.power the smallest acceptable power
#' @param verbose boolean to indicate whether to print out algorithmic information for the sample size calculation
#' @param sample_size_init_power a sample size to initialize the algorithm
#' (not necessary, defaults to 10) to find the sample size for power
#' @param sample_size_init_fi a sample size to initialize the algorithm
#' (not necessary, defaults to the sample size for power) to find the sample size for fi
#' @param alpha a numeric for the significance cutoff
#' @param nsim the number of simulated draws to consider when estimating the power
#' and expected fragility index
#' @param eps a parameter to control the error. The smaller is it, the more precise the
#' output but the longer the function will take to run.
#' @param gamma the power of n^{-1} in the gradient descent in the Polyak-Ruppert averaging.
#' The default is 0.60.
#' @param niters The number of times to repeat the sample size calculations. The median of the
#' repitions is then reported.
#'
#' @return a matrix with 3 columns and row length equal to the length of min.fi. Each row contains the
#' min.fi value and the output of general.fi.samplesize.
#'
#' @examples
#' out <- general.fi.sscurve(min.fi = seq(0, 10, by = 5), get.p.val = function(m) fisher.test(m)$p.value,
#' get.replacements = get.replacements, get.sample = function(ss) draw.binom(ss, matrix=TRUE), niters = 1)
#'
#' @export
general.fi.sscurve <- function(min.fi, get.p.val, get.replacements, get.sample, cl = NULL,
                         min.power = .8, alpha = .05, verbose = FALSE, niters = 5,
                         sample_size_init_power = 10L, sample_size_init_fi = NULL,
                         nsim = 30, eps = .1, tau = 1 / 2, algorithm = "walsh", gamma = .6) {
  sample_sizes <- c()
  last_sample_size_fi <- sample_size_init_fi
  last_sample_size_power <- sample_size_init_power
  for (min.fi.val in min.fi) {
    sample_size <- general.fi.samplesize(
      get.p.val = get.p.val, get.replacements = get.replacements, get.sample = get.sample,
      cl = cl, verbose = verbose,
      min.fi = min.fi.val, min.power = min.power, alpha = alpha,
      sample_size_init_power = last_sample_size_power, sample_size_init_fi = last_sample_size_fi,
      nsim = nsim, eps = eps, algorithm = algorithm, tau = tau, gamma = gamma, niters = niters
    )

    last_sample_size_power <- unname(sample_size["power_ss"])
    last_sample_size_fi <- unname(sample_size["fi_ss"])
    min.power <- NULL

    sample_sizes <- rbind(sample_sizes, sample_size)
  }
  rownames(sample_sizes) <- NULL
  return(cbind("min.fi" = min.fi, "n" = sample_sizes))
}


#' Get equivalent parameters in traditional power calculation
#'
#' This function is designed to help understand the role of the fragility index based sample size calculation in
#' terms of parameters involved in usual p value based sample size calculations. The primary inputs into the function
#' are the usual significance cutoff, power, and effect size, together with the fragility index based sample size. The
#' function then considers the question "how low must the significance cutoff be to reach the input sample size, when all
#' other inputs are fixed at the same value". Then it repeats the question for each input in turn.
#'
#' @param alpha a numeric for the significance cutoff
#' @param pi a numeric for the power of the p value based test
#' @param delta a numeric for the effect size, the difference between the event probability in
#' the first and second group
#' @param n a numeric for the designed sample size
#' @param get.sample.alt.f a function to get a two group sample under the alternative, depending on the sample size
#' and the effect size
#' @param get.p.val a function which inputs a contingency table and outputs a p value
#' @param verbose a boolean for whether to print status updates while running
#' @param cl a cluster from the `parallel` package, not currently supported
#' @param nsim the number of draws used to estimate the distribution of the p value under the alternative, by default 100000
#' @param mc.iters the number of draws used to estimate the distribution of the p value under the alternative in a noisy root finding
#' algorithm, by default 100
#' @param delta.iters the number of iterations used of the root finding algorithm to find the effect size delta which designs the same
#' sample size. The overall output is the median over each of the delta.iters iterations.
#'
#' @return a length 3 numeric vector containing the alpha, pi, and delta for which if the other two are at the input values, the input
#' sample size would be designed.
#'
#' @examples
#' get.p.val <- function(mat) {
#'   pearson.test.stata <- function(mat) {
#'     n1 <- sum(mat[1, ])
#'     n2 <- sum(mat[2, ])
#'     p1 <- mat[1, 1] / n1
#'     p2 <- mat[2, 1] / n2
#'     pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
#'     ts <- (p1 - p2) / sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
#'     p_value <- 2 * pnorm(abs(ts), lower.tail = FALSE)
#'     return(ifelse(is.nan(p_value), 1, p_value))
#'   }
#'   suppressWarnings(p <- pearson.test.stata(mat))
#'   return(ifelse(is.na(p), 1, p))
#' }
#' get.sample.alt.f <- function(ss, delta) {
#'   draw.binom(ss,
#'     theta1 = 0.08, theta2 = 0.08 + delta,
#'     row.prop = 1 / 2, matrix = TRUE
#'   )
#' }
#' get.traditional.ss.params(0.05, .8, .06, 1850, get.sample.alt.f, get.p.val = get.p.val)
#'
#' @export
get.traditional.ss.params <- function(alpha, pi, delta, n,
                                        get.sample.alt.f, get.p.val,
                                        verbose = FALSE, cl = NULL, limits = c(0, 1),
                                        nsim = 10000, mc.iters = 100, delta.iters = 100) {
  ## specify n. select 2 of alpha, pi, delta. find the value of the other that gives n
  ## TO DO: put in parallel support...
  p_vals <- replicate(nsim, get.p.val(get.sample.alt.f(n, delta)))

  # get delta
  delta_to_beta <- function(delta.) {
    p.vals <- replicate(mc.iters, get.p.val(get.sample.alt.f(n, delta.)))
    mean(p.vals < alpha)
  }
  delta_vec <- replicate(delta.iters, find_zero(function(delta.) delta_to_beta(delta.) - pi,
                                                x.init = delta,
                                                limits = limits, D = 1 / 5, eps = 0.001, fz.verbose = verbose
  )$x)

  # return
  c(alpha = quantile(p_vals, pi, names = FALSE), pi = mean(p_vals < alpha), delta = median(delta_vec))
}


#' Get rejection rates of the p value and fragility index based test
#'
#' This function is used to calculate a monte carlo estimate of the size and power of the p value based and
#' fragility index based statistical test.
#'
#' @param get.p.val a function which inputs a contingency table and outputs a p value
#' @param get.sample.null a function which inputs a sample size and outputs a sample under the null hypothesis.
#' By default NULL, which indicates to skip calculations under the null hypothesis
#' @param get.sample.alt a function which inputs a sample size and outputs a sample under the alternative hypothesis.
#' By default NULL, which indicates to skip calculations under the alternative hypothesis
#' @param phi a numeric for the fragility index quantile
#' @param n a numeric for the sample size
#' @param alpha a numeric for the significance cutoff
#' @param algorithm a string indicating the algorithm to use to calculate the fragility index. By default 'greedy'.
#' See bin.fi for more options.
#' @param cl a cluster from the `parallel` package, not currently supported
#' @param nsim a numeric for the number of monte carlo simulations.
#'
#' @return a length 4 vector with the size and power of both tests.
#'
#' @examples
#' get.p.val <- function(mat) {
#'   pearson.test.stata <- function(mat) {
#'     n1 <- sum(mat[1, ])
#'     n2 <- sum(mat[2, ])
#'     p1 <- mat[1, 1] / n1
#'     p2 <- mat[2, 1] / n2
#'
#'     pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
#'     ts <- (p1 - p2) / sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
#'     p_value <- 2 * pnorm(abs(ts), lower.tail = FALSE)
#'     return(ifelse(is.nan(p_value), 1, p_value))
#'   }
#'   suppressWarnings(p <- pearson.test.stata(mat))
#'
#'   return(ifelse(is.na(p), 1, p))
#' }
#'
#' get.sample.null <- function(ss) draw.binom(ss, theta1 = .14, theta2 = .14, row.prop = 1 / 2, matrix = TRUE)
#' get.sample.alt <- function(ss) draw.binom(ss, theta1 = .14, theta2 = .08, row.prop = 1 / 2, matrix = TRUE)
#'
#' get.rejection.rates(get.p.val, get.sample.null, get.sample.alt, phi=2, n=100, algorithm = "walsh")
#'
#' @export
get.rejection.rates <- function(get.p.val, get.sample.null = NULL, get.sample.alt = NULL,
                                phi, n, alpha = 0.05, algorithm = "greedy", cl = NULL, nsim = 1000) {
  # get function to do simulation
  sim_null_fis <- function(empty) {
    dat <- get.sample.null(n)

    if (algorithm == "walsh") {
      mat <- dat
      bin.fi.walsh(mat, get.p = get.p.val, alpha = alpha)$FI
    } else if (algorithm == "greedy") {
      X <- dat[[1]]
      Y <- dat[[2]]

      greedy.fi(X, Y, get.replacements = get.replacements, get.p.val = get.p.val, alpha = alpha)$FI
    }
  }
  sim_alt_fis <- function(empty) {
    dat <- get.sample.alt(n)

    if (algorithm == "walsh") {
      mat <- dat
      bin.fi.walsh(mat, get.p = get.p.val, alpha = alpha)$FI
    } else if (algorithm == "greedy") {
      X <- dat[[1]]
      Y <- dat[[2]]

      greedy.fi(X, Y, get.replacements = get.replacements, get.p.val = get.p.val, alpha = alpha)$FI
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
    parallel::clusterExport(cl, varlist = c("get.sample.null", "get.sample.alt", "get.replacements", "get.p.val", "n", "alpha", "algorithm"), envir = environment())
    parallel::clusterExport(cl, varlist = c("draw.binom", "greedy.fi", "bin.fi.walsh"), envir = environment()) # setdiff(ls("package:frglty"), c('system.file', 'library.dynam.unload', 'get.rejection.rates'))

    if (!is.null(get.sample.null)) null.FIs <- parallel::parSapply(X = 1:nsim, FUN = sim_null_fis, cl = cl)
    if (!is.null(get.sample.alt)) alt.FIs <- parallel::parSapply(X = 1:nsim, FUN = sim_alt_fis, cl = cl)
  }

  # return rejection proportion... could be size or power depending on model.size
  return(c(
    "size1" = mean(null.FIs > 0), "power1" = mean(alt.FIs > 0),
    "size2" = mean(null.FIs > phi), "power2" = mean(alt.FIs > phi) #### playing with > vs >=
  ))
}






















#' Calculate all of the incidence fragility indices
#'
#' The function is used to calculate the incidence fragility indices for every possible sufficiently likely threshold.
#' Note, the incidence fragility indices are formally defined in the article [Incidence fragility index]. The function bin.fi
#' in this package will calculate an incidence fragility index for a given threshold, and this function loops over bin.fi
#' for relevant grid of sufficiently likely thresholds min.p (i.e. q). The function `incidence.plot` allows for a convenient
#' visualization of these values.
#'
#' @param crosstab a 2x2 non negative integer matrix (contingency table) representing the output of a clinical trial.
#' The rows correspond to treatment/control and the columns to event/nonevent.
#' @param alg a string specifying the FI algorithm, 'exact', 'greedy', 'hybrid'.
#' See the documentation of bin.fi for more information. Note that 'original' and 'original.bothdir'
#' do not support the per-patient sufficiently likely threshold so cannot be used here.
#' @param test a string specifying the test, defaulting to 'fisher' for the Fisher exact test.
#' Alternatives include 'fisher.midp' and 'pearson.chisq'.
#' @param verbose a logical value for whether to print status updates while running
#' @param fi.start the starting fragility index considered in the exact algorithm
#' @param alpha a numeric for the significance threshold of the statistical test
#'
#' @return a matrix containing the sufficiently likely threshold values q at which there is
#' a changepoint of the corresponding values of the incidence fragility index
#'
#' @examples
#' x <- matrix(nrow=2,byrow=TRUE,c(5, 100, 30, 70),
#' dimnames = list(c('treatment', 'control'), c('event', 'nonevent')))
#' out <- bin.fi.incidence(crosstab=x, alg='exact')
#'
#' @export
bin.fi.incidence <- function(crosstab = NULL,
                             alg = "exact", test = "fisher",
                             fi.start = NULL, alpha = 0.05, verbose = FALSE,
                             delta = NULL,  X = NULL, Y = NULL) {
  data.table <- crosstab # renamed

  if (alg=='original' | alg=='original.bothdir') stop("Please do not use Walsh's algorithm, as it does
                                                      not support the sufficiently likely threshold")

  # find min.p values
  if (is.null(data.table)) {
    df <- cbind(X, Y)
    df.short <- df[!duplicated(df), ]
    get_prob_from_row <- function(rr) {
      dfr <- df[df[[1]] == rr[[1]], ]
      sum(dfr[[2]] != rr[[2]]) / nrow(dfr)
    }
    minp.grid <- apply(df.short, 1, get_prob_from_row)
  } else if (is.null(X) & is.null(Y)) {
    mat.of.probs <- data.table
    mat.of.probs <- mat.of.probs[, c(2, 1)] / apply(mat.of.probs, 1, sum)
    minp.grid <- c(mat.of.probs)
  } else {
    stop("Please provide either data.table or both X and Y")
  }
  minp.grid <- sort(unique(minp.grid))
  eps <- min(as.numeric(dist(unique(c(0, 1, minp.grid))), method = "manhattan")) / 1000
  minp.grid <- minp.grid + eps
  minp.grid <- c(0, minp.grid)

  fi.grid <- vector(length = length(minp.grid))
  ind <- 1
  for (i in 1:length(minp.grid)) {
    if (verbose) print(paste0('Working on: ', minp.grid[ind]))
    # if (!is.null(fi.start)) fi.start <- 1#ifelse(is.infinite(fi.start), 1, fi.start)
    if (!is.null(fi.start) && is.infinite(fi.start)) {
      fi.grid[ind] <- fi.start
    } else {
      fi.grid[ind] <- suppressWarnings(bin.fi(
        crosstab = data.table, X = X, Y = Y,
        alg = alg, test = test,
        fi.start = fi.start, alpha = alpha, verbose = verbose,
        delta = delta, q = minp.grid[ind]
      )$FI)
    }
    fi.start <- fi.grid[ind]
    ind <- ind + 1
  }
  return(cbind("min.p" = c(0, minp.grid[-1] - eps), "fi" = fi.grid))
}


#' Plot the output of bin.fi.incidence
#'
#' The plot follows the convention that the not filled-in circled points are not on the line
#' and the points with a filled-in circle are on the line.
#'
#' @param out the output of bin.fi.incidence
#' @param ylab a string or LaTeX expression for the vertical axis label
#'
#' @return a plot visualizing the incidence fragility indices
#'
#' @examples
#' x <- matrix(nrow=2,byrow=TRUE,c(5, 100, 30, 70),
#' dimnames = list(c('treatment', 'control'), c('event', 'nonevent')))
#' out <- bin.fi.incidence(crosstab=x, alg='exact')
#' incidence.plot(out)
#'
#' @export
incidence.plot <- function(out, ylab=latex2exp::TeX("Incidence fragility index $FI_q$")) {
  max.fi <- max(out[!is.infinite(out[, 2]), 2])
  min.fi <- min(out[!is.infinite(out[, 2]), 2])
  out <- as.data.frame(out)

  out_seg <- out # segments
  out_seg$minp.right <- c(out_seg$min.p[2:nrow(out_seg)], 1)
  out_seg <- out_seg[!is.infinite(out_seg$fi), ]

  out_fill <- out # right end points
  out_fill$fi <- c(NA, out_fill$fi[1:(nrow(out_fill) - 1)])
  out_fill <- out_fill[!is.infinite(out_fill$fi), ]
  out_fill <- out_fill[-1, ]

  out_unfill <- out # left end points
  out_unfill <- out_unfill[!is.infinite(out_unfill$fi), ]
  out_unfill <- out_unfill[-1, ]

  if (max.fi > 0) {
    plt.lims <- c(0, 1 + max.fi)
  } else {
    plt.lims <- c(min.fi - 1, 0)
  }
  ggplot() +
    geom_segment(aes(x = min.p, y = fi, xend = minp.right, yend = fi), data = out_seg) +
    geom_point(aes(x = min.p, y = fi), data = out_unfill, shape = 1) +
    geom_point(aes(x = min.p, y = fi), data = out_fill) +
    labs(x = latex2exp::TeX("Modification likelihood threshold $q$"),
         y = ylab) + # put in latex2exp::TeX
    xlim(c(0, 1)) +
    scale_y_continuous(limits = plt.lims)
}




















#' LTFU-aware fragility index calculator
#'
#' This function implements a calculation of the LTFU-aware fragility indices described in [LTFU fragility index].
#' These fragility indices consider the potential impact of patients who are lost to follow up.
#'
#'
#' @param crosstab a 2x2 real non-negative integer matrix (contingency table) with control/treatment on the rows and
#' event/nonevent on the columns.
#' @param ltfu a length 2 integer vector counting the number of lost patients in both groups.
#' The first entry corresponds to the first row of mat.
#' @param get.p a function which inputs a 2x2 matrix and returns a p value
#' @param q a numeric vector representing the sufficiently likely threshold. A LTFU-aware fragility
#' index will be calculated for each.
#' @param alpha a number for the significance threshold. The default is 0.05
#' @param sampling.control a list providing the ndraws and mult arguments for sampling. It is NULL by default.
#' When NULL, p_o = p_l is taken by default and the distribution of X_l | X_o is available exactly. ndraws is the
#' number of posterior draws, and mult is the multiplier s (see the LTFU-aware paper).
#' @param xlab a string for the horizontal axis label of the output plot
#' @param ylab a string for the vertical axis label of the output plot
#' @param fig.size A number, by default 1.1. This gives the size of the tiles in the plot. It may need to be
#' manually tuned.
#' @param gradientn.scale A number, by default 0.99. It determines where the tick marks will be for the posterior
#' probability scale. It may need to be manually tuned.
#'
#' @return Returns a list with two entries. The first `FI`, shows the LTFU-aware fragility index
#' for each supplied entry of q, The last row shows the largest possible q and the corresponding
#' LTFU-aware fragility index. The second `info``, shows the outcome modifications that make each
#' LTFU-aware fragility index possible.
#'
#' @examples
#' mat <- matrix(nrow=2, byrow=TRUE, c(315, (1198+314)-315, 350, (1189+328)-350))
#' ltfu <- c(24+3, 27+5)
#' out <- ltfu.fi(mat, ltfu, function(m) fisher.test(m)$p.value)
#'
#' @export
ltfu.fi <- function(crosstab, ltfu, get.p,
                    q=0, alpha=.05,
                    sampling.control=NULL,
                    xlab='Control group LTFU event rate', ylab='Treatment group LTFU event rate',
                    fig.size=1.1, gradientn.scale=.99) {
  mat <- crosstab # renamed

  ltfu.grid <- as.matrix(expand.grid(g1=0:ltfu[1], g2=0:ltfu[2]))
  incidences <- apply(mat, 1, function(v)v[1]/sum(v))
  no <- apply(mat, 1, sum)

  # establish direction to reverse significance
  dir.sig <- get.p(mat)<alpha

  # get posterior probabilities
  if (is.null(sampling.control)) { # use closed form, with s=\infty
    get_pr <- function(rr) {
      m1 <- extraDistr::dbbinom(rr[1], ltfu[1], mat[1,1], no[1]-mat[1,1], log=TRUE) # it would be nice to be able to log
      m2 <- extraDistr::dbbinom(rr[2], ltfu[2], mat[2,1], no[2]-mat[2,1], log=TRUE)
      return(m1+m2)
    }
    prs <- apply(ltfu.grid, 1, get_pr)
    prs <- prs-log(sum(exp(prs))) # the division helps numerical stability here
  } else { # sample with sampling.control...
    #sampling.control <- list(ndraws=100000, mult=1.3)
    ndraws <- sampling.control$ndraws
    mult <- sampling.control$mult
    get_draws <- function(j) {
      Xo <- mat[j,1]
      no <- sum(mat[j,])

      po <- rbeta(ndraws, Xo+0.5, no-Xo+0.5)
      s <- get_beta_parameters(Xo/no, mult)
      pl <- rbeta(ndraws, s*po+1, s*(1-po)+1)
      Xl <- rbinom(ndraws, ltfu[j], pl)
      return(Xl)#/ltfu[j])
    }
    m1 <- get_draws(1)
    m2 <- get_draws(2)

    prs <- data.frame('g1'=m1, 'g2'=m2) # get log pr vector in same order as ltfu.grid
    prs <- dplyr::count(prs, g1, g2)
    prs <- plyr::join(as.data.frame(ltfu.grid), prs, by=c('g1','g2'), type='left')[[3]]
    prs[is.na(prs)] <- 0
    prs <- prs/sum(prs)
    prs <- log(prs)
  }

  # sort by posterior probability
  df <- cbind(ltfu.grid, 'log_prs'=prs)
  df <- df[order(prs, decreasing=TRUE),]
  #ltfu.grid <- ltfu.grid[order(prs, decreasing=TRUE),]
  #prs <- prs[order(prs, decreasing=TRUE)]

  # get best guess (imputation)
  expected <- df[1,1:2]
  #the_best_guess <- wm_prs==(1:nrow(ltfu.grid))

  # get all p values
  get.p.from.new.mat <- function(rr) {
    newmat <- mat + matrix(nrow=2,byrow=TRUE,c(rr[1], ltfu[1]-rr[1], rr[2], ltfu[2]-rr[2]))
    get.p(newmat)
  }
  pvals <- apply(df[,1:2], 1, get.p.from.new.mat)
  df <- cbind(df, 'pval'=pvals)

  # plot
  plt.alpha <- round(alpha, 2)
  maxpr <- signif(exp(max(df[,'log_prs']))*gradientn.scale, 2)
  plt.dat <- data.frame(df, mode=as.numeric(1==(1:nrow(df))))#data.frame(ltfu.grid, p=pvals, mode=the_best_guess, post_pr=prs)
  plt.dat$log_prs <- exp(plt.dat$log_prs)
  plt.dat$sig <- as.factor(plt.dat$pval<alpha)
  levels(plt.dat$sig) <- c(FALSE, TRUE)
  plt <- ggplot(data=plt.dat[order(plt.dat$mode),],
                aes(x=g1/ltfu[1], y=g2/ltfu[2], fill=sig)
  )+
    geom_tile(aes(color=log_prs), size=fig.size)+#, width=ww, height=hh), size=ss)+#width=.182, height=.182), size=1.5)+
    scale_color_gradientn(name='Posterior \nprobability',
                          colours=c("white","blue","black"),
                          breaks=c(maxpr/2, maxpr))+#n.breaks=3)+
    scale_fill_manual(values=c("#999999", "#E69F00"),
                      name="Statistical \nsignificance",
                      labels=c(paste0("p > ", plt.alpha), #>=
                               paste0("p < ", plt.alpha)
                      ),
                      drop=FALSE)+
    labs(x=xlab, y=ylab)+
    theme_bw()+
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position='bottom')

  # posterior probability of being significant
  pp.sig <- unname(sum(exp(df[df[,'pval']<alpha,'log_prs'])))
  #print(paste0('the posterior prob of being significant: ', pp.sig))

  # set up the LTFU-aware FI calculation
  fi.grid <- apply(df[,c('g1', 'g2')], 1, function(rr) sum(abs(unlist(rr)-expected)))
  df <- cbind(df, 'fi'=fi.grid, 'cume_prob'=cumsum(exp(df[,'log_prs'])))

  # do calculation
  sl_lowest_ind_included <- sapply(q, function(qq) max(which(df[,'cume_prob']<=1-qq)))
  get.fi.row <- function(m.ind) {
    df2 <- df[1:m.ind,]
    if (dir.sig) {
      df2 <- df2[df2[,'pval'] >= alpha,]
    } else {
      df2 <- df2[df2[,'pval'] < alpha,]
    }
    if (nrow(df2)==0) return(df2)
    return(df2)
  }

  FI.mat <- matrix(nrow=length(q)+1, ncol=2)
  colnames(FI.mat) <- c('q', 'FI')
  FI.info <- vector(mode='list', length=length(q)+1)
  for (q.ind in 1:length(q)) { # calculate for supplied q
    FI.mat[q.ind,1] <- q[q.ind]
    fi.row <- get.fi.row(sl_lowest_ind_included[q.ind])
    if (nrow(fi.row)==0) {
      FI.mat[q.ind,2] <- NA
      FI.info[[q.ind]] <- NA
      next
    }
    fi.row <- fi.row[fi.row[,'fi']==min(fi.row[,'fi']),,drop=FALSE] # only look at rows with same # of outcome mods
    FI.mat[q.ind,2] <- fi.row[1,'fi']
    FI.info[[q.ind]] <- fi.row
  }
  # then get biggest possible q
  fi.row <- get.fi.row(nrow(df))
  if (nrow(fi.row)==0) { # if cannot reverse significance
    FI.mat[length(q)+1,1:2] <- c(1, NA)
    FI.info[[length(q)+1]] <- NA
  } else {
    FI.mat[length(q)+1,1] <- 1-min(fi.row[,'cume_prob'])
    FI.mat[length(q)+1,2] <- fi.row[1,'fi']
    FI.info[[length(q)+1]] <- fi.row[1,,drop=FALSE] # the minimum
  }
  FI.mat[,'FI'] <- (-1)^(!dir.sig)*FI.mat[,'FI']
  #print(df)

  return(list('FI'=FI.mat, 'info'=FI.info, 'reverse.p'=pp.sig, 'imputation'=expected, 'plot'=plt)) # use the negative if reverse convention
}
































#' Calculate a fragility index using a greedy algorithm
#'
#' This is a very general function which approximately calculates a fragility index
#' using covariates, response(s), and any p value function. It's described in detail
#' in the article [Generalized fragility index].
#'
#' This is a general function which uses a greedy algorithm to compute a fragility
#' index.The function arguments accept data from two sources: covariates in X and
#' responses in Y. Covariates are unchanged by the algorithm, while responses are
#' iteratively changed. The type of each response is specified to determine which
#' substitute outcomes to consider.
#'
#' @param X a data frame of covariates which are not subject to modification.
#' @param Y a data frame of responses which are subject to modification.
#' @param get.replacements a function which outputs a data frame containing all possible row replacements
#' for Y which are to be considered. The functions inputs the row of Y under consideration,
#' the row of X under consideration, the row name, and the full original data frames X and Y.
#' @param get.p.val a function that inputs X and Y and returns a p value
#' @param alpha a numeric for the significance cutoff
#' @param verbose a logical value for whether to print status updates while running
#' @param cl a cluster from the `parallel` package, used to compute fragility index over
#' each modified observation at each stage of the greedy algorithm
#' @param only.consider a vector of row names to only consider, by default NULL which considers all patients
#' @param dont.consider a vector of row names to not consider
#'
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
#'   \item{Zmod}{A data frame of the modified cbind(X,Y)}
#' }
#'
#' @examples
#' n <- 100
#' X <- data.frame("tr_group" = sample(c("treated", "not treated"), n, TRUE))
#' Y <- data.frame("outcome" = sample(c("sick", "healthy"), n, TRUE))
#' get.p.val <- function(X, Y) fisher.test(table(X[[1]], Y[[1]]))$p.value
#' get.replacements <- function(y, x, rn, Y, X) data.frame(setdiff(unique(Y[[1]]), y))
#' greedy.fi(X, Y, get.p.val = get.p.val, get.replacements = get.replacements)$FI
#'
#' @export
greedy.fi <- function(X, Y,
                      get.replacements, get.p.val, alpha = 0.05,
                      verbose = FALSE, cl = NULL,
                      only.consider = NULL, dont.consider = c()) {
  # alpha=.05; verbose=TRUE; only.consider=c(); dont.consider=c();

  # get num_patients
  if(is.null(only.consider)) {
    num_patients <- nrow(Y) - length(dont.consider)
  } else {
    num_patients <- length(setdiff(only.consider, dont.consider))
  }
  # if (length(only.consider) > 0) { # removed on june 8 2021
  #   num_patients <- length(only.consider)
  # } else {
  #   num_patients <- nrow(Y) - length(dont.consider)
  # }

  # fix up only.consider
  if (is.null(only.consider)) only.consider <- rownames(X)

  # edit new copies of X and Y so that init versions are stored get_replacements
  XX <- X
  YY <- Y

  # # turn get.replacements into list if only a function
  # if (!is.list(get.replacements)) get.replacements <- list(get.replacements)

  # loop while same_significance
  starting.p.val <- get.p.val(XX, YY)
  if (verbose) print(paste0("0. Starting p val: ", round(starting.p.val, 3)))
  all.p.vals <- c(starting.p.val)
  old.p.val <- 999 # not possible value
  same.significance <- TRUE
  frag.ind.counter <- 0
  changed.patients <- c()
  old.resp.df <- YY[c(), ] # empty data frame with same column types as YY
  new.resp.df <- YY[c(), ]
  while (same.significance) {
    # get reduced data values to search through
    rows.to.take <- !duplicated(cbind(XX, YY))
    rows.to.take <- rows.to.take & !(rownames(XX) %in% dont.consider) & (rownames(XX) %in% only.consider)
    XX.red <- subset(XX, rows.to.take)
    YY.red <- subset(YY, rows.to.take)

    if (nrow(XX.red) == 0) {
      warning("did not converge since ran out of patients (ie dont.consider got too big)")
      frag.ind.counter <- Inf
      break
    }

    # get possible replacements list (same length as num of rows of reduced)
    repl <- vector("list", length = nrow(YY.red))
    names(repl) <- rownames(YY.red)
    # print(length(repl)) ### ##### ##### #############
    for (i in 1:length(repl)) {
      # # init
      # holder <- vector('list', length=ncol(YY.red)) # temporary storage
      # names(holder) <- colnames(YY.red)
      #
      # ## get the possible responses for each column, stored as list
      # for (k in 1:ncol(YY)) {
      #   holder[[k]] <- get.replacements[[k]](YY.red[[k]][i], XX.red[i,], row.names(YY.red)[i], Y, X)
      #
      #   ### check if there's no possible responses
      #   if (length(holder[[k]])==0) stop(paste0('Outcome ', k, ' had no possible alternative values.'))
      # }
      #
      # ## reorganize to a dataframe with # of columns = # of cols of YY, and number of rows = number of combinations
      # holder <- expand.grid(holder, stringsAsFactors=FALSE)
      #
      # print(holder)
      # print(class(holder))
      # break
      holder <- get.replacements(YY.red[i, ], XX.red[i, ], row.names(YY.red)[i], Y, X)

      # print(i) ###########################################
      # print(holder)

      ## make each row (ie each combination) a separate entry of a list
      repl[[i]] <- split(holder, seq(nrow(holder)))
      # print(repl[[i]]) ####################################################################################
      for (j in 1:length(repl[[i]])) { # fix the row names to the same as original
        rownames(repl[[i]][[j]]) <- rownames(YY.red)[i]
      }
    } # end for loop over i(repl)

    # init matrix to store output p values
    max.num.resp.options <- max(unlist(lapply(repl, FUN = length)))
    out.p <- matrix(NA, nrow = length(repl), ncol = max.num.resp.options)

    # define function to get p vals
    ## import: repl, YY.red, XX, YY, get.p.val
    get.out.p.row <- function(i) {
      repl.i <- repl[[i]]
      out.p.row <- vector(length = length(repl.i))
      for (j in 1:length(repl.i)) {
        if (!any(is.na(repl.i[[j]]))) { # added "any" on jun 29 2021 while debugging surv.fi
          YYY <- YY # just changed this to YYY while copy and pasting through
          YYY[rownames(YY.red)[i], ] <- repl.i[[j]]
          out.p.row[j] <- get.p.val(XX, YYY)
        } else { ## 4/20/2021 added a feature that an NA replacement gives an NA p value
          out.p.row[j] <- NA
        }
      }
      return(out.p.row)
    }

    # get p value for all possible changes
    if (!is.null(cl)) {
      parallel::clusterExport(cl, varlist = c("repl", "YY.red", "XX", "YY", "get.p.val"), envir = environment())
      out.p.lst <- parallel::parSapply(cl = cl, X = 1:nrow(YY.red), FUN = get.out.p.row, simplify = FALSE)
    } else {
      out.p.lst <- sapply(X = 1:nrow(YY.red), FUN = get.out.p.row, simplify = FALSE)
    }

    # turn list of vectors into matrix (ragged, padded with NAs if necessary)
    for (ind in 1:length(out.p.lst)) {
      out.p.lst.ind <- out.p.lst[[ind]]
      out.p.lst.ind <- c(out.p.lst.ind, rep(NA, ncol(out.p)-length(out.p.lst.ind)))
      out.p[ind,] <- out.p.lst.ind
    }
    #out.p <- matrix(unlist(out.p), nrow = length(repl), byrow = TRUE)

    # find best p value
    if (starting.p.val < alpha) { # this function seems to work fine with NA values
      best.inds <- arrayInd(which.max(out.p), dim(out.p))[1, ]
    } else {
      best.inds <- arrayInd(which.min(out.p), dim(out.p))[1, ]
    }
    best.i <- best.inds[1]
    best.j <- best.inds[2]

    # update variables
    my_row <- rownames(YY) == rownames(YY.red)[best.i]
    current.p.val <- out.p[best.i, best.j]
    all.p.vals <- c(all.p.vals, current.p.val)
    dont.consider <- c(dont.consider, rownames(YY.red)[best.i]) # keeps fragility index definition honest--cant change same person twice, for restricted case
    changed.patients <- c(changed.patients, rownames(YY.red)[best.i])

    frag.ind.counter <- frag.ind.counter + 1
    old.resp <- subset(data.frame(XX, YY), my_row) # YY[my_row,]
    old.resp.YY <- subset(YY, my_row)
    # old.resp <- subset(YY, my_row)# YY[my_row,]
    new.resp.YY <- repl[[best.i]][[best.j]]
    YY[my_row, ] <- new.resp.YY#repl[[best.i]][[best.j]]

    old.resp.df <- rbind(old.resp.df, old.resp)
    # new.resp.df <- rbind(new.resp.df, repl[[best.i]][[best.j]])
    new.resp.df <- rbind(new.resp.df, cbind(subset(XX, my_row), repl[[best.i]][[best.j]]))

    if (verbose) {
      print(paste0(
        frag.ind.counter, ".",
        " patient:", rownames(XX.red)[best.i],
        ", p val:", round(get.p.val(XX, YY), 3),
        ", replace:(", toString(old.resp[1, -(1:ncol(X))]),
        ") with (", toString(YY[rownames(XX.red)[best.i], ]),
        ")"
      ))
    }

    # check for stopping
    if (starting.p.val < alpha) {
      if (current.p.val >= alpha) same.significance <- FALSE
    } else {
      if (current.p.val < alpha) same.significance <- FALSE
    }

    # check if spinning wheels
    if (!is.null(only.consider) & frag.ind.counter==length(only.consider)) { # (due to only.consider)
      warning("did not converge since only.consider was too small: ie, all patients had their response change without altering significance of test")
      frag.ind.counter <- Inf
      break
    } # how does this interact with dont consider???
    spinning.oldnew <- TRUE
    for (i.son in 1:dim(old.resp.YY)[2]) {
      #print(old.resp.YY[[i.son]])
      #print(new.resp.YY[[i.son]])
      #print(spinning.oldnew)
      #print(isTRUE(all.equal(old.resp.YY[[i.son]], new.resp.YY[[i.son]])))
      spinning.oldnew <- spinning.oldnew & isTRUE(all.equal(old.resp.YY[[i.son]], new.resp.YY[[i.son]]))
    }
    #print(spinning.oldnew)
    if (spinning.oldnew) { # old_response == new_response
      warning("did not converge: the best outcome modification for each remaining patient was their original outcome")
      frag.ind.counter <- Inf
      break
    }
    # if (old.p.val == current.p.val) { # this is actually too numerically stable
    #   warning("did not converge since only.consider was too small: ie, all patients had their response change without altering significance of test")
    #   frag.ind.counter <- Inf
    #   break
    # }
    old.p.val <- current.p.val
  } ## end while loop

  rev <- (starting.p.val > alpha)
  list(
    "FI" = frag.ind.counter * (-1)^rev, "p_value_sequence" = unname(all.p.vals),
    "reverse" = unname(rev), "num_patients" = num_patients,
    "patients" = changed.patients,
    "old_responses" = old.resp.df, "new_responses" = new.resp.df,
    "Zmod" = cbind(XX,YY) # perhaps remove this later?? or add a print.greedyfi statement
  )
}


#' Calculate a fragility index for survival data with right-censored survival response
#'
#' This is a function which is a wrapper around internal functions which calculate
#' the fragility index for different data types and test specifications. It is used
#' for survival data. The function uses a log rank test by default.
#'
#' @param time the time of either an event or right censoring
#' @param status a 0,1 variable, with 1 if event and 0 if right-censored
#' @param group a factor with levels representing group membership
#' @param test a string specifying the test type, default is 'logrank' but could also use
#' 'rmst.diff' for a restricted mean survival test.
#' @param cl a cluster from the `parallel` package, used to compute fragility index over
#' each modified observation at each stage of the greedy algorithm
#' @param alpha a number for the size of test
#' @param tau an optional parameter for the rmst difference test, by default NULL
#' @param verbose a logical value for whether to print status updates while running
#' @param q the maximum per-patient probability of permitted modifications
#'
#' @return the output of greedy.fi for the time to event test
#'
#' @examples
#' dat <- get.survival.data(100, 6)
#' cl <- parallel::makeCluster(parallel::detectCores() - 2)
#' out <- surv.fi(dat$time, dat$status, dat$group, q=0.3, cl = cl)
#' parallel::stopCluster(cl)
#'
#' @export
surv.fi <- function(time, status, group, test = "logrank", q = 0.5,
                    cl = NULL, alpha = .05,
                    tau = NULL, verbose = FALSE) {
  max.likelihood <- q # renamed
  #if (is.null(max.time)) max.time <- max(time)
  #min.time <- min(time)

  Y <- data.frame("time" = time, "status" = status)
  X <- data.frame("group" = group)

  if (test == "rmst.diff") {
    X$group <- (X$group == X$group[1]) * 1 # requires group to be 0,1
    get.p.val <- function(X, Y) survRM2::rmst2(Y$time, Y$status, X$group, tau = tau)$unadjusted.result[1, 4]
  }
  if (test == "logrank") {
    get.p.val <- function(X, Y) {
      # this commented out way is faster but sometimes gives a wrong answer
      # groups <- levels(as.factor(X$group))
      # p <- fastlogranktest::logrank_test(Y$time[X$group==groups[1]], Y$time[X$group==groups[2]],
      #                  Y$status[X$group==groups[1]], Y$status[X$group==groups[2]])[3]
      # if (p < 10^(-8)) { # logrank_test is faster, but sometimes it rounds too hard
      p <- pchisq(survival::survdiff(survival::Surv(Y$time, Y$status) ~ X$group, rho = 0)$chisq, lower.tail = FALSE, df = 1)
      # }
      return(p)
    }
  }

  # old replacements before forcing get.replacements give a df
  # get.replacements <- list(function(y,x,rn,Y,X) pmin(pmax(min.time, y+c(-1,1)*max.step), max.time), # doesn't let times go negative
  #                          function(y,x,rn,Y,X) setdiff(unique(Y[[2]]), y)
  # )

  # new replacements after forcing get.replacements to give a df
  # get.replacements <- function(y,x,rn,Y,X) {
  #   time.repl <- pmin(pmax(min.time, y$time+c(-1,1)*max.step), max.time)
  #   status.repl <- setdiff(unique(Y[[2]]), y$status)
  #   return(expand.grid(time.repl, status.repl))
  # }

  # new new replacements after updating to new scheme involving median
  # get event preds
  mod.event <- survreg(Surv(time, status) ~ group, dist = "weibull")
  pred.event <- predict(mod.event, type = "quantile", p = 1 / 2)

  # get censor preds
  pred.censor <- tryCatch(
    {
      mod.censor <- survreg(Surv(time, 1 - status) ~ 1, dist = "weibull")
      conv.censor <- TRUE
      predict(mod.censor, type = "quantile", p = 1/2)
    },
    warning = function(w) { # find a better way to handle the warnings and errors...
      stop('ahh weibull did not work')
      conv.censor <- FALSE
      rep(max(Y$time), nrow(Y))
    },
    error = function(e) {
      stop('ahh weibull did not work')
      conv.censor <- FALSE
      rep(max(Y$time), nrow(Y))
    }
  )

  # ensure observed agrees with expected if observed
  pred.event[Y$status == 1] <- Y$time[Y$status == 1]
  pred.censor[Y$status == 0] <- Y$time[Y$status == 0]

  # determine which patients need the below swap
  needed.swap <- rep(FALSE, length(status))
  needed.swap[Y$status == 1 & (pred.censor < pred.event)] <- TRUE
  needed.swap[Y$status == 0 & (pred.event < pred.censor)] <- TRUE

  # ensure that observed outcomes are not before predicted outcomes
  pred.censor[Y$status == 1] <- pmax(pred.event[Y$status == 1], pred.censor[Y$status == 1])
  pred.event[Y$status == 0] <- pmax(pred.event[Y$status == 0], pred.censor[Y$status == 0])

  # some helper functions
  survreg2weib_shape <- function(sr_scale) 1/sr_scale
  survreg2weib_scale <- function(sr_lp) exp(sr_lp) # from ?survreg.distributions
  weib_mode <- function(k, lam) {
    if (k <= 1) return(0)
    else return(lam*((k-1)/k)^(1/k))
  }
  weib_meansd <- function (shape, scale) {
    nu <- 1/shape
    if (nu < 1e-06) {
      mu <- scale * (1 + nu * digamma(1) + nu^2 * (digamma(1)^2 + trigamma(1))/2)
      sigma <- scale^2 * nu^2 * trigamma(1)
    }
    else {
      mu <- gamma(1 + (nu)) * scale
      sigma <- sqrt(gamma(1 + 2 * nu) - (gamma(1 + nu))^2) * scale
    }
    c(mu, sigma)
  }

  get.replacements <- function(y, x, rn, Y, X) { # recall that y = list(time, status)
    # get weibull parameters
    # event
    k.event <- survreg2weib_shape(mod.event$scale)
    lam.event <- survreg2weib_scale(predict(mod.event, type='lp')[rownames(X)==rn])
    # censor
    k.censor <- survreg2weib_shape(mod.censor$scale) # assumes the weibull fit
    lam.censor <- survreg2weib_scale(predict(mod.censor, type='lp')[rownames(X)==rn])

    # get the appropriate outcomes
    y.event <- pred.event[rownames(X)==rn]
    y.censor <- pred.censor[rownames(X)==rn]

    # split based on k, get limits
    q <- sqrt(max.likelihood)
    cdf.event <- pweibull(y.event, k.event, lam.event)
    cdf.censor <- pweibull(y.censor, k.censor, lam.censor)
    # event
    if (k.event <= 1) {
      if (cdf.event > q) { # if do not fill up all the leftern area
        lower.lim.event <- qweibull(cdf.event-q, k.event, lam.event)
        upper.lim.event <- y.event
      } else {
        lower.lim.event <- 0 # will this cause problems?
        upper.lim.event <- qweibull(q, k.event, lam.event)
      }
    } else {
      # find the equal likelihood point
      mode.event <- weib_mode(k.event, lam.event)
      if (y.event == mode.event) {
        el.event <- mode.event
      } else if (y.event > mode.event) { # to the left
        #print(c(0, mode.event))
        el.event <- uniroot(f=function(xx) dweibull(xx, k.event, lam.event)-dweibull(y.event, k.event, lam.event),
                            interval=c(0, mode.event))$root
      } else {
        #print(c(mode.event, 1000))
        msd.event <- weib_meansd(k.event, lam.event)
        el.event <- uniroot(f=function(xx) dweibull(xx, k.event, lam.event)-dweibull(y.event, k.event, lam.event),
                            interval=c(mode.event, msd.event[1]+5*msd.event[2]))$root # hardcoded!!
      }
      cdf.el.event <- pweibull(el.event, k.event, lam.event)

      if (abs(cdf.event - cdf.el.event) >= q) { # actual region is inside the interval found
        if (y.event > mode.event) { # to the left
          lower.lim.event <- qweibull(cdf.event-q, k.event, lam.event)
          upper.lim.event <- y.event
        } else {
          lower.lim.event <- y.event
          upper.lim.event <- qweibull(cdf.event+q, k.event, lam.event)
        }
      } else { # have to do an actual HDR
        hdr.event <- unlist(as.data.frame(stat.extend::HDR.weibull(q, k.event, lam.event))[1:2])
        lower.lim.event <- hdr.event[1]
        upper.lim.event <- hdr.event[2]
      }
    }

    # same for censor
    if (k.censor <= 1) {
      if (cdf.censor > q) { # if do not fill up all the leftern area
        lower.lim.censor <- qweibull(cdf.censor-q, k.censor, lam.censor)
        upper.lim.censor <- y.censor
      } else {
        lower.lim.censor <- 0 # will this cause problems?
        upper.lim.censor <- qweibull(q, k.censor, lam.censor)
      }
    } else {
      # find the equal likelihood point
      mode.censor <- weib_mode(k.censor, lam.censor)
      if (y.censor == mode.censor) {
        el.censor <- mode.censor
      } else if (y.censor > mode.censor) { # cover to the left
        #print(c(0, mode.censor))
        el.censor <- uniroot(f=function(xx) dweibull(xx, k.censor, lam.censor)-dweibull(y.censor, k.censor, lam.censor),
                             interval=c(0, mode.censor))$root
        cvrg.censor <- cdf.censor - pweibull(el.censor, k.censor, lam.censor)
      } else {
        #print(c(mode.censor, 1000))
        msd.censor <- weib_meansd(k.event, lam.event)
        el.censor <- uniroot(f=function(xx) dweibull(xx, k.censor, lam.censor)-dweibull(y.censor, k.censor, lam.censor),
                             interval=c(mode.censor, msd.censor[1]+5*msd.censor[2]))$root # hardcoded!!
      }
      cdf.el.censor <- pweibull(el.censor, k.censor, lam.censor)

      if (abs(cdf.censor - cdf.el.censor) >= q) { # actual region is inside the interval found
        if (y.censor > mode.censor) { # to the left
          lower.lim.censor <- qweibull(cdf.censor-q, k.censor, lam.censor)
          upper.lim.censor <- y.censor
        } else {
          lower.lim.censor <- y.censor
          upper.lim.censor <- qweibull(cdf.censor+q, k.censor, lam.censor)
        }
      } else { # have to do an actual HDR
        hdr.censor <- unlist(as.data.frame(stat.extend::HDR.weibull(q, k.censor, lam.censor))[1:2])
        lower.lim.censor <- hdr.censor[1]
        upper.lim.censor <- hdr.censor[2]
      }
    }

    # put together
    if (!needed.swap[rownames(X)==rn]) { # if did not need to swap out a lower imputation
      # highest and lowest of each component
      if (upper.lim.censor <= upper.lim.event) {
        Cmin <- lower.lim.censor
        Cmax <- upper.lim.censor
      } else if (lower.lim.censor < upper.lim.event & upper.lim.event < upper.lim.censor) {
        Cmin <- lower.lim.censor
        Cmax <- upper.lim.event
      } else {
        Cmin <- NA
        Cmax <- NA
      }
      if (upper.lim.event <= upper.lim.censor) {
        Emin <- lower.lim.event
        Emax <- upper.lim.event
      } else if (lower.lim.event < upper.lim.censor & upper.lim.censor < upper.lim.event) {
        Emin <- lower.lim.event
        Emax <- upper.lim.censor
      } else {
        Emin <- NA
        Emax <- NA
      }
      mods.to.return <- matrix(c(Cmin, Cmax, Emin, Emax,
                                 0, 0, 1, 1),
                               ncol=2)
      colnames(mods.to.return) <- c('time', 'status')
      mods.to.return <- mods.to.return[!is.na(mods.to.return[,1]),]
      mods.to.return <- as.data.frame(mods.to.return)
      mods.to.return <- mods.to.return[!duplicated(mods.to.return),]
    } else { # if this patient did, restrict to only modifying within their observed status
      if (y[[2]]==1) { # event observed
        mods.to.return <- data.frame('time'=c(lower.lim.event, upper.lim.event), 'status'=c(1,1))
      } else { # censor observed
        mods.to.return <- data.frame('time'=c(lower.lim.censor, upper.lim.censor), 'status'=c(0,0))
      }
    }

    # four corners
    # mods.to.return <- matrix(c(lower.lim.event, lower.lim.event, upper.lim.event, upper.lim.event,
    #                            lower.lim.censor, upper.lim.censor, lower.lim.censor, upper.lim.censor),
    #                          ncol=2)
    # mods.to.return <- t(apply(mods.to.return, 1, function(rr) {
    #   c(min(rr), 2-which.min(rr))
    # }))
    # mods.to.return <- as.data.frame(mods.to.return)
    # mods.to.return <- mods.to.return[!duplicated(mods.to.return),]
    return(mods.to.return)
  }

  # get output
  out <- greedy.fi(X, Y, get.p.val = get.p.val, get.replacements = get.replacements, cl = cl, verbose = verbose)
  return(out)

  # # get distribution functions
  # F_event_inv <- function(p, group_name = "control") {
  #   # predict(mod.event, type='quantile', p=p, newdata=list(group=group_name))
  #   if (group_name == "control") {
  #     event.scale <- exp(mod.event$coefficients[1])
  #   } else if (group_name == "treatment") {
  #     event.scale <- exp(sum(mod.event$coefficients))
  #   } else {
  #     stop("aaahhhhh wrong specification of group in internal function!")
  #   }
  #
  #   qweibull(p, scale = event.scale, shape = 1 / mod.event$scale)
  # }
  # F_event <- function(t, group_name = "control") {
  #   if (group_name == "control") {
  #     event.scale <- exp(mod.event$coefficients[1])
  #   } else if (group_name == "treatment") {
  #     event.scale <- exp(sum(mod.event$coefficients))
  #   } else {
  #     stop("aaahhhhh wrong specification of group in internal function!")
  #   }
  #
  #   pweibull(t, scale = event.scale, shape = 1 / mod.event$scale)
  # }
  #
  # if (conv.censor) { # if mod.censor exists use it like mod.event
  #   F_censor_inv <- function(p, group_name = "control") {
  #     # predict(mod.censor, type='quantile', p=p, newdata=list(group=group_name))
  #     if (group_name == "control") {
  #       censor.scale <- exp(mod.censor$coefficients[1])
  #     } else if (group_name == "treatment") {
  #       censor.scale <- exp(sum(mod.censor$coefficients))
  #     } else {
  #       stop("aaahhhhh wrong specification of group in internal function!")
  #     }
  #
  #     qweibull(p, scale = censor.scale, shape = 1 / mod.censor$scale)
  #   }
  #   F_censor <- function(t, group_name = "control") {
  #     if (group_name == "control") {
  #       censor.scale <- exp(mod.censor$coefficients[1])
  #     } else if (group_name == "treatment") {
  #       censor.scale <- exp(sum(mod.censor$coefficients))
  #     } else {
  #       stop("aaahhhhh wrong specification of group in internal function!")
  #     }
  #
  #     pweibull(t, scale = censor.scale, shape = 1 / mod.censor$scale)
  #   }
  # } else { # otherwise assume uniform but for mass at end of study
  #   pr.end.con <- mean(dat$time[dat$group == "control"] == max.time &
  #                        dat$status[dat$group == "control"] == 0)
  #   pr.end.trt <- mean(dat$time[dat$group == "treatment"] == max.time &
  #                        dat$status[dat$group == "treatment"] == 0)
  #
  #   F_censor <- function(t, group_name = "control") {
  #     if (t <= min.time) {
  #       return(0)
  #     }
  #     if (min.time <= t & t <= max.time) {
  #       return(t / max.time * ifelse(group_name == "control", 1 - pr.end.con, 1 - pr.end.trt))
  #     }
  #     if (max.time <= t) {
  #       return(1)
  #     }
  #   }
  #   F_censor_inv <- function(p, group_name = "control") {
  #     if (p < 0 | 1 < p) stop("bad probability")
  #     pr.end <- ifelse(group_name == "control", pr.end.con, pr.end.trt)
  #
  #     if (p <= 1 - pr.end) {
  #       return(p * max.time / (1 - pr.end))
  #     }
  #     if (p >= 1 - pr.end) {
  #       return(max.time)
  #     }
  #   }
  # }
  #
  # F_event_plus_censor_inv <- function(p, group_name = "control") {
  #   uniroot(
  #     function(t) {
  #       F_censor(t, group_name) + F_event(t, group_name) - p
  #     },
  #     c(min(0, min.time), max(max.time, F_censor_inv(p), F_event_inv(p)))
  #   )$root
  # }
  #
  #
  # # get change lookup
  # get_repl_mat <- function(type = "event", extrem = "max", i) {
  #   get_entry <- function(type = "event", extrem = "max", i) {
  #     if (extrem == "max") {
  #       s <- 1
  #     } else if (extrem == "min") {
  #       s <- -1
  #     }
  #     if (type == "event") {
  #       F_ <- F_event
  #       F_inv <- F_event_inv
  #       pred <- pred.event
  #       pred.opp <- pred.censor
  #     } else if (type == "censor") {
  #       F_ <- F_censor
  #       F_inv <- F_censor_inv
  #       pred <- pred.censor
  #       pred.opp <- pred.event
  #     }
  #
  #     # print(min.p) #####
  #     value_extrem_inner <- F_(pred[i], dat$group[i]) + s * min.p
  #     # print(pred[i]) ######
  #     # print(F_(pred[i], dat$group[i])) #####
  #     value_extrem_inner <- max(0, min(1, value_extrem_inner))
  #     value_extrem <- F_inv(value_extrem_inner, dat$group[i])
  #     # print(c(value_extrem, pred.opp[i]))  #####
  #     if (value_extrem <= pred.opp[i]) {
  #       v <- value_extrem
  #     } else {
  #       value_out_inner <- s * min.p + F_event(pred.event[i], dat$group[i]) + F_censor(pred.censor[i], dat$group[i])
  #       # print(value_out_inner) ######
  #       value_out_inner <- max(0, min(1, value_out_inner))
  #       # print(value_out_inner) #######
  #       v <- F_event_plus_censor_inv(value_out_inner, dat$group[i])
  #     }
  #     return(c("time" = v, "status" = as.numeric(type == "event")))
  #   }
  #
  #   rbind(
  #     get_entry("event", "max", i),
  #     get_entry("event", "min", i),
  #     get_entry("censor", "max", i),
  #     get_entry("censor", "min", i)
  #   )
  # }
  #
  # # second take... first checking if set is empty
  # # for if hatE <= hatC
  # uniroot(
  #   function(t) F_event(t) - F_censor(t) - (min.p + F_event(pred.event[1]) - F_censor(pred.censor[1])),
  #   c(pred.event[1], pred.censor[1])
  # )
}


#' Calculate a fragility index for the one sample t test
#'
#' This function returns a fragility index (and accomponying information) for the one sample t test,
#' when each patients outcome is restricted to only be modified by a certain amount. We observe the sufficiently likely
#' convention which is described in the generalized fragility index article.
#'
#' @param y a numeric vector of outcomes
#' @param q a numeric for the minimum probability of outcome changes, by default 1
#' @param mu0 the null mean, by default 0
#' @param alpha a numberic for the significance threshold, by default 0.05
#' @param verbose A boolean for whether to print greedy.fi steps while running,
#' by default FALSE
#' @param cl A parallel cluster for faster calculation in greedy.fi, by default NULL
#'
#' @return The output of greedy.fi (a list) with an additional element which has the
#' overall data (relative) likelihood, as described in the article [Generalized fragility index].
#'
#' @examples
#' set.seed(123456790)
#' y <- rnorm(100, mean = 0)
#' p.grid <- seq(.1, .9, by = .1)
#' fi.grid <- sapply(p.grid, function(p) unlist(ttest.fi(y, q=p)[c('FI', 'sl')]))
#' ggplot2::qplot(p.grid, fi.grid[1,], geom = c('point', "line"), xlim = c(min(p.grid), max(p.grid)),
#'    xlab = "Within-patient Likelihood bound", ylab = "Fragility index", main = "t test fragility indices (n=100)")
#' ggplot2::qplot(log10(fi.grid[2,]), fi.grid[1,], xlab = "Full data (between-patient) Likelihood, log scale",
#'    ylab = "Fragility index", main = "t test fragility indices (n=100)")
#'
#' @export
ttest.fi <- function(y, q = 1, mu0 = 0, alpha = .05,
                     verbose = FALSE, cl=NULL) {
  max.likelihood <- q # renamed

  get.p.val <- function(X, Y) {
    y <- Y[[1]]
    if (any(!is.finite(y))) {
      2 * pt(-1, df = length(y) - 1) # from calculation of limit
    } else {
      t.test(y, mu = mu0)$p.value
    }
  }

  mu.est <- mean(y)
  sigma.est <- sd(y)

  # define get.replacements
  #if (outcome.alg == 'hdr') {
  get.replacements <- function(y, x, rn, Y, X) {
    Fy <- pnorm((y-mu.est)/sigma.est)#pt((y-mu.est)/sigma.est, nrow(Y)-1)
    if (max.likelihood <= 1-2*min(Fy, 1-Fy)) { # one end point is y
      if (y < mu.est) {
        new.outcomes <- mu.est+sigma.est*qnorm(Fy+max.likelihood)#qt(Fy+max.likelihood, nrow(Y)-1)
      } else {
        new.outcomes <- mu.est+sigma.est*qnorm(Fy-max.likelihood)#qt(Fy-max.likelihood, nrow(Y)-1)
      }
      new.outcomes <- c(y, new.outcomes)
    } else { # actual HDR
      new.outcomes <- mu.est+sigma.est*(qt(0.5+c(-1,1)*max.likelihood/2, nrow(Y)-1))
    }

    yred <- Y[[1]][row.names(Y)!=rn]
    n <- length(Y[[1]])

    sm <- sum(yred)
    c1 <- sm/n - mu0
    c2 <- -sm/n
    c3 <- sum((yred - sm/n)^2)
    c4 <- -2*sm/n^2
    o <- (2*c1*c2*n^2 - 2*c1*c2*n + c1*c4*n^2 - 2*c2^2*n - 2*c3*n)/(-2*c1*n^2 + 2*c1*n + 2*c2*n - 2*c2 + c4*n) # the extremum of the test statistic
    if (min(new.outcomes) < o & o < max(new.outcomes)) {
      new.outcomes <- c(new.outcomes, o)
    } else {
      new.outcomes <- c(new.outcomes, NA)
    }
    # I should figure out whether o is a max or min.. then narrow down new.outcomes accordingly

    return(data.frame(new.outcomes))
  }
  #}
  # if (outcome.alg == "data") {
  #   left.probs <- pnorm(y, mu.est, sigma.est)
  #   if (is.null(max.likelihood)) max.likelihood <- 1 - .99 * min(c(left.probs, 1 - left.probs))
  #
  #   change.lookup <- matrix(
  #     ncol = 2, byrow = FALSE,
  #     c(
  #       "left" = qnorm(pmax(0, left.probs - (1 - max.likelihood)), mu.est, sigma.est),
  #       "right" = qnorm(pmin(1, left.probs + (1 - max.likelihood)), mu.est, sigma.est)
  #     )
  #   )
  #   get.replacements <- function(y, x, rn, Y, X) {
  #     mat <- change.lookup[rownames(Y) == rn, ]
  #     as.data.frame(mat)
  #   }
  # }
  # if (outcome.alg == "order") {
  #   if (is.null(max.likelihood)) stop("Please specify max.likelihood if using the order option.")
  #
  #   # estimate order statistic distribution
  #   n.sims <- 100000
  #   OS <- matrix(rnorm(n.sims * length(y), mu.est, sigma.est), nrow = n.sims)
  #   OS <- t(apply(OS, 1, sort)) # gets order statistics in the rows
  #   OS <- apply(OS, 2, sort) # gets order statistics in increasing order
  #
  #   # get position of y relative to order statistics
  #   y <- sort(y)
  #   inds <- c()
  #   for (i in 1:length(y)) {
  #     inds[i] <- which.min(abs(OS[, i] - y[i]))
  #   }
  #
  #   # get change lookup table
  #   change.lookup <- matrix(ncol = 2, nrow = length(y))
  #   for (i in 1:length(y)) {
  #     change.lookup[i, 1] <- OS[max(1, inds[i] - (1 - max.likelihood) * (n.sims - 1)), i]
  #     change.lookup[i, 2] <- OS[min(nrow(OS), inds[i] + (1 - max.likelihood) * (n.sims - 1)), i]
  #   }
  #   get.replacements <- function(y, x, rn, Y, X) as.data.frame(change.lookup[rownames(Y) == rn, ])
  # }

  X <- data.frame("constant_covariate" = rep(0, length(y)))
  Y <- data.frame("outcome" = y)
  out <- greedy.fi(
    X = X, Y = Y,
    get.replacements = get.replacements, get.p.val = get.p.val,
    alpha = alpha, verbose = verbose, cl=cl
  )
  out$FI <- (-1)^(get.p.val(X,Y)>=alpha)*out$FI
  out[[length(out) + 1]] <- max.likelihood
  names(out)[length(out)] <- "max.likelihood"

  # get sufficiently likely info
  selected.inds <- as.numeric(rownames(out$old_responses))
  y.mod <- y
  y.mod[selected.inds] <- out$new_responses[,2]
  out$sl <- normal.sl(y, y.mod)
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
#' consuctive dots. The order should be event.e, n.e, event.c, n.c, as in meta::binmeta
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
#' @param q the sufficiently likely threshold, defaults to 0, permitting all modifications. Any
#' empirical conditional outcome probability strictly lower than q will not be permitted to
#' receive outcome modifications.
#'
#' @return the output of greedy.fi specialized to the meta analysis case
#'
#' @examples
#' data(Olkin95, package = "meta")
#' out <- binmeta.fi(data = Olkin95[1:5, -2], restr.study = c("Dewar", "Lippschutz"))
#' out$FI
#'
#' @export
binmeta.fi <- function(data, restr.study = NULL, method = "MH", sm = "RR",
                       cl = NULL, alpha = .05, verbose = FALSE, q = 0) {
  # restr.study=NULL; method='MH'; sm='RR'; cl=NULL; alpha=.05; verbose=FALSE
  # data=olkin95

  # convert to X, Y format
  ## remove all unnecessary columns
  ### assume first column is names, 2:5 columns are relative info

  ## get long format df
  row.to.df <- function(dat.row) {
    study.name <- dat.row[[1]]
    dat.row <- dat.row[-1] # remove study name

    dat.row <- unlist(dat.row)

    n <- sum(dat.row[c(2, 4)])

    Z <- matrix(ncol = 2, nrow = n)
    colnames(Z) <- c("X", "Y")
    #rownames(Z) <- 1:n

    row.ind <- 1
    for (i in 1:dat.row[1]) {
      if (dat.row[1]==0) break
      Z[row.ind, ] <- c(1, 1)
      row.ind <- row.ind + 1
    }
    for (i in 1:(dat.row[2] - dat.row[1])) {
      if (dat.row[2] - dat.row[1]==0) break
      Z[row.ind, ] <- c(1, 0)
      row.ind <- row.ind + 1
    }
    for (i in 1:dat.row[3]) {
      if (dat.row[3]==0) break
      Z[row.ind, ] <- c(0, 1)
      row.ind <- row.ind + 1
    }
    for (i in 1:(dat.row[4] - dat.row[3])) {
      if (dat.row[4] - dat.row[3]==0) break
      Z[row.ind, ] <- c(0, 0)
      row.ind <- row.ind + 1
    }

    Z <- as.data.frame(Z)
    return(cbind("study" = study.name, Z))
  }
  all.dfs <- plyr::alply(data, 1, row.to.df)
  df <- dplyr::bind_rows(all.dfs)

  ## extract X and Y from df
  X <- subset(df, select = c(1, 2))
  Y <- subset(df, select = 3)

  # get p value function
  get.p.val <- function(X, Y) {
    df.to.row <- function(df) {
      a <- sum(df$X == 1 & df$Y == 1)
      b <- sum(df$X == 1 & df$Y == 0)
      c <- sum(df$X == 0 & df$Y == 1)
      d <- sum(df$X == 0 & df$Y == 0)
      return(list("event.e" = a, "n.e" = a + b, "event.c" = c, "n.c" = c + d))
    }
    Z <- cbind(X, Y)

    Z <- data.table::as.data.table(Z)
    dat <- dplyr::summarize(dplyr::group_by(Z, study), "event.e" = sum(X == 1 & Y == 1), "n.e" = sum(X == 1), "event.c" = sum(X == 0 & Y == 1), "n.c" = sum(X == 0))
    dat <- as.data.frame(dat) # remove data.table class

    p <- meta::metabin(event.e, n.e, event.c, n.c,
                       data = dat, method = method,
                       sm = sm, comb.random = TRUE, comb.fixed = FALSE,
                       overall.hetstat = FALSE
    )$pval.random
    return(p)
  }

  # setup get.replacements
  #get.replacements <- list(function(y, x, rn, Y, X) setdiff(unique(Y[[1]]), y))
  get.replacements <- function(y, x, rn, Y, X) data.frame("outcome" = setdiff(unique(Y[[1]]), y))

  # setup only.consider through restr.study
  if (is.null(restr.study)) {
    only.consider <- NULL
  } else {
    only.consider <- rownames(X)[X$study %in% restr.study]
  }

  # finish setting up  only.consider through the choice of the sufficiently likely threshold q
  get_phatsl <- function(dat.row) c(rep(1-dat.row[1]/dat.row[2], dat.row[1]), rep(dat.row[1]/dat.row[2], dat.row[2]-dat.row[1]),
                                    rep(1-dat.row[3]/dat.row[4], dat.row[3]), rep(dat.row[3]/dat.row[4], dat.row[4]-dat.row[3]))
  phat_sl <- apply(as.matrix(data[,-1]), 1, get_phatsl)
  phat_sl <- Reduce(c, phat_sl)
  only.consider.sl <- rownames(X)[phat_sl >= q]
  #print(sum(phat_sl >= q))
  if (q>0) { # ie skip if only.consider will not be changed since q=0
    if (is.null(only.consider)) {
      only.consider <- only.consider.sl
    } else {
      only.consider <- intersect(only.consider, only.consider.sl)
    }
  }

  # get ouput
  out <- greedy.fi(X, Y, get.p.val = get.p.val, get.replacements = get.replacements, only.consider = only.consider, cl = cl, verbose = verbose)
  # what's a good way to handle our not having good patient IDs? I should put in study name at this point.. or into rownames before (in addition to having column)?
  return(out)
}


#' Calculate a fragility index for the coefficient tests in a GLM
#'
#' This is a function which outputs a GLM coefficient table with an extra column
#' which shows the fragility index corresponding to each p-value. Please note that the
#' outcome modifications in the gaussian case are determined to be symmetric around the
#' observation with radius max.step, which is distinct from the sufficiently likely approach
#' considered throughout the package.
#'
#' @param formula a formula, see the documention of 'glm' for more details
#' @param family a description of the error distribution, see the documention of 'glm' for more details
#' @param data a data frame with measurements on the columns and cases on the rows
#' @param max.step an atomic vector for the max step of each response, for when the response type
#' is restricted
#' @param alpha a number for the size of test
#' @param cl a cluster from the `parallel` package, used in `greedy.fi`
#' @param verbose a logical value for whether to print status updates while running
#'
#' @return the coefficient table from 'glm' with an extra column containing the fragility index of each coefficient
#'
#' @examples
#' data(PimaIndiansDiabetes2, package = "mlbench")
#' primdat <- PimaIndiansDiabetes2[complete.cases(PimaIndiansDiabetes2), ][1:100, ]
#' glm.fi(diabetes ~ ., binomial(), primdat, verbose = FALSE)
#'
#' @export
glm.fi <- function(formula, family, data, max.step = 1, alpha = .05, cl = NULL, verbose = FALSE) {
  if (family$family != gaussian()$family & family$family != binomial()$family) stop("Family must be binomial or gaussian.")

  if (family$family == binomial()$family) {
    get.replacements <- list(function(y, x, rn, Y, X) setdiff(unique(Y[[1]]), y)) # recall Y is left unchanged in the algorithm
  } else if (family$family == gaussian()$family) {
    get.replacements <- list(function(y, x, rn, Y, X) y + c(-1, 1) * max.step) # recall Y is left unchanged in the algorithm
  }

  out <- glm(formula, family, data)
  coef.table <- summary(out)$coefficients

  Z <- model.frame(formula, data = data)
  X <- Z[-1]
  Y <- Z[1]
  fi.coef <- c()
  for (i in 1:nrow(coef.table)) {
    get.p <- function(X, Y) summary(glm(formula, family, cbind(X, Y)))$coefficients[i, 4]

    this.fi <- greedy.fi(X = X, Y = Y, get.replacements = get.replacements, get.p.val = get.p, cl = cl, verbose = verbose, alpha = alpha)
    fi.coef <- c(fi.coef, this.fi$FI) # multiply by sign
  }
  return(cbind(coef.table, "Fragility index" = fi.coef))
}


#' Calculate a fragility index for a coefficient test in a GLM when modifying another covariate
#'
#' This function returns a fragility index (and accomponying information) for an interesting fragility index
#' which only modified a randomly observed covariate when testing the coefficient of another covariate such
#' as an intervention. This is the only example in the package which modified a covariate instead of an
#' outcome (or response). We assume that the distribution of the covariate is some Gaussian. We accomplish the
#' fragility measure by putting the covariate in the `Y` argument of `greedy.fi` and the
#' outcome in the `X` argument together with the intervention status. The function iteratively performs
#' optimization using the Brent algorithm to find the best single modification for each patient.
#'
#' @param X.regr a numeric matrix giving the covariates in the regression. The first column must
#' contain the covariate which is subject to modification.
#' @param y.regr a numeric vector giving the response in the regression, with length equal to the number of rows in X.regr
#' @param fam the family in a glm for the regression, by default binomial()
#' @param cl A parallel cluster for faster calculation in greedy.fi, by default NULL
#' @param verbose a boolean indicating whether to print status updates while running, by default TRUE
#' @param q a numeric for the minimum probability of outcome changes, by default .9
#' @param alpha a numberic for the significance threshold, by default 0.05
#'
#' @return The output of greedy.fi (a list) with an additional element which has the
#' per-patient modification likelihood bound, as in the article [Generalized fragility index].
#'
#' @examples
#' set.seed(1234567890)
#' n.age <- 200
#' beta.age <- rep(.2, 3)
#' x.age <- rnorm(n.age)
#' z.age <- rbinom(n.age, 1, 1/2)
#' eta.age <- apply(t(beta.age*t(cbind(1,x.age,z.age))),1,sum)
#' y.age <- rbinom(n.age, 1, binomial()$linkinv(eta.age))
#'
#' out <- glm.gaussian.covariate.fi(cbind(x.age, z.age), y.age, q = .7)
#'
#' @export
glm.gaussian.covariate.fi <- function(X.regr, y.regr, fam=binomial(), cl=NULL,
                                      verbose=TRUE, q=.9, alpha=0.05) {
  max.likelihood <- q # renamed

  glm.init.ggcf <- glm(y.regr~X.regr, family=fam)
  mu.start.ggcf <- glm.init.ggcf$fitted.values
  init.p.ggcf <- summary(glm.init.ggcf)$coefficients[3,4]

  get.p.val <- function(X,Y) {
    summary(glm(X[[2]]~Y[[1]]+X[[1]], # y~x+z
                family=fam,
                mustart=mu.start.ggcf))$coefficients[3,4]
  }

  mu.est <- mean(X.regr[,1]) # hard coding the first as subject to modification
  sigma.est <- sd(X.regr[,1])

  get.replacements <- function(y, x, rn, Y, X) {
    # get left and right end point
    Fy <- pnorm((y-mu.est)/sigma.est)#pt((y-mu.est)/sigma.est, nrow(Y)-1)
    if (max.likelihood <= 1-2*min(Fy, 1-Fy)) { # one end point is y
      if (y < mu.est) {
        new.outcomes <- mu.est+sigma.est*qnorm(Fy+max.likelihood)
      } else {
        new.outcomes <- mu.est+sigma.est*qnorm(Fy-max.likelihood)
      }
      endpoints <- c(y, new.outcomes)
    } else {
      endpoints <- mu.est+sigma.est*(qt(0.5+c(-1,1)*max.likelihood/2, nrow(Y)-1))
    }
    endpoints <- sort(endpoints)
    #print(endpoints)

    # set up function (with correct sign... need to minimize)
    f_greedy <- function(xx) {
      Y[rn,1] <- xx
      valtoreturn <- get.p.val(X,Y)*(-1)^(init.p.ggcf < alpha)
      return(valtoreturn)
    }

    # do L-BFGS-B
    out.optim <- optim(par=y, fn=f_greedy,
                       method='Brent',#'L-BFGS-B', gr=NULL,
                       lower=endpoints[1], upper=endpoints[2])#,
    #control=list(factr=1e5)) # i hope not extremely slow...
    return(data.frame(unique(c(y, out.optim$par))))
  }

  X <- data.frame(X.regr[,2], y.regr)
  Y <- data.frame(X.regr[,1])
  out <- greedy.fi(
    X = X, Y = Y,
    get.replacements = get.replacements, get.p.val = get.p.val,
    alpha = alpha, verbose = verbose, cl=cl
  )
  out$FI <- (-1)^(get.p.val(X,Y)>=alpha)*out$FI
  out[[length(out) + 1]] <- max.likelihood
  names(out)[length(out)] <- "max.likelihood"

  # # get sufficiently likely info
  # selected.inds <- as.numeric(rownames(out$old_responses))
  # y.mod <- y
  # y.mod[selected.inds] <- out$new_responses[,2]
  # out$sl <- normal.sl(y, y.mod)
  return(out)
}
