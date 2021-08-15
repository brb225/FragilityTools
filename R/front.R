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

  # organize
  mat <- matrix(c(xx, yy, m - xx, n - yy), nrow = 2)
  rownames(mat) <- c('group1', 'group2')
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

#' Sample size calculator for 2x2 tables, taking into account power and fragility index
#'
#' This function is a wrapper around general.fi.samplesize, see the documentation there.
#'
#' @param min.fi the smallest acceptable QUANTILE fragility index. When NULL, the FI calculation is skipped
#' and sample_size_init_fi is taken to produce the desired FI.
#' @param min.power the smallest acceptable power. When NULL, the power calculation is skipped and
#' sample_size_init_power is taken to produce the desired power.
#' @param alpha a numeric for the size of the p value based test
#' @param tau the quantile of FI to bound, default 1/2
#' @param algorithm A string specifying the algorithm to use to calculate fragility indices.
#' The default is "walsh" (or 'original'). Alternatives include "greedy" for using greedy.fi
#' @param test a string specifying which p value based test to use. By default 'fisher' for Fisher's exact test.
#' Can also specify `chisq.prop` for Pearson's chi square test.
#' @param theta1 a numeric for the event rate in the first group, for simulating the alternative model
#' @param theta2 a numeric for the event rate in the second group, for simulating the alternative model
#' @param row.prop a numeric for the proportion of patients in the first group, for simulating the alternative model
#' @param niters the number of iterations to run the algorithm. The final output is the median across iterations.
bin.fi.samplesize <- function(min.fi = 10, min.power = .8, alpha = .05, tau=.5,
                              test = 'fisher', algorithm='original', theta1=.3, theta2=.7, row.prop=1/2,
                              niters=5, cl=NULL) {
  if (test=='fisher') {
    get.p.val <- function(m) fisher.test(m)$p.value
  } else if (test == "chisq.prop") {
    get.p.val <- function(tab) {
      n1 <- sum(tab[1, ])
      n2 <- sum(tab[2, ])
      p1 <- tab[1, 1] / n1
      p2 <- tab[2, 1] / n2

      pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
      ts <- (p1 - p2) / sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
      p_value <- 2 * pnorm(abs(ts), lower.tail = FALSE)
      return(ifelse(is.nan(p_value), 1, p_value))
    }
  } else {
    stop('Please choose an appropriate test option')
  }

  if (algorithm=='original') algorithm <- 'walsh'

  if (algorithm=='walsh') {
    get.sample <- function(ss) draw.binom(ss, theta1=theta1, theta2=theta2, row.prop = row.prop, matrix=TRUE)
  } else if (algorithm=='greedy') {
    get.sample <- function(ss) draw.binom(ss, theta1=theta1, theta2=theta2, row.prop = row.prop, matrix=FALSE)
  }

  general.fi.samplesize(min.fi = min.fi, min.power = min.power, alpha = alpha, tau=tau,
                        get.p.val = get.p.val, get.sample = get.sample,
                        algorithm=algorithm, niters=niters, cl=cl)
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
#' @param get.p.val a fucntion that inputs a matrix and returns a p value, otherwise a function
#' that inputs X and Y and returns a p value when alg='greedy'.
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
#' The default is "walsh". Alternatives include "greedy" for using greedy.fi
#' @param get.replacements a function which outputs a data frame containing all possible row replacements
#' for Y which are to be considered. The functions inputs the row of Y under consideration,
#' the row of X under consideration, the row name, and the full original data frames X and Y.
#' This gets used in the greedy algorithm.
#'
#' @return the length two numeric vector of the calculated sample sizes for the desired power and fragility index
#'
#' @examples
#' ss <- general.fi.samplesize(min.fi = 10, min.power = .8, alpha = .05, tau=.5,
#' get.p.val = function(m)fisher.test(m)$p.value, get.sample = function(ss) draw.binom(ss, matrix=TRUE),
#' nsim = 10, niters = 2, verbose = TRUE)
#'
#' ss <- general.fi.samplesize(min.fi = 10, min.power = .8, alpha = .05, tau=.5,
#' get.p.val = function(X,Y) fisher.test(table(X[,1], Y[,1]))$p.value,
#' get.sample = function(ss) draw.binom(ss, matrix=FALSE),
#' alg='greedy', niters=1, verbose = TRUE,
#' get.replacements=function(y, x, rn, Y, X) data.frame(setdiff(c("event", "nonevent"), y)))
#'
#' @export
general.fi.samplesize <- function(min.fi = 10, min.power = .8,
                                  sample_size_init_power = 100L, sample_size_init_fi = NULL,
                                  get.p.val, get.replacements, get.sample, gamma = 0.6, niters = 50,
                                  cl = NULL, verbose = FALSE, alpha = .05, tau = 1 / 2, nsim = 30,
                                  eps = .1, algorithm = "walsh") {

  if (algorithm=='original') algorithm <- 'walsh'

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
          fi_vals_i <- bin.fi.walsh(mat, get.p.val, alpha, dir='left', group='event')$FI
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

    #print(ggplot2::qplot(ss.s))
    sample_size1 <- ceiling(quantile(ss.s, .5, names = FALSE)) # the points of the niters is to loop over the whole thing to stabilize..
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

    #print(ggplot2::qplot(ss.s))
    sample_size2 <- ceiling(quantile(ss.s, .5, names = FALSE))
  } else {
    if (verbose) print("skipped fi calculations")
    sample_size2 <- sample_size_init_fi
  }

  # take largest sample size
  sample_size <- c("p_ss" = sample_size1, "fi_ss" = sample_size2)
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

    last_sample_size_power <- unname(sample_size["p_ss"])
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
#' @param effect.size.plt A numeric, by default FALSE. It determines whether the plot will color tiles by the effect size
#' (when TRUE) or by statistical significance (when FALSE) of the augmented data including both the observed and missing
#' patients. Assumes events are in the first column.
#' @param extreme a numeric, either 0, 1, or 2. 1 leads to a standard prior specification. 0 shrinks towards 0, 2 shrinks
#' towards 1/2. Experimental--please do not change.
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
                    fig.size=1.1, gradientn.scale=.99, effect.size.plt = FALSE, extreme=1) {
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
      print(paste0('Group ',j,': s=',s))
      if (extreme==1) {
        pl <- rbeta(ndraws, s*po+1, s*(1-po)+1)
      } else if (extreme==0) {
        pl <- rbeta(ndraws, s*po/2+1, s*(1-po/2)+1)
      } else if (extreme==2) {
        pl <- rbeta(ndraws, s*(po/2+1/4)+1, s*(1-(po/2+1/4))+1)
      } else {
        stop("You set a wrong value of 'extreme', but you weren't supposed to edit the value in the first place!")
      }
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

  aug_effect_size <- function(rr) {
    g1 <- rr[1]; g2 <- rr[2];

    return('es'=(crosstab[1,1]+g1)/(sum(crosstab[1,])+ltfu[1]) - (crosstab[2,1]+g2)/(sum(crosstab[2,])+ltfu[2]))
  }
  if (effect.size.plt) { # get augmented effect sizes
    df <- cbind(df, 'es'=apply(df, 1, aug_effect_size))
  }

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

  if (effect.size.plt) {
    plt.dat$es <- as.numeric(plt.dat$es)
    plt <- ggplot(data=plt.dat[order(plt.dat$mode),],
                  aes(x=g1/ltfu[1], y=g2/ltfu[2], fill=es))+
      geom_tile(aes(color=log_prs), size=fig.size)+
      scale_color_gradientn(name='Posterior \nprobability',
                            colours=c("white","blue","black"),
                            breaks=c(maxpr/2, maxpr))+#n.breaks=3)+
      scale_fill_continuous(name='Effect size', breaks=trunc(c(min(plt.dat$es), max(plt.dat$es))*1000)/1000,
                            low = "grey5", high = "orange")+
      labs(x=xlab, y=ylab)+
      theme_bw()+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(legend.position='bottom')
  } else {
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
  }

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




























