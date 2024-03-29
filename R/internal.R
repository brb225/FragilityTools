#' Exactly calculate a fragility index for 2x2 data (2.0)
#'
#' Do not call this function: use bin.fi with option 'exact' instead
#' The algorithm starts at a given FI value and radiates outward from there.
#' Note this can use an increasing or decreasing search, depending on the input.
#'
#' Do not call this function: use bin.fi with option 'exact' instead. If bin.fi.exact2 is
#' initialized with too large of a value, the algorithm could malfunction. It's best to
#' leave the default behavior unless sure.
#'
#' @param crosstab a 2x2 contingency table with interventions on rows and outcomes on columns
#' @param get.p a function inputting a table and outputting a p-value
#' @param alpha a scalar numeric for the significance cutoff
#' @param fi.start The starting fi considered, by default 1.
#' @param can.vary a 2x2 boolean matrix for whether each entry can decrease. See the
#' [incidence fragility] article for explanation.
#' @param start.reversed a boolean for whether the start.fi value already reverses
#' statistical significance. Note, FALSE just means unknown. (note: in a future update, this
#' will be true/false/na)
#' @param start.p The p value of the table with reversed statistical significance, paired with
#' the parameter start.reversed. Defaults to NA when knowledge of reversal is not known
#'
#' @return a list containing the signed fragility index and other accompanying values,
#' similar to `greedy.fi`
#'
#' @examples
#' x <- matrix(nrow = 2, byrow = TRUE, rgeom(4, 1 / 50))
#' colnames(x) <- c("event", "nonevent")
#' rownames(x) <- c("control", "treatment")
#' get.p <- function(tab) fisher.test(tab)$p.value
#' FragilityTools:::bin.fi.exact2(x, get.p)
bin.fi.exact2 <- function(crosstab, get.p, alpha = 0.05, fi.start = 1,
                          can.vary = matrix(rep(TRUE, 4), nrow = 2),
                          start.reversed=FALSE, start.p=NA) {
  x <- crosstab # renamed
  fi.start <- abs(fi.start) # force fi.start to not be negative

  highest.f <- sum(apply(x, 1, max))
  fi.start <- min(fi.start, highest.f) # force fi.start to not exceed number of patients

  if (is.null(rownames(x))) {
    rownames(x) <- 1:2
    #stop("Please include rownames on the data table")
  }
  # get initial p value to check if frag.index = 0
  p0 <- get.p(x)

  # get function which computes p value for given f1 and f2
  get_p_from_frow <- function(frow) {
    x2 <- x + outer(frow, c(1, -1))
    return(get.p(x2))
  }

  # get function that computes min and max p value for given fi
  # is there a smarter way to do this rather than calculuate all possible changes?
  get_extreme_p <- function(fi) {
    #print(paste0('start get extreme p: ', fi))
    f_vec <- (-fi):fi
    f <- cbind(f_vec, fi - abs(f_vec))
    f <- rbind(f, cbind(f_vec, abs(f_vec) - fi))
    colnames(f) <- c("f1", "f2")

    f[, 1] <- pmin(f[, 1], x[1, 2])
    f[, 1] <- pmax(f[, 1], -x[1, 1])
    f[, 2] <- pmin(f[, 2], x[2, 2])
    f[, 2] <- pmax(f[, 2], -x[2, 1])

    f <- subset(f, !duplicated(f)) # f[!duplicated(f),]

    if (sum(can.vary) < 4) {
      cond <- rep(TRUE, nrow(f))
      if (!can.vary[1, 1]) {
        cond <- cond & (f[, 1] >= 0)
      }
      if (!can.vary[1, 2]) {
        cond <- cond & (f[, 1] <= 0)
      }
      if (!can.vary[2, 1]) {
        cond <- cond & (f[, 2] >= 0)
      }
      if (!can.vary[2, 2]) {
        cond <- cond & (f[, 2] <= 0)
      }
      f <- subset(f, cond) # f[cond,]
    }

    #print(paste0('the dim of f is: ', nrow(f)))
    if (nrow(f) > 0) {
      ps <- apply(f, 1, get_p_from_frow)
    } else {
      ps <- Inf
    }
    #print('stop get extreme p')
    return(list(c(min(ps), max(ps)), f, c(which.min(ps), which.max(ps))))
  }

  # search
  m.ind <- ifelse(p0 < alpha, 2, 1)
  fi.current <- fi.start
  if (start.reversed) { # can just set values that work...
    p.current <- start.p
    # find a way to skip f.mat.current and whichs.current..
    f.mat.current <- NA
    whichs.current <- NA
  } else { # have to look through to understand what's going on
    extreme.calc.current <- get_extreme_p(fi.current)
    p.current <- extreme.calc.current[[1]][m.ind]
    f.mat.current <- extreme.calc.current[[2]]
    whichs.current <- extreme.calc.current[[3]]
  }

  # print(c('p.start'=p0, 'fi.start'=fi.current, 'current.p'=p.current)) ###

  all.p.vals <- c(fi.current = p.current)
  direction <- ifelse((p0 < alpha) + (p.current < alpha) != 1, 1, -1) # if started on correct side...
  # print(paste0('direction: ', direction))
  repeat {
    if (is.infinite(p.current)) { # gets around bug where fi.start=1 and p value is infinite (would give fi=0 at next step)
      fi.final <- Inf
      patients <- NA
      warning("Ran out of patients to change outcome. p val start at inf.")
      break
    }
    fi.next <- fi.current + direction
    extreme.calc.next <- get_extreme_p(fi.next)
    p.next <- extreme.calc.next[[1]][m.ind]
    f.mat.next <- extreme.calc.next[[2]]
    whichs.next <- extreme.calc.next[[3]]

    # print(c('current'=p.current, 'next'=p.next, 'fi.next'=fi.next), digits=17) ##

    all.p.vals <- c(all.p.vals, fi.next = p.next)

    if ((p.current < alpha) + (p.next < alpha) == 1) { # if crossed
      if ((p.next < alpha) + (p0 < alpha) != 1) { # if started on correct side
        fi.final <- fi.current
        if (fi.current==fi.start & start.reversed) { # if on first iteration of hacky starting reversed
          patients <- NA
        } else {
          patients <- f.mat.current[whichs.current[m.ind],]
        }
      } else {
        fi.final <- fi.next
        patients <- f.mat.next[whichs.next[m.ind],]
      }
      break
    }

    # if(is.infinite(p.next) | p.current==p.next) { # this is the problem!!
    matequal <- function(x, y) is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y)
    if (fi.next == 0 | fi.next > highest.f | (direction == 1 & matequal(f.mat.current, f.mat.next))) { # (m.ind==1&p.next>p.current) | (m.ind==2&p.next<p.current)) {
      fi.final <- Inf
      patients <- NA
      warning("Ran out of patients to change outcome.")
      break
    }

    fi.current <- fi.next
    p.current <- p.next
    f.mat.current <- f.mat.next
    whichs.current <- whichs.next
  } # end of repeat

  #print('returning from exact calc')

  rev <- (p0 >= alpha)
  list(
    "FI" = fi.final * (-1)^rev, "p_value_sequence" = unname(all.p.vals),
    "reverse" = unname(rev), "num_patients" = sum(x),
    "patients" = patients, "old_responses" = NA, "new_responses" = NA
  )
}


#' Exactly calculate a fragility index for 2x2 data (deprecated)
#'
#' Exactly calculate a fragility index using a 2x2 table.
#' Do a grid search over all FI values <= some value. This has
#' been deprecated and will stop with an error to use
#' bin.fi.exact2 instead, as it is more efficient.
#'
#' @param crosstab a 2x2 contingency table with treatments on rows and outcomes on columns
#' @param alpha a string whose value being "exact" leads to an exact calculation
#' @param get.p a function inputting a table and outputting a p-value
#' @param max.f the maximum fragility index considered
#' @param verbose a boolean value for whether to print progress
#' percentages (in increments of roughly ten percent)
#'
#' @return a list containing the signed fragility index and other accompanying values,
#' similar to `greedy.fi`.
#'
#' @examples
#' x <- matrix(nrow = 2, byrow = TRUE, rgeom(4, 1 / 50))
#' colnames(x) <- c("event", "nonevent")
#' rownames(x) <- c("control", "treatment")
#' get.p <- function(tab) fisher.test(tab)$p.value
#' FragilityTools:::bin.fi.exact(x, get.p)
bin.fi.exact <- function(crosstab, get.p, alpha = 0.05, max.f = Inf, verbose = FALSE) {
  warning("Please use bin.fi.exact2 instead") # added when depracating

  x <- crosstab # renamed
  if (is.null(rownames(x))) stop("Please include rownames on the data table")
  # get initial p value to check if frag.index = 0
  pp <- get.p(x)

  # set up
  # get p value after modifying with f values
  get.p.mod <- function(f1, f2, x) {
    mod.x <- matrix(nrow = 2, byrow = TRUE, c(
      x[1, 1] + f1, x[1, 2] - f1,
      x[2, 1] + f2, x[2, 2] - f2
    ))
    get.p(mod.x)
  }

  # grid of f's to consider
  f1.vals <- seq(-min(x[1, 1], max.f), min(x[1, 2], max.f))
  f2.vals <- seq(-min(x[2, 1], max.f), min(x[2, 2], max.f))

  counter <- 0
  tot <- length(f1.vals) * length(f2.vals)

  # calculate grid of p values from each (f1, f2) pair
  f.mat <- matrix(nrow = length(f1.vals), ncol = length(f2.vals))
  p.mat <- matrix(nrow = length(f1.vals), ncol = length(f2.vals))
  for (i in 1:length(f1.vals)) {
    for (j in 1:length(f2.vals)) {
      counter <- counter + 1
      # if (counter %% round(tot/10) == 0 & verbose) print (counter/tot)

      f1 <- f1.vals[i]
      f2 <- f2.vals[j]

      f.mat[i, j] <- abs(f1) + abs(f2)
      p.mat[i, j] <- get.p.mod(f1, f2, x)
    }
  }
  rownames(p.mat) <- f1.vals
  colnames(p.mat) <- f2.vals

  # extract best (f1, f2) pair
  ## rearrange data to convenient format
  d <- data.frame(
    f1.val = rep(row.names(p.mat), ncol(p.mat)),
    f2.val = rep(colnames(p.mat), each = nrow(p.mat)),
    f.index = as.vector(f.mat),
    p.value = as.vector(p.mat),
    stringsAsFactors = FALSE
  )

  ## remove outside of max.f consideration
  d <- d[d$f.index <= max.f, ]

  ## select (f1,f2)'s which give closest p val which is on other side
  rel.rows <- if (pp < alpha) which(d$p.value >= alpha) else which(d$p.value < alpha)
  if (length(rel.rows) != 0) { # if there are (f1,f2) vals which give the opposite significance
    d <- d[rel.rows, ] # only look at (f1,f2) which give opposite
    d <- d[min(d$f.index) == d$f.index, ] # only look at smallest f index
  }
  if (pp < alpha) { # picks most extreme p val if more than one with smallest f index
    d <- d[which.max(d$p.value), ]
  } else {
    d <- d[which.min(d$p.value), ]
  }

  row.frag.index <- cbind(d, "reverse" = (pp >= alpha))
  colnames(row.frag.index)[1:2] <- rownames(x)

  rev <- (pp >= alpha)
  # get in same format as greedy.fi
  return(list(
    "FI" = row.frag.index[1, 3] * (-1)^rev, "p_value_sequence" = c(pp, row.frag.index[1, 4]),
    "reverse" = rev, "num_patients" = sum(x),
    "patients" = row.frag.index[, 1:2],
    "old_responses" = NA, "new_responses" = NA
  ))
}


#' Original fragility index calculation, for binary data
#'
#' This function implements the original algorithm to calculate a fragility index,
#' as it was proposed by Walsh et al. (2014). The algorithm determines a single group
#' in which to make outcome modifications. We also implement a variant which considers
#' outcome modifications in both directions (i.e. increasing and decreasing event counts),
#' which is a hybrid of the algorithm due to Walsh et al. (2014) and the algorithm
#' commonly used for the "reverse" fragility index.
#'
#' @param crosstab a 2x2 contingency table with groups on the rows
#' @param get.p a function which accepts a 2x2 matrix and outputs a p value
#' @param alpha a numeric for the significance cutoff, 0.05 by default
#' @param dir a character, either "left", "right", or "both". The default is "left". Walsh originally used
#' 'left' which increases the outcome count in the left column (which is typically events)
#' @param group a character specifying how to choose the single group in which to make modifications.
#' The options are: 'event' for the fewest events (default), 'nonevent' for the fewest nonevents, or
#' 'both' for the fewest overall. We assume that events are in the first column.
#'
#' @return a list containing the signed fragility index and other accompanying values,
#' similar to `greedy.fi`.
#'
#' @examples
#' FragilityTools:::bin.fi.walsh(matrix(c(100, 96, 13, 28), nrow = 2), function(mat) fisher.test(mat)$p.value)
bin.fi.walsh <- function(crosstab, get.p, alpha = 0.05, dir='both', group='event') {
  mat <- crosstab # renamed

  start.p <- get.p(mat)
  if (group=='event') {
    group.sizes <- mat[,1] # event counts in each group
  } else if (group=='nonevent') {
    group.sizes <- mat[,2] # nonevent counts in each group
  } else if (group=='both') {
    group.sizes <- apply(mat, 1, min) # smallest entry in each group
  }
  smallest.grs <- unname(which(group.sizes == min(group.sizes)))[1]

  get.mod.p <- function(f, smallest.gr) {
    mod.mat <- mat
    mod.mat[smallest.gr, ] <- mod.mat[smallest.gr, ] + c(1, -1) * f
    #print(mod.mat)
    return(get.p(mod.mat))
  }
  do.grid.search <- function(smallest.gr, dir) {
    end.p <- NA
    for (f in 1:max(mat[smallest.gr, ])) {
      if ((dir != 'right') && -mat[smallest.gr, 1] <= f && f <= mat[smallest.gr, 2]) { # check if rightward search is in bounds
        p1 <- get.mod.p(f, smallest.gr)
        if ((p1 - alpha) * (start.p - alpha) < 0) {
          end.p <- p1
          end.f <- f
          break
        }
      }
      if ((dir != 'left') && -mat[smallest.gr, 1] <= -f && -f <= mat[smallest.gr, 2]) { # check if leftward search is in bounds
        p2 <- get.mod.p(-f, smallest.gr)
        if ((p2 - alpha) * (start.p - alpha) < 0) {
          end.p <- p2
          end.f <- -f
          break
        }
      }
    }
    if (is.na(end.p)) {
      end.f <- ifelse(start.p<alpha, Inf, -Inf)
    }
    return(list(f = end.f, end.p = end.p))
  }

  pats <- c(0, 0)
  out <- do.grid.search(smallest.grs, dir)
  pats[smallest.grs] <- out$f

  # if (length(smallest.grs) == 1) {
  #   out <- do.grid.search(smallest.grs)
  #   pats[smallest.grs] <- out$f # does this correctly handle the sign of f?
  # } else {
  #   out1 <- do.grid.search(smallest.grs[1])
  #   out2 <- do.grid.search(smallest.grs[2])
  #   if (abs(out1$f) <= abs(out2$f)) {
  #     out <- out1
  #     pats[1] <- out1$f # does this correctly handle the sign of f?
  #   } else {
  #     out <- out2
  #     pats[2] <- out2$f # does this correctly handle the sign of f?
  #   }
  # }
  f <- out$f
  end.p <- out$end.p

  return(list(
    "FI" = abs(f) * (-1)^(start.p >= alpha), "p_value_sequence" = c(start.p, end.p),
    "reverse" = (start.p >= alpha), "num_patients" = sum(mat),
    "patients" = pats,
    "old_responses" = NA, "new_responses" = NA
  ))
}


#' Greedy calculation of fragility indices for 2x2 contingency tables
#'
#' This is an internal function, please use bin.fi instead
#'
#' @param crosstab a 2x2 contingency table with interventions on rows and outcomes on columns
#' @param get.p a function inputting a table and outputting a p-value
#' @param alpha a scalar numeric for the significance cutoff, by default 0.05
#' @param can.vary a 2x2 boolean matrix for whether each entry can decrease. See the
#' [incidence fragility index] article for explanation.
#' @param crosstab.offset a 2x2 contingency table with interventions on rows and outcomes on columns.
#' The default has all entries equal to 0. These patients are not subject to modification, unlike crosstab,
#' but they do contribute to the evaluation of the p value, like crosstab. Note, this is a variant of the
#' 'dont.consider' argument in greedy.fi
#'
#' @return a list containing the signed fragility index, in addition to the modified contingency
#' table which has reversed statistical significance and the p value sequence (as in `greedy.fi`)
#'
#' @examples
#' x <- matrix(nrow = 2, byrow = TRUE, rgeom(4, 1 / 50))
#' colnames(x) <- c("event", "nonevent")
#' rownames(x) <- c("control", "treatment")
#' get.p <- function(tab) fisher.test(tab)$p.value
#' FragilityTools:::bin.fi.greedy(x, get.p)
bin.fi.greedy <- function(crosstab, get.p, alpha=.05,
                          can.vary = matrix(rep(TRUE, 4), nrow = 2),
                          crosstab.offset = matrix(rep(0,4), nrow = 2)) {
  x <- crosstab # renamed

  get.p.val <- get.p
  get.p <- function(mm) get.p.val(mm+crosstab.offset) # redefine get.p to always augment data with fixed patients

  init.p <- get.p(x)
  signFI <- (-1)^(init.p>=alpha)
  FIcount <- 0

  p.seq <- init.p

  mats <- list(outer(c(1,0),c(-1,1)),
               outer(c(0,1),c(-1,1)),
               outer(c(1,0),-c(-1,1)),
               outer(c(0,1),-c(-1,1)))
  mats <- mats[c(can.vary)] # order of mats picked so this works

  get.p.mod <- function(mod) {
    mat <- x+mod
    if (any(mat<0)) return(NA)
    return(get.p(mat))
  }

  while (TRUE) {
    FIcount <- FIcount+1
    out <- plyr::laply(mats, get.p.mod)

    if (all(is.na(out))) { # no feasible modifications
      FIcount <- Inf
      break
    }
    if (FIcount > sum(x)) { # ran out of patients to modify, will also prevent loops due to p being constant
      FIcount <- Inf
      break
    }
    mout <- ifelse(init.p<alpha, max(out, na.rm=TRUE), min(out, na.rm=TRUE))
    if (FIcount > 1 & (alpha - init.p)*(p.seq[FIcount] - mout) > 0) { # cannot make progress
      FIcount <- Inf
      break
    }

    if (init.p < alpha) {
      ind <- which.max(out)
    } else {
      ind <- which.min(out)
    }
    x <- x + mats[[ind]]

    p.seq <- c(p.seq, get.p(x))
    if ((init.p < alpha & get.p(x) >= alpha) | (init.p >= alpha & get.p(x) < alpha)) break
  }

  return(list(FI=FIcount*signFI, mat=x, p_value_sequence = p.seq, rev = init.p>=alpha))
}




#' Stochastic generalized fragility indices for 2x2 tables
#'
#' Use bin.fi instead of this function.
#' This is an internal function which calculates stochastic generalized fragility indices for 2x2 tables
#' in a more optimized way than stochastic.fi.
#'
#' @param crosstab a 2x2 contingency table with interventions on rows and outcomes on columns
#' @param can.vary a 2x2 boolean matrix for whether each entry can decrease. See the
#' [incidence fragility index] article for explanation.
#' @param get.p a function that inputs a matrix and returns a p value
#' @param r the index of the stochastic fragility index, by default 0.5. Having r=0 is equivalent to the generalized fragility
#' index and having r=1 means that all patient combinations of the output size can reverse significance.
#' @param nsim The number of simulations in the root finding algorithm, by default 10
#' @param gfi.init An initialization of the output size, by default 10
#' @param alpha a numeric for the significance cutoff
#' @param verbose a logical value for whether to print status updates while running
#' @param cl a cluster from the `parallel` package, used to compute fragility index over
#' each modified observation at each stage of the greedy algorithm
#' @param D a parameter of Polyak-Ruppert averaging, by default 40
#' @param gamma a parameter of Polyak-Ruppert averaging, by default .2
#'
#' @return a length 2 list, with the first entry giving the stochastic generalized fragility index and the
#' last entry giving the history of the root finding algorithm.
#'
#' @examples
#' x <- matrix(nrow = 2, byrow = TRUE, rgeom(4, 1 / 50))
#' colnames(x) <- c("event", "nonevent")
#' rownames(x) <- c("control", "treatment")
#' get.p <- function(tab) fisher.test(tab)$p.value
#' FragilityTools:::bin.fi.stochastic(x, get.p)
#'
#' @export
bin.fi.stochastic <- function(crosstab, get.p,
                              can.vary = matrix(rep(TRUE, 4), nrow = 2),
                              r=0.5, nsim=10, qfi.init = 10L,
                              alpha=.05, verbose=FALSE, cl = NULL,
                              D=40, gamma=.2) {
  # note, this algorithm only gives some root, not the smallest root.. for eg it will likely "fail" for r=1

  #do some checks
  if (r==1) warning('The output cannot be trusted for this r: the output may be larger than the minimum such SGFI')

  out.fi <- bin.fi.greedy(crosstab, get.p, alpha, can.vary)
  classical.fi <- abs(out.fi$FI)
  if (is.infinite(classical.fi)) return(Inf) # stop if reversing is infeasible with all data
  # any way to modify find_zero to avoid this check?

  # the noisy function to find root of expectation
  pat.long <- c(rep('a', crosstab[1,1]), rep('b', crosstab[1,2]), rep('c', crosstab[2,1]), rep('d', crosstab[2,2]))
  centered_prob_of_rev <- function(ss) { # should parallelize over this instead of greedy.fi
    ss <- floor(ss)
    did_reverse <- c()
    for (sim in 1:nsim) {
      # get patients who have permitted modifications
      x.v <- sample(pat.long, ss)
      x.v <- matrix(nrow=2,byrow=TRUE, c(sum(x.v=='a'), sum(x.v=='b'), sum(x.v=='c'), sum(x.v=='d')))

      suppressWarnings(
        out <- bin.fi.greedy(x.v, get.p, alpha, can.vary, crosstab-x.v)
      )
      did_reverse <- c(did_reverse, is.finite(out$FI))
    }
    return(mean(did_reverse) - r)
  }

  fz <- find_zero(centered_prob_of_rev, x.init = qfi.init,
                  D=D, burnin_dur = 10, gamma=gamma, eps=.05, fz.verbose=verbose,
                  limits=c(classical.fi, sum(crosstab)), proj = function(a) max(min(a, sum(crosstab)), classical.fi))
  fz$x <- ceiling(fz$x)*(-1)^out.fi$rev
  names(fz)[1] <- 'FI'
  return(fz)
}








#' Find when an increasing function observed with noise crosses x-axis
#'
#' Find the smallest positive integer x-value larger than in initial
#' value for which an approximately monotonic function f is nonnegative.
#' This is an internal function used in `general.fi.samplesize`. The
#' argument is increased exponentially until it becomes nonnegative, then
#' a binary search is performed.
#'
#' @param f a function with one argument
#' @param x.init the initial value towards finding the zero of f
#' @param fz.verbose a logical value for whether to print status updates while running
#' @param D an internal argument controlling the step size in Ruppert-Polyak averaging
#' @param burnin_dur the number of iterations of SGD to burn in before there's checks for
#' the average to converge
#' @param eps a parameter to control the error. The smaller is it, the more precise the
#' output but the longer the function will take to run.
#' @param proj a function to project/map the argument of the function onto the functions domain
#' @param limits a length two numeric vector with the first entry being the lower bound and the
#' second being an upper bound.
#' @param gamma the power of n^{-1} in gradient descent
#'
#' @return the argument of f which satisfies the desired conditions.
#'
#' @examples
#' FragilityTools:::find_zero(function(x) x - 100 + rnorm(1), fz.verbose = FALSE)
find_zero <- function(f, x.init = 10L, fz.verbose = FALSE, D = 1, gamma = .6,
                      burnin_dur = 16, eps = .1, proj = function(a) a, limits = c(1, 9999999),
                      init.step=TRUE) {
  # updates proj to take into account the upper and lower limits
  proj. <- function(a) {
    out <- proj(a)
    out <- min(out, limits[2])
    out <- max(out, limits[1])
    return(out)
  }

  # tries to get lowest value with f>=0
  x <- x.init
  history <- matrix(c("x" = 0, "y" = 0, "avg" = 0),nrow=1)[0,]
  colnames(history) <- c('x', 'y', 'avg')

  if (init.step) {
    y <- f(x)
    history <- rbind(history, c("x" = x, "y" = y, "avg" = Inf))
    if (fz.verbose) print(c(x, y))

    ## init
    if (y < 0) {
      # increase exponentially until changing sign of f
      while (y < 0) {
        old_y <- y
        old_x <- x

        x <- proj.(2 * old_x)
        y <- f(x)

        if (fz.verbose) print(c(x, y))
        history <- rbind(history, c("x" = x, "y" = y, "avg" = Inf))

        if (old_x == limits[2] & x == limits[2]) break
      }
      if (fz.verbose) print("finished increasing")

      x <- (y * old_x - old_y * x) / (y - old_y)
      x <- proj.(x)
    } else if (y == 0) {
      x <- x
    } else {
      # decrease exponentially until changing sign of f
      while (y > 0) {
        old_y <- y
        old_x <- x

        x <- proj.(old_x / 2)
        y <- f(x)

        if (fz.verbose) print(c(x, y))
        history <- rbind(history, c("x" = x, "y" = y, "avg" = Inf))

        if (old_x == limits[1] & x == limits[1]) break
      }
      if (fz.verbose) print("finished decreasing")

      if (!(x %in% limits)) {
        x <- (y * old_x - old_y * x) / (y - old_y)
        x <- proj.(x)
      }
    }
  }

  # do ruppert-polyak averaging
  ## burn in
  y <- c(f(x[1]))
  if (fz.verbose) print(c(x[1], y[1]))
  history <- rbind(history, c("x" = x, "y" = y, "avg" = Inf))

  old_avg <- x[1]
  if (fz.verbose) print("start ruppert polyak averaging burn in")
  for (t in 2:burnin_dur) {
    x[t] <- proj.(x[t - 1] - (D * (t - 1)^(-gamma)) * y[t - 1])
    y[t] <- f(x[t])

    if (fz.verbose) print(c(x[t], y[t]))
    history <- rbind(history, c("x" = x[t], "y" = y[t], "avg" = Inf))
  }
  old_x <- x[burnin_dur]
  old_y <- y[burnin_dur]

  ## iterate until convergence
  if (fz.verbose) print("start iterating until convergence")

  x <- c(old_x)
  y <- c(old_y)

  if (fz.verbose) print(c(x[1], y[1]))
  history <- rbind(history, c("x" = x, "y" = y, "avg" = Inf))

  t <- 2 # burnin_dur %/% 2
  old_avg <- x[1]
  while (TRUE) {
    x[t] <- proj.(x[t - 1] - (D * (t - 1)^(-gamma)) * y[t - 1])
    y[t] <- f(x[t])

    if (fz.verbose) print(c(x[t], y[t], mean(x)))
    history <- rbind(history, c("x" = x[t], "y" = y[t], "avg" = mean(x)))

    if (abs(old_avg - mean(x)) < eps) break
    old_avg <- mean(x)
    t <- t + 1
  }
  if (fz.verbose) print("finished")

  final_x <- proj.(mean(x))
  return(list("x" = final_x, "history" = history))
}


#' Convert a 2x2 table into X,Y data frame format
#'
#' This is a helper function since FragilityTools uses X,Y data frames by default,
#' but uses 2x2 tables when possible for improved computational efficiency.
#'
#' @param mat 2x2 numeric matrix
#' @param event.col either 1 or 2, specifying which column has events, by default 1
#'
#' @return list containing two one-column data frames X and Y
#'
#' @examples
#' FragilityTools:::mat_to_xy(matrix(1:4, nrow = 2))
mat_to_xy <- function(mat, event.col = 1) {
  # get names
  r1 <- rownames(mat)[1]
  r2 <- rownames(mat)[2]
  c1 <- colnames(mat)[event.col]
  c2 <- colnames(mat)[3 - event.col]

  if (is.null(r1)) r1 <- 'group1'
  if (is.null(r2)) r2 <- 'group2'
  if (is.null(c1)) c1 <- 'event'
  if (is.null(c2)) c2 <- 'nonevent'

  # get relevant info
  xx <- mat[1, event.col] # "treatment"
  yy <- mat[2, event.col] # "control"

  m <- sum(mat[1, ])
  n <- sum(mat[2, ])

  # build X's
  x1 <- rep(r1, m)
  x2 <- rep(r2, n)
  X <- c(x1, x2)

  # build Y's
  y11 <- rep(c1, xx)
  y12 <- rep(c2, m - xx)
  y21 <- rep(c1, yy)
  y22 <- rep(c2, n - yy)
  Y <- c(y11, y12, y21, y22)

  # make dataframe with factors
  X <- data.frame(factor(X, levels = c(r1, r2))) # c(0,1)))
  Y <- data.frame(factor(Y, levels = c(c1,c2))) # c(0,1)))

  rownames(X) <- 1:nrow(X)
  colnames(X) <- "X"
  rownames(Y) <- 1:nrow(Y)
  colnames(Y) <- "Y"

  # return
  return(list("X" = X, "Y" = Y))
}


#' Finds beta inverse dispersion parameter when given a center and a multiplier to determine the scale
#'
#' This is an internal function for use in the LTFU-aware
#' fragility index calculation functions.
#'
#' Gets parameters for sampling in ltfu.fi
#' Chooses them so that the mode is po and the 75th
#' quantile is mult*po
#'
#' @param po the location parameter of the beta distribution
#' @param mult the multiplier which is used to determine the scale. The 75th
#' quantile of the beta distribution will be mult*po.
#'
#' @return a length one numeric with the beta inverse disperion parameter
#'
#' @examples
#' params <- FragilityTools:::get_beta_parameters(.3)
get_beta_parameters <- function(po, mult=1.3) {
  alpha <- function(s) s*po+1
  beta <- function(s) s-s*po+1

  get75q <- function(s) qbeta(.875, alpha(s), beta(s))
  shat <- uniroot(function(s) get75q(s)-mult*po, c(0,999999))$root
  return(shat)
}


#' Simulate right-censored survival responses in two groups
#'
#' @param n the number of observations
#' @param num.times the number of time points at which events or right-censors can occur
#' @param p1 the proportion of events in the first group, by default NULL which takes p1 to be random
#' @param p2 the proportion of events in the second group, by default NULL which takes p1 to be random
#'
#' @return a dataframe with three columns: group, time, and status
#'
#' @examples
#' dat <- FragilityTools:::get.survival.data(n = 100, num.times = 6)
get.survival.data <- function(n, num.times, p1 = NULL, p2 = NULL) {
  # roughly 50/50 in control and treatment
  dat <- data.frame(person_id = 1:n, group = sample(c("control", "treatment"), n, TRUE), stringsAsFactors = FALSE)

  # generate event times depending on group
  all.times <- sort(sample(1:(5 * num.times), num.times)) # pick 10 possible event times times
  # all.times <- all.times - min(all.times) #makes 0 the first time, for plotting purposes
  if (is.null(p1) & is.null(p2)) {
    p <- list("control" = runif(num.times, 0, min(1, 1 / (num.times / 4))), "treatment" = runif(num.times, 0, min(1, 1 / (num.times / 2))))
  } else {
    p <- list("control" = rep(p1, num.times), "treatment" = rep(p2, num.times))
  }
  event.time <- function(group) {
    p <- p[[group]]

    time.of.death <- Inf # for survival censoring at end of study
    for (i in 1:num.times) {
      if (runif(1) < p[i]) {
        time.of.death <- all.times[i]
        break
      }
    }

    censored <- sample(c("yes", "no"), 1, TRUE, c(.1, .9))
    if (time.of.death == Inf) {
      time.of.death <- max(all.times) + 1 # edits to remove Inf
      censored <- "yes"
    }
    list(censored, time.of.death)
  }

  info <- sapply(dat$group, event.time)
  data.frame(
    "group" = dat[, 2],
    "time" = unlist(info[2, ]),
    "status" = (unlist(info[1, ]) == "no") + 0,
    stringsAsFactors = FALSE
  )
}


#' Calculate level of smallest HDR which contains the modified observation
#'
#' Used for internal functions in the study of generalized fragility indices.
#'
#' @param y a numeric vector of outcomes
#' @param y.mod a numeric vector of outcomes after the fragility index modification
#'
#' @return a scalar numeric: a ratio of the likilihoods in y and y.mod, estimated
#' using the MLE with the data y
#'
#' @examples
#' y <- rnorm(100, 0, 1)
#' y.mod <- rnorm(100,.1, 1.1) # should be returned from an FI procedure
#' FragilityTools:::normal.sl(y, y.mod)
normal.sl <- function(y, y.mod) {
  n <- length(y)
  sdy <- sd(y)#*sqrt(n/(n-1))

  logp_real <- sum(dnorm(y, mean(y), sdy, log=TRUE))
  logp_mod <- sum(dnorm(y.mod, mean(y), sdy, log=TRUE))
  return(exp(logp_real-logp_mod)) # this can return > 1 sometimes... should we only consider sufficient statistics?
  #return(abs(exp(logp_real)-exp(logp_mod))) # this returns huge numbers


  #
  #   # get test statistics
  #   T.mod <- (mean(y.mod)-mu0)/sd(y.mod)*sqrt(n)
  #   T.real <- (mean(y)-mu0)/sd(y)*sqrt(n)
  #
  #   # pick out specific non central t distribution
  #   df <- n-1
  #   ncp <- sqrt(n)/sd(y)*(mean(y)-mu0)
  #
  #   # short circuit and return cdf difference
  #   return(abs(pt(T.mod, df, ncp)-pt(T.real, df, ncp)))
  #
  #
  #
  #
  #
  #
  #   nct.exp.val <- ncp*sqrt(df/2)*gamma((df-1)/2)/gamma(df/2)
  #   lb <- nct.exp.val*sqrt(df/(df+5/2))
  #   ub <- nct.exp.val*sqrt(df/(df+1)) # bound from wikipedia
  #   T.high <- suppressWarnings( # find mode
  #     optimize(function(x) dt(x, df, ncp, log=TRUE),
  #              interval=c(-lb, ub),
  #              maximum=TRUE)$maximum
  #     )
  #
  #   # if T.mod has higher density (based on pre-hdi area)
  #   if (dt(T.mod, df, ncp)>dt(T.real, df, ncp)) return(abs(pt(T.mod, df, ncp)-pt(T.real, df, ncp)))
  #
  #   # if greater than mode, search
  #   if (T.mod > T.high) {
  #     rt <- uniroot(function(TT) dt(TT, df, ncp)-dt(T.mod, df, ncp),
  #                    lower=T.high-2*(T.mod-T.high),
  #                    upper=T.high)$root
  #   }
  #   if (T.mod < T.high) {
  #     rt <- uniroot(function(TT) dt(TT, df, ncp)-dt(T.mod, df, ncp),
  #                    lower=T.high,
  #                    upper=T.high+2*(T.high-T.mod))$root
  #   }
  #   return(abs(pt(rt, df, ncp)-pt(T.mod, df, ncp)))
  #
  #   stop('error!! got through all ifs without returning..')
}


#' Calculate level of smallest HDR which contains the modified observation
#'
#' Used for internal functions in the study of generalized fragility indices.
#'
#' @param y a numeric vector of outcomes
#' @param y.mod a numeric vector of outcomes after the fragility index modification
#' @param hatbeta a numeric vector of estimated coefficients
#' @param Xregr a numeric matrix of covariates used to explain the outcomes
#'
#' @return a scalar numeric: a ratio of the conditional likilihoods in y and y.mod, estimated
#' using the MLE with the data y
#'
#' @examples
#' y <- rpois(100, 1)
#' y.mod <- rpois(100, 1.2) # should be returned from an FI procedure
#' FragilityTools:::poisson.regr.sl(y, y.mod, 0, matrix(rnorm(100), ncol=100))
poisson.regr.sl <- function(y, ymod, hatbeta, Xregr) { # we have to do glm because knowing the median (i.e. qr) doesn't tell u about the whole distribution
  eta.pois <- apply(t(hatbeta*t(cbind(1,Xregr))),1,sum)
  empirical.mean <- poisson()$linkinv(eta.pois)

  logp_real <- sum(dpois(y, empirical.mean, log=TRUE))
  logp_mod <- sum(dpois(y.mod, empirical.mean, log=TRUE))
  return(exp(logp_real-logp_mod))
}
