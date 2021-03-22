#' Calculate an exact fragility index for the fisher exact test (2.0)
#'
#' Exactly alcululate a fragility index using for the fisher exact test using a table.
#' Start at a given FI value and radiate out from there.
#'
#' Do not call this function: use bin.fi with option 'exact' instead. If bin.fi.exact2 is
#' initialized with too large of a value, the algorithm could malfunction. It's best to
#' leave the default behavior unless sure.
#'
#' @param x a 2x2 matrix with treatments on rows and outcomes on columns
#' @param alpha a string whose value being "exact" leads to an exact calculation
#' @param get.p a function inputting a table and outputting a p-value
#' @param fi.start The starting fi considered, by default 1.
#' @param can.vary a 2x2 boolean matrix for whether each entry can change.
#' @param start.reversed a boolean for whether the start.fi value already reverses
#' statistical significance. Note, FALSE just means unknown.
#' @param start.p The p value of the table with reversed statistical singificance, paired with
#' the parameter start.reversed
#' @return a fragility index
#'
#' @examples
#' x <- matrix(nrow = 2, byrow = TRUE, rgeom(4, 1 / 50))
#' colnames(x) <- c("event", "nonevent")
#' rownames(x) <- c("control", "treatment")
#' get.p <- function(tab) fisher.test(tab)$p.value
#' bin.fi.exact2(x, get.p)
bin.fi.exact2 <- function(x, get.p, alpha = 0.05, fi.start = 1, can.vary = matrix(rep(TRUE, 4), nrow = 2),
                          start.reversed=FALSE, start.p=NA) {
  fi.start <- abs(fi.start) # force fi.start to not be negative

  highest.f <- sum(apply(x, 1, max))
  fi.start <- min(fi.start, highest.f) # force fi.start to not exceed number of patients

  if (is.null(rownames(x))) stop("Please include rownames on the data table")
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

#' Calculate an exact fragility index for the fisher exact test
#'
#' Exactly alcululate a fragility index using for the fisher exact test using a table.
#' Do a grid search over all FI values <= some value.
#'
#' Do not call this function: use bin.fi with option 'exact' instead
#'
#' @param x a 2x2 matrix with treatments on rows and outcomes on columns
#' @param alpha a string whose value being "exact" leads to an exact calculation
#' @param get.p a function inputting a table and outputting a p-value
#' @param max.f the maximum fragility index considered
#' @param verbose a boolean value for whether to print progress
#' percentages (in increments of roughly ten percent)
#' @return a fragility index
#'
#' @examples
#' x <- matrix(nrow = 2, byrow = TRUE, rgeom(4, 1 / 50))
#' colnames(x) <- c("event", "nonevent")
#' rownames(x) <- c("control", "treatment")
#' get.p <- function(tab) fisher.test(tab)$p.value
#' bin.fi.exact(x, get.p)
bin.fi.exact <- function(x, get.p, alpha = 0.05, max.f = Inf, verbose = FALSE) {
  stop("using bin.fi.exact2 instead") # added when depracating

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

#' Original FI, for binary data
#'
#' @param mat a 2x2 matrix with groups on the rows
#' @param get.p a function which accepts a 2x2 matrix and outputs a p value
#' @param alpha the p-value cutoff to reject, 0.05 by default
#' @param dir a character, either "left", "right", or "both". The default is "both".
#'
#' @examples
#' bin.fi.walsh(matrix(c(100, 96, 13, 28), nrow = 2), function(mat) fisher.test(mat)$p.value)
bin.fi.walsh <- function(mat, get.p, alpha = 0.05, dir='both') {
  start.p <- get.p(mat)
  group.sizes <- apply(mat, 1, min) # smallest entry in each group
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
#' @return the argument of f which satisfies the desired conditions.
#'
#' @examples
#' find_zero(function(x) x - 100 + rnorm(1), fz.verbose = FALSE)
find_zero <- function(f, x.init = 10L, fz.verbose = FALSE, D = 1, gamma = .6,
                      burnin_dur = 16, eps = .1, proj = function(a) a, limits = c(1, 9999999)) {
  # updates proj to take into account the upper and lower limits
  proj. <- function(a) {
    out <- proj(a)
    out <- min(out, limits[2])
    out <- max(out, limits[1])
    return(out)
  }

  # tries to get lowest value with f>=0
  x <- x.init
  y <- f(x)
  history <- c("x" = x, "y" = y, "avg" = NA)
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
      history <- rbind(history, c("x" = x, "y" = y, "avg" = NA))

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
      history <- rbind(history, c("x" = x, "y" = y, "avg" = NA))

      if (old_x == limits[1] & x == limits[1]) break
    }
    if (fz.verbose) print("finished decreasing")

    if (!(x %in% limits)) {
      x <- (y * old_x - old_y * x) / (y - old_y)
      x <- proj.(x)
    }
  }

  # do ruppert-polyak averaging
  ## burn in
  y <- c(f(x[1]))
  if (fz.verbose) print(c(x[1], y[1]))
  history <- rbind(history, c("x" = x, "y" = y, "avg" = NA))

  old_avg <- x[1]
  if (fz.verbose) print("start ruppert polyak averaging burn in")
  for (t in 2:burnin_dur) {
    x[t] <- proj.(x[t - 1] - (D * (t - 1)^(-gamma)) * y[t - 1]) # proj.(a) max(1,round(a)) #########################
    y[t] <- f(x[t])

    if (fz.verbose) print(c(x[t], y[t]))
    history <- rbind(history, c("x" = x[t], "y" = y[t], "avg" = NA))
  }
  old_x <- x[burnin_dur]
  old_y <- y[burnin_dur]

  ## iterate until convergence
  if (fz.verbose) print("start iterating until convergence")

  x <- c(old_x)
  y <- c(old_y)

  if (fz.verbose) print(c(x[1], y[1]))
  history <- rbind(history, c("x" = x, "y" = y, "avg" = NA))

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

#' Convert a 2x2 table to a X,Y data frames
#'
#' @param mat 2x2 numeric matrix
#' @return list containing two one-column data frames X and Y
#'
#' @examples
#' mat_to_xy(matrix(1:4, nrow = 2))
mat_to_xy <- function(mat, event.col = 1) {
  # get names
  r1 <- rownames(mat)[1]
  r2 <- rownames(mat)[2]
  c1 <- colnames(mat)[event.col]
  c2 <- colnames(mat)[3 - event.col]

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
  X <- data.frame(factor(X, levels = rownames(mat))) # c(0,1)))
  Y <- data.frame(factor(Y, levels = colnames(mat))) # c(0,1)))

  rownames(X) <- 1:nrow(X)
  colnames(X) <- "X"
  rownames(Y) <- 1:nrow(Y)
  colnames(Y) <- "Y"

  # return
  return(list("X" = X, "Y" = Y))
}

#' Gets parameters for sampling in ltfu_fi
#' Chooses them so that the mode is po and the 75th
#' quantile is mult*po
get_beta_parameters <- function(po, mult=1.3) {
  alpha <- function(s) s*po+1
  beta <- function(s) s-s*po+1

  get75q <- function(s) qbeta(.875, alpha(s), beta(s))
  shat <- uniroot(function(s) get75q(s)-mult*po, c(0,999999))$root
  return(shat)
}
