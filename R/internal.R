#' Calculate an exact fragility index for the fisher exact test
#'
#' Exactly alcululate a fragility index using for the fisher exact test using a table
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
#' x <- matrix(nrow=2,byrow=TRUE,rgeom(4, 1/50))
#' colnames(x) <- c('event', 'nonevent')
#' rownames(x) <- c('control', 'treatment')
#' get.p <- function(tab) fisher.test(tab)$p.value
#' bin.fi.exact(x, get.p)
#'
bin.fi.exact <- function(x, get.p, alpha = 0.05, max.f = Inf, verbose=FALSE) {
  if (is.null(rownames(x))) stop("Please include rownames on the data table")
  # get initial p value to check if frag.index = 0
  pp <- get.p(x)

  # set up
  # get p value after modifying with f values
  get.p.mod <- function(f1, f2, x) {
    mod.x <- matrix(nrow = 2, byrow=TRUE, c(x[1,1]+f1,x[1,2]-f1,
                                            x[2,1]+f2,x[2,2]-f2))
    get.p(mod.x)
  }

  # grid of f's to consider
  f1.vals <- seq(-min(x[1,1], max.f), min(x[1,2], max.f))
  f2.vals <- seq(-min(x[2,1], max.f), min(x[2,2], max.f))

  counter <- 0
  tot <- length(f1.vals)*length(f2.vals)

  # calculate grid of p values from each (f1, f2) pair
  f.mat <- matrix(nrow = length(f1.vals), ncol=length(f2.vals))
  p.mat <- matrix(nrow = length(f1.vals), ncol=length(f2.vals))
  for (i in 1:length(f1.vals)) {
    for (j in 1:length(f2.vals)) {
      counter <- counter+1
      if (counter %% round(tot/10) == 0 & verbose) print (counter/tot)

      f1 <- f1.vals[i]
      f2 <- f2.vals[j]

      f.mat[i,j] <- abs(f1) + abs(f2)
      p.mat[i,j] <- get.p.mod(f1, f2, x)
    }
  }
  rownames(p.mat) <- f1.vals
  colnames(p.mat) <- f2.vals

  # extract best (f1, f2) pair
  ## rearrange data to convenient format
  d <- data.frame(f1.val=rep(row.names(p.mat),ncol(p.mat)),
                  f2.val=rep(colnames(p.mat),each=nrow(p.mat)),
                  f.index=as.vector(f.mat),
                  p.value=as.vector(p.mat),
                  stringsAsFactors = FALSE)

  ## remove outside of max.f consideration
  d <- d[d$f.index <= max.f,]

  ## select (f1,f2)'s which give closest p val which is on other side
  rel.rows <- if(pp < alpha) which(d$p.value >= alpha) else which(d$p.value < alpha)
  if (length(rel.rows)!=0) { # if there are (f1,f2) vals which give the opposite significance
    d <- d[rel.rows,] # only look at (f1,f2) which give opposite
    d <- d[min(d$f.index)==d$f.index,] # only look at smallest f index
  }
  if (pp < alpha) { # picks most extreme p val if more than one with smallest f index
    d <- d[which.max(d$p.value),]
  } else {
    d <- d[which.min(d$p.value),]
  }

  row.frag.index <- cbind(d, 'reverse'= (pp >= alpha))
  colnames(row.frag.index)[1:2] <- rownames(x)
  #print(row.frag.index)

  rev <- (pp>=alpha)
  # get in same format as greedy.fi
  return(list('FI'=row.frag.index[1,3]*(-1)^rev, 'p_value_sequence'=c(pp, row.frag.index[1,4]), 'reverse'=rev,
              'num_patients'=sum(x), 'patients'=row.frag.index[,1:2], 'old_responses'=NA, 'new_responses'=NA))
}

#' Original FI, for binary data
#'
#' @param mat a 2x2 matrix with groups on the rows
#' @param get.p a function which accepts a 2x2 matrix and outputs a p value
#' @param alpha the p-value cutoff to reject, 0.05 by default
#'
#' @examples
#' bin.fi.walsh(matrix(c(100,96,13,28),nrow=2), function(mat)fisher.test(mat)$p.value)
#'
bin.fi.walsh <- function(mat, get.p, alpha=0.05) {
  start.p <- get.p(mat)
  small.gr <- which.min(apply(mat, 1, sum))

  get.mod.p <- function(f) {
    mod.mat <- mat
    mod.mat[small.gr,] <- mod.mat[small.gr,]+c(-1,1)*f
    return(get.p(mod.mat))
  }

  end.p <- NA
  for (f in 1:max(mat[small.gr,])) {
    if (-mat[small.gr,2] <= f && f <= mat[small.gr,1]) {
      p1 <- get.mod.p(f)
      if ((p1-alpha)*(start.p-alpha)<0) {
        end.p <- p1
        break
      }
    }
    if (-mat[small.gr,2] <= -f && -f <= mat[small.gr,1]) {
      p2 <- get.mod.p(-f)
      if ((p2-alpha)*(start.p-alpha)<0) {
        end.p <- p2
        break
      }
    }
  }
  if (is.na(end.p)) f <- Inf

  return(list('FI'=abs(f)*(-1)^(start.p>=alpha),'p_value_sequence'=c(start.p,end.p),
              'reverse'=(start.p>=alpha),'num_patients'=sum(mat),
              'patients'=NA,'old_responses'=NA,'new_responses'=NA))
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
#' find_zero(function(x) x-100+rnorm(1), fz.verbose=FALSE)
#'
find_zero <- function(f, x.init=10L, fz.verbose=FALSE, D=1, gamma=.6,
                      burnin_dur=16, eps=.1, proj=function(a)a, limits=c(1,9999999)) {
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
  history <- c('x'=x, 'y'=y, 'avg'=NA)
  if (fz.verbose) print(c(x, y))

  ## init
  if (y<0) {
    # increase exponentially until changing sign of f
    while(y<0) {
      old_y <- y
      old_x <- x

      x <- proj.(2*old_x)
      y <- f(x)

      if (fz.verbose) print(c(x, y))
      history <- rbind(history, c('x'=x, 'y'=y, 'avg'=NA))

      if (old_x==limits[2] & x==limits[2]) break
    }
    if (fz.verbose) print('finished increasing')

    x <- (y*old_x-old_y*x)/(y-old_y)
    x <- proj.(x)
  } else if (y==0) {
    x <- x
  } else {
    # decrease exponentially until changing sign of f
    while(y>0) {
      old_y <- y
      old_x <- x

      x <- proj.(old_x/2)
      y <- f(x)

      if (fz.verbose) print(c(x, y))
      history <- rbind(history, c('x'=x, 'y'=y, 'avg'=NA))

      if (old_x==limits[1] & x==limits[1]) break
    }
    if (fz.verbose) print('finished decreasing')

    if (!(x %in% limits)) {
      x <- (y*old_x-old_y*x)/(y-old_y)
      x <- proj.(x)
    }
  }

  # do ruppert-polyak averaging
  ## burn in
  y <- c(f(x[1]))
  if (fz.verbose) print(c(x[1], y[1]))
  history <- rbind(history, c('x'=x, 'y'=y, 'avg'=NA))

  old_avg <- x[1]
  if (fz.verbose) print('start ruppert polyak averaging burn in')
  for (t in 2:burnin_dur) {
    x[t] <- proj.(x[t-1]-(D*(t-1)^(-gamma))*y[t-1]) #proj.(a) max(1,round(a)) #########################
    y[t] <- f(x[t])

    if (fz.verbose) print(c(x[t], y[t]))
    history <- rbind(history, c('x'=x[t], 'y'=y[t], 'avg'=NA))
  }
  old_x <- x[burnin_dur]
  old_y <- y[burnin_dur]

  ## iterate until convergence
  if (fz.verbose) print('start iterating until convergence')

  x <- c(old_x)
  y <- c(old_y)

  if (fz.verbose) print(c(x[1], y[1]))
  history <- rbind(history, c('x'=x, 'y'=y, 'avg'=NA))

  t <- 2#burnin_dur %/% 2
  old_avg <- x[1]
  while(TRUE) {
    x[t] <- proj.(x[t-1]-(D*(t-1)^(-gamma))*y[t-1])
    y[t] <- f(x[t])

    if (fz.verbose) print(c(x[t], y[t], mean(x)))
    history <- rbind(history, c('x'=x[t], 'y'=y[t], 'avg'=mean(x)))

    if (abs(old_avg - mean(x)) < eps) break
    old_avg <- mean(x)
    t <- t+1
  }
  if (fz.verbose) print('finished')

  final_x <- proj.(mean(x))
  return(list('x'=final_x, 'history'=history))
}

#' Convert a 2x2 table to a X,Y data frames
#'
#' @param mat 2x2 numeric matrix
#' @return list containing two one-column data frames X and Y
#'
#' @examples
#' mat_to_xy(matrix(1:4, nrow=2))
#'
mat_to_xy <- function(mat, event.col=1) {
  # get info out
  xx <- mat[1,event.col]
  yy <- mat[2,event.col]

  m <- sum(mat[1,])
  n <- sum(mat[2,])

  # get sequences, taking into account that some are empty
  seq11 <- 1:xx
  if (xx<1) seq11 <- c()
  seq12 <- (xx+1):m
  if (m<xx+1) seq12 <- c()

  seq21 <- m+(1:yy)
  if (yy<1) seq21 <- c()
  seq22 <- m+((yy+1):n)
  if (n<yy+1) seq22 <- c()

  # init
  X <- c()
  Y <- c()

  # loop through and fill entries
  for (i in seq11) {
    X <- c(X,0)
    Y <- c(Y,1)
  }
  for (i in seq12) {
    X <- c(X,0)
    Y <- c(Y,0)
  }
  for (i in seq21) {
    X <- c(X,1)
    Y <- c(Y,1)
  }
  for (i in seq22) {
    X <- c(X,1)
    Y <- c(Y,0)
  }

  X <- data.frame(factor(X, levels=c(0,1)))
  Y <- data.frame(factor(Y, levels=c(0,1)))

  rownames(X) <- 1:nrow(X); colnames(X) <- 'X';
  rownames(Y) <- 1:nrow(Y); colnames(Y) <- 'Y';

  return(list('X'=X,'Y'=Y))
}
