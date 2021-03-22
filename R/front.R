#' LTFU-aware fragility index calculator
#'
#' @param mat a 2x2 real non-negative integer matrix with control/treatment on the rows and
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
#' out <- ltfu_fi(mat, ltfu, function(m) fisher.test(m)$p.value)
ltfu_fi <- function(mat, ltfu, get.p,
                    q=0, alpha=.05,
                    sampling.control=NULL,
                    xlab='Control group LTFU event rate', ylab='Treatment group LTFU event rate',
                    fig.size=1.1, gradientn.scale=.99) {
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
    prs <- prs-log(sum(exp(prs))) # the division shouldn't be necessary lol
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
  print(
    ggplot(data=plt.dat[order(plt.dat$mode),],
           aes(x=g1/ltfu[1], y=g2/ltfu[2], fill=sig)
    )+
      geom_tile(aes(color=log_prs), size=fig.size)+#, width=ww, height=hh), size=ss)+#width=.182, height=.182), size=1.5)+
      scale_color_gradientn(name='Posterior \nprobability',
                            colours=c("white","blue","black"),
                            breaks=c(maxpr/2, maxpr))+#n.breaks=3)+
      scale_fill_manual(values=c("#999999", "#E69F00"),
                        name="Statistical \nsignificance",
                        labels=c(paste0("p <= ", plt.alpha),
                                 paste0("p < ", plt.alpha)
                                 ),
                        drop=FALSE)+
      labs(x=xlab, y=ylab)+
      theme_bw()+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(legend.position='bottom')
  )

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

  return(list('FI'=FI.mat, 'info'=FI.info, 'reverse.p'=pp.sig, 'imputation'=expected)) # use the negative if reverse convention
}

#' Plot the output of bin.fi.incidence
#'
#' @param out the output of bin.fi.incidence
#'
#' @return a plot visualizing out
#'
#' @examples
#' x <- matrix(nrow=2,byrow=TRUE,c(5, 100, 30, 70),
#' dimnames = list(c('treatment', 'control'), c('event', 'nonevent')))
#' out <- bin.fi.incidence(data.table=x, alg='exact')
#' incidence.plot(out)
incidence.plot <- function(out) {
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
         y = latex2exp::TeX("Incidence fragility index $FI_q$")) + # put in latex2exp::TeX
    xlim(c(0, 1)) +
    scale_y_continuous(limits = plt.lims)
}

#' Loop over bin.fi for relevant grid of min.p values
#'
#' @param data.table a 2x2 non negative integer matrix representing the output of a clinical trial.
#' The rows correspond to treatment/control and the columns to event/nonevent.
#' @param X a data frame of covariates
#' @param Y a data frame of responses
#' @param alg a string specifying the FI algorithm, 'greedy' (default), 'exact', 'original', or 'original.bothdir'
#' @param test a string specifying the test, defaulting to 'fisher' for the Fisher exact test.
#' An alternative is 'fisher.midp' and 'pearson.chisq'.
#' @param fi.start the starting fragility index considered in the exact algorithm.
#' @param alpha a number for the size of test
#' @param verbose a logical value for whether to print status updates while running
#' @param warm.start If using the exact algorithm without specifying a `max.f` argument,
#' a greedy algorithm will be ran first to upper bound the search for the exact algorithm.
#' Set to FALSE to only run the exact algorithm.
#' @param delta the noninferiority margin for when test='ni.normal'
#'
#' @return a matrix containing the p values at which there's a changepoint and the corresponding
#' values of the incidence fragility index
#'
#' @example
#' x <- matrix(nrow=2,byrow=TRUE,c(5, 100, 30, 70),
#' dimnames = list(c('treatment', 'control'), c('event', 'nonevent')))
#' out <- bin.fi.incidence(data.table=x, alg='exact')
#' @export
bin.fi.incidence <- function(data.table = NULL, X = NULL, Y = NULL, alg = "exact",
                        test = "fisher", fi.start = NULL, alpha = 0.05, verbose = FALSE,
                        warm.start = TRUE, delta = NULL) {
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
        data.table = data.table, X = X, Y = Y, alg = alg,
        test = test,
        fi.start = fi.start,
        alpha = alpha, verbose = verbose,
        warm.start = warm.start, delta = delta, min.p = minp.grid[ind]
      )$FI)
    }
    fi.start <- fi.grid[ind]
    ind <- ind + 1
  }
  return(cbind("min.p" = c(0, minp.grid[-1] - eps), "fi" = fi.grid))
}

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
#' @param alg a string specifying the FI algorithm, 'greedy' (default), 'exact', 'original', or 'original.bothdir'
#' @param test a string specifying the test, defaulting to 'fisher' for the Fisher exact test.
#' An alternative is 'fisher.midp' and 'pearson.chisq'.
#' @param fi.start the starting fragility index considered in the exact algorithm.
#' @param warm.start If using the exact algorithm without specifying a `max.f` argument,
#' a greedy algorithm will be ran first to upper bound the search for the exact algorithm.
#' Set to FALSE to only run the exact algorithm.
#' @param alpha a number for the size of test
#' @param verbose a logical value for whether to print status updates while running
#' @param delta the noninferiority margin for when test='ni.normal'
#' @param min.p the minimum probability of allowable outcome changes, defaults to 0 (allowing all changes).
#' Note that alg='original' or alg='original.bothdir' does not support this option.
#' @return the fragility index and other values, depending on which algorithm is specified.
#'
#' @examples
#' x <- matrix(nrow = 2, byrow = TRUE, c(5, 100, 30, 70))
#' colnames(x) <- c("event", "nonevent")
#' rownames(x) <- c("treatment", "control")
#' p.grid <- seq(0, 1, by = .01)
#' fi.grid <- sapply(p.grid, function(p) bin.fi(data.table = x, alg = "greedy", min.p = p, verbose = FALSE)$FI)
#' ggplot2::qplot(p.grid[fi.grid < Inf], fi.grid[fi.grid < Inf],
#'   xlab = "Likelihood bound",
#'   ylab = "Fragility index", xlim = c(0, 1), main = "Fisher exact test fragility indices (n=205)"
#' )
#' @export
bin.fi <- function(data.table = NULL, X = NULL, Y = NULL,
                   alg = "greedy", test = "fisher", fi.start = NULL, alpha = .05,
                   verbose = FALSE, warm.start = TRUE, delta = NULL, min.p = 0) {
  # for greedy algorithm or to warm start exact algorithm
  if (alg == "greedy") { # | (max.f==Inf & warm.start & alg=='exact')) {
    # convert data.table to X,Y
    if (is.null(Y)) {
      dat <- mat_to_xy(data.table, event.col = 2)
      X <- dat[[1]]
      colnames(X) <- "group"
      Y <- dat[[2]]
      colnames(Y) <- "resp"
    }
    # get p value function
    if (test == "fisher") {
      get.p.val <- function(X, Y) fisher.test(table(cbind(X, Y)))$p.value
    } else if (test == 'greater.fisher') {
      get.p <- function(tab) fisher.test(table(cbind(X, Y)), alternative = 'g')$p.value
    } else if (test == 'less.fisher') {
      get.p <- function(tab) fisher.test(table(cbind(X, Y)), alternative = 'l')$p.value
    } else if (test == 'cambell') {
      get.p.val <- function(X,Y) {
        x <- table(cbind(X, Y))
        r.marg <- rowSums(x)
        c.marg <- colSums(x)
        expected.mat <- outer(r.marg, c.marg, '*')/sum(x)

        if (min(expected.mat)>=1) {
          #N-1 pearson chi square
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
      get.p.val <- function(X, Y) {
        m <- sum(Y[[1]] == Y[[1]][1])
        n <- nrow(Y) - m
        k <- sum(X[[1]] == X[[1]][1])
        x <- sum(X[[1]] == X[[1]][1] & Y[[1]] == Y[[1]][1])

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
      get.p.val <- function(X, Y) {
        p <- stats::chisq.test(table(cbind(X, Y)))$p.value
        if (is.nan(p)) return(1)  # if a marginal is 0, then table is rank one
        return(p)
      }
    } else if (test == "chisq.prop") {
      get.p.val <- function(X, Y) {
        n1 <- sum(X[[1]] == 0)
        n2 <- sum(X[[1]] == 1)
        p1 <- sum(X[[1]] == 0 & Y[[1]] == 1) / n1
        p2 <- sum(X[[1]] == 1 & Y[[1]] == 1) / n2

        pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
        ts <- (p1 - p2) / sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
        p_value <- 2 * pnorm(abs(ts), lower.tail = FALSE)
        return(ifelse(is.nan(p_value), 1, p_value))
      }
    } else if (test == "ni.normal") {
      get.p.val <- function(X, Y) {
        n1 <- sum(X[[1]] == 0)
        n2 <- sum(X[[1]] == 1)
        p1 <- sum(X[[1]] == 0 & Y[[1]] == 1) / n1
        p2 <- sum(X[[1]] == 1 & Y[[1]] == 1) / n2

        ts <- (p1 - p2 - delta) / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
        if (is.nan(ts)) ts <- Inf # if p1=p2=0 or p1=p2=1, then I should reject!

        return(pnorm(ts, lower.tail = FALSE))
      }
    } else {
      stop("Please select an available test option")
    }

    # get replacements function.. taking into account min.p
    get.replacements <- function(y, x, rn, Y, X) data.frame("outcome" = setdiff(unique(Y[[1]]), y)) # recall Y is left unchanged in the algorithm

    # get only.consider according to mid.p
    df <- cbind(X, Y)
    df.short <- df[!duplicated(df), ]
    get_prob_from_row <- function(rr) {
      dfr <- df[df[[1]] == rr[[1]], ]
      sum(dfr[[2]] != rr[[2]]) / nrow(dfr)
    }
    df.short$prob <- apply(df.short, 1, get_prob_from_row)
    # if (verbose) print(df.short)
    df <- plyr::join(df, df.short, by = c(colnames(X), colnames(Y)))
    only.consider <- rownames(X)[df$prob >= min.p]

    # run alg
    out <- greedy.fi(
      X = X, Y = Y,
      get.replacements = get.replacements, only.consider = only.consider,
      get.p.val = get.p.val, alpha = alpha, verbose = verbose
    )

    # post process "patients" if started as table?
  } # end df format

  # for exact or walsh algorithm (both want a table, not data frame format)
  if (alg == "exact" || alg == "original" || alg == 'original.bothdir') {
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

    # run alg if walsh or warmstarting exact
    start.reversed<-FALSE
    start.p <- NA
    if (alg == "original" | alg=="original.bothdir" | warm.start & alg == "exact" & (is.null(fi.start) || is.infinite(fi.start))) {
      if(alg=='original.bothdir' | alg=='exact') dir <- 'both'
      if(alg=='original') dir <- ifelse(get.p(data.table)<alpha, 'left', 'right')

      #print('start warm start') ###
      out <- bin.fi.walsh(mat = data.table, get.p = get.p, alpha = alpha, dir = dir)
      #print('finish warm start') ###

      out.fi <- abs(out$FI)
      if (!is.infinite(out.fi)) {
        start.reversed <- TRUE
        start.p <- out$p_value_sequence[2]
      }
      if (is.infinite(out.fi)) out.fi <- 1 # bin.fi.exact2 cant handle inf start.fi
    }

    if (alg == "exact") {
      # get can.vary according to min.p
      mat.of.probs <- data.table
      mat.of.probs <- mat.of.probs[, c(2, 1)] / apply(mat.of.probs, 1, sum)
      can.vary <- mat.of.probs >= min.p

      # get output
      if (warm.start) {
        if (is.null(fi.start)) {
          out <- bin.fi.exact2(
            x = data.table, alpha = alpha, get.p = get.p,
            fi.start = out.fi, can.vary = can.vary, start.reversed=start.reversed, start.p=start.p
          )
        } else {
          out <- bin.fi.exact2(
            x = data.table, alpha = alpha, get.p = get.p,
            fi.start = fi.start, can.vary = can.vary, start.reversed=start.reversed, start.p=start.p
          )
        }
      } else {
        out <- bin.fi.exact2(
          x = data.table, alpha = alpha, get.p = get.p,
          fi.start = fi.start, can.vary = can.vary, start.reversed=start.reversed, start.p=start.p
        )
      }
    }
  } # end matrix format
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
#'   mat <- matrix(nrow = 2, ncol = 2)
#'   mat[1, 1] <- sum((x == 0) & (y == 0))
#'   mat[1, 2] <- sum((x == 0) & (y == 1))
#'   mat[2, 1] <- sum((x == 1) & (y == 0))
#'   mat[2, 2] <- sum((x == 1) & (y == 1))
#'   fisher.test(mat)$p.value
#' }
#' get.replacements <- list(function(y, x, rn, Y, X) setdiff(unique(Y[[1]]), y))
#'
#' ss <- general.fi.samplesize(
#'   min.fi = 10, min.power = .8, get.p.val = get.p.val, niters = 1,
#'   get.replacements = get.replacements, get.sample = draw.binom, nsim = 25, verbose = TRUE
#' )
#' @export
general.fi.samplesize <- function(min.fi = 10, min.power = .8,
                                  sample_size_init_power = 100L, sample_size_init_fi = NULL,
                                  get.p.val, get.replacements, get.sample, gamma = 0.6, niters = 50,
                                  cl = NULL, verbose = FALSE, alpha = .05, tau = 1 / 2, nsim = 30, eps = .1, algorithm = "greedy") {
  # helper functions
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

    sample_size1 <- ceiling(quantile(ss.s, .5, names = FALSE))
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
#'   mat <- matrix(nrow = 2, ncol = 2)
#'   mat[1, 1] <- sum((x == 0) & (y == 0))
#'   mat[1, 2] <- sum((x == 0) & (y == 1))
#'   mat[2, 1] <- sum((x == 1) & (y == 0))
#'   mat[2, 2] <- sum((x == 1) & (y == 1))
#'   fisher.test(mat)$p.value
#' }
#' get.replacements <- list(function(y, x, rn, Y, X) setdiff(unique(Y[[1]]), y))
#' get.sample <- draw.binom
#'
#' out <- min.fi.curve(
#'   min.fi = seq(0, 10, by = 5), get.p.val = get.p.val, get.replacements = get.replacements,
#'   get.sample = get.sample, niters = 1, algorithm = "greedy"
#' )
#' @export
min.fi.curve <- function(min.fi, get.p.val, get.replacements, get.sample, cl = NULL,
                         min.power = .8, alpha = .05, verbose = FALSE, niters = 5,
                         sample_size_init_power = 10L, sample_size_init_fi = NULL,
                         nsim = 30, eps = .1, tau = 1 / 2, algorithm = "greedy", gamma = .6) {
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
#' get.rejection.rates(get.p.val, NULL, get.sample.null, get.sample.alt, 5, 100, algorithm = "walsh")
#' @export
get.rejection.rates <- function(get.p.val, get.replacements = NULL, get.sample.null = NULL, get.sample.alt = NULL,
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

#' Get equivalent parameters in traditional power calculation
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
#' equivalent_usual_parameters(0.05, .8, .06, 850, get.sample.alt.f, get.p.val = get.p.val)
#' @export
equivalent_usual_parameters <- function(alpha, pi, delta, n,
                                        get.sample.alt.f, get.p.val,
                                        verbose = FALSE, cl = NULL, limits = c(0, 1),
                                        nsim = 10000, mc.iters = 100, delta.iters = 100) {
  ## fix n. select 2 of alpha, pi, delta. find the value of the other that gives n
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
