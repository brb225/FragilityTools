# note, this file is not well organized, sorry.

# looking at proportion of times the fragility index has unique finite value.
library(tidyverse)
out_30_80 <- vector(mode='list', length=1000)
n <- 500
p1 <- .30
p2 <- .80
for (i in 1:1000) {
  if(i%%100==0) print(i)
  a <- rbinom(1, n, p1)
  c <- rbinom(1, n, p2)
  x <- matrix(nrow=2,byrow=TRUE,c(a, n-a, c, n-c))
  colnames(x) <- c('event', 'nonevent')
  rownames(x) <- c('treatment', 'control')
  out <- bin.fi.incidence(data.table=x, alg='exact')

  out_30_80[[i]] <- out
  #lo <- length(unique(out[!is.infinite(out[,2]), 2]))
  #num_finite <- c(num_finite, lo)
}

get.unique.finite.val.cnt <- function(out) length(unique(out[!is.infinite(out[,2]), 2]))
med.length.sum <- function(out) {
  if (length(unique(out[!is.infinite(out[,2]), 2]))>1) {
    out2 <- out[out[,2]!=out[1,2],]
    out2[1,1]
  } else {
    NA
  }
}
finite.length <- function(out) {
  subset(out, is.infinite(out[,2]))[1,1]
}
third.length <- function(out) {
  if (length(unique(out[,2]))<4) return(NA)
  k <- which.max(is.infinite(out[,2]))
  out[k,1]-out[k-1,1]
}

plyr::laply(out_30_80, .fun=get.unique.finite.val.cnt) %>% table
plyr::laply(out_30_80, .fun=med.length.sum) %>% summary























#### get every trial matrix with 100 patients in each arm ####
# set up functions
get_mat_list <- function(k1, k2) {
  mat.list <- expand.grid(0:k1, 0:k2)

  #if(k1==k2) mat.list <- mat.list[mat.list[,1]<=mat.list[,2],]
  mat.list <- mat.list[mat.list[,1]<=ceiling(k1/2),] # restricts row permutations

  # a <- cbind(mat100.list, (k-mat100.list[,2:1]))
  # da <- duplicated(rbind(a[,1:2], a[,3:4]))
  # da <- da[-(1:nrow(a))] # if transform is duplicated
  # mat100.list <- mat100.list[!da,] # this isn't right.. but idk a better way

  mat.list <- mat.list[sample(1:nrow(mat.list)),] # shuffle to get better time estimate
  mat.list <- apply(mat.list, 1, function(v)
    list(matrix(byrow=TRUE, nrow=2, c(v[1], k1-v[1], v[2], k2-v[2]), dimnames=list(1:2, 1:2))))
  lapply(mat.list, function(l) l[[1]])
}
get_generalized <- function(mat.list) {
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, varlist=c("bin.fi", "bin.fi.exact2", "bin.fi.walsh"), envir=environment())
  out <- pbapply::pblapply(mat.list, bin.fi.incidence, cl=cl)
  parallel::stopCluster(cl)
  return(out)
}
get_walsh <- function(mat.list) {
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, varlist=c("bin.fi", "bin.fi.walsh"), envir=environment())
  out <- pbapply::pblapply(mat.list, function(mat)bin.fi(mat, alg='walsh'), cl=cl)
  parallel::stopCluster(cl)
  return(out)
}

# get matrices and FIs
mat_blah.list <- get_mat_list(20, 150)#(25, 75)#(200,20)#
gen_blah <- get_generalized(mat_blah.list)
walsh_blah <- get_walsh(mat_blah.list)

# find whether orig alg output is not gen FI anywhere
all_bad_orig <- rep(FALSE, length(mat100.list))
for (i in 1:length(mat100.list)) {
  if (!(walsh100[[i]]$FI %in% gen100[[i]][,2])) all_bad_orig[i] <- TRUE
}

nice_score <- function(m) {
  if (m[1,1]<=sum(m[1,])/2) return(m[2,1]<sum(m[2,])/2)
  else return(m[2,1]>sum(m[2,])/2)
}
interesting.inds <- (1:length(mat100.list))[all_bad_orig &
                                              plyr::laply(mat100.list, nice_score) &
                                              plyr::laply(walsh100, function(ll) ll$FI>0)
                                            ]
errors <- plyr::laply(interesting.inds, function(k) min(abs(walsh100[[k]]$FI-gen100[[k]])))
interesting.inds[errors==max(errors)] #603

# get figure 2 replacement (long q support, mostly agrees with q=0)
finite.length <- function(out) {
  subset(out, is.infinite(out[,2]))[1,1]
}
first.length <- function(out) {
  fi0 <- out[1,2]
  k <- which(out[,2]!=fi0)[1]
  out[k,1]
}
finite.cnt <- function(out) {
  v <- unique(out[,2])
  v <- v[!is.infinite(v)]
  length(v)
}
nice_score <- function(m) {
  if (m[1,1]<=sum(m[1,])/2) return(m[2,1]<sum(m[2,])/2)
  else return(m[2,1]>sum(m[2,])/2)
}

mat_list <- mat_blah.list #
gen_list <- gen_blah #
walsh_list <- walsh_blah #
data.frame(id=1:length(mat_blah.list),
           finiteL=plyr::laply(gen_list, finite.length),
           firstL=plyr::laply(gen_list, first.length),
           nice=plyr::laply(mat_blah.list, nice_score),
           num_finite=plyr::laply(gen_list, finite.cnt),
           positive=plyr::laply(walsh_list, function(o)o$FI)>0) %>%
  mutate(ratio=firstL/finiteL) %>%
  filter(num_finite==2, nice, ratio>.75, finiteL>.75) %>%
  arrange(desc(ratio), desc(finiteL)) %>%
  head

# find opposite example---short finite and short near 0
## relies on above code
data.frame(id=1:length(mat_blah.list),
           finiteL=plyr::laply(gen_list, finite.length),
           firstL=plyr::laply(gen_list, first.length),
           nice=plyr::laply(mat_blah.list, nice_score),
           num_finite=plyr::laply(gen_list, finite.cnt),
           positive=plyr::laply(walsh_list, function(o)o$FI)>0) %>%
  mutate(ratio=firstL/finiteL) %>%
  filter(num_finite==2, nice, ratio<.25, finiteL<.5) %>%
  arrange(desc(ratio), desc(finiteL)) %>%
  head #5482

# find nice groups of 3 finite
## relies on code above
data.frame(id=1:length(mat_blah.list),
           num_finite=plyr::laply(gen_list, finite.cnt),
           positive=plyr::laply(walsh_list, function(o)o$FI)>0) %>%
  filter(num_finite==3) %>%
  head #25, 75 -- id=113











# do comparison 1
finitelength100 <- plyr::laply(out100, .fun=finite.length)
mat100.list[finitelength100==0]

# do comparison 2
out100 %>% plyr::laply(third.length) %>% summary #0.4 max from 80,50

# do comparison 3
uniqueval100 <- plyr::laply(out100, .fun=get.unique.finite.val.cnt)
uniqueval100 %>% table
mat100.list[uniqueval100==3] %>% plyr::laply(function(m) fisher.test(m)$p.value) %>% summary
mat100.list[uniqueval100==3] %>% lapply(function(m) incidence.plot(bin.fi.incidence(m)))
#plyr::laply(out100, .fun=med.length.sum) %>% summary






#50, 100 example
mat50.150.list <- get_mat_list(50, 150)
out50.150 <- get_output(mat50.150.list)
plyr::laply(out50.150, .fun=get.unique.finite.val.cnt) %>% table
