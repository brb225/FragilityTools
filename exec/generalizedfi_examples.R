library(FragilityTools) #devtools::load_all()
library(ggplot2)

#### binary: 2x2 ####
get.p <- function(mat) fisher.test(mat)$p.value
get.mat <- function(n, p1, p2) {
  n1 <- rpois(1, n)
  n2 <- rpois(1, n)
  x1 <- rbinom(1, n, p1)
  x2 <- rbinom(1, n, p2)

  mat <- matrix(nrow=2, byrow=TRUE, c(x1, n1-x1, x2, n2-x2))
  rownames(mat) <- 1:2
  colnames(mat) <- 1:2
  return(mat)
}

iters <- 1000

# run code for signal setting
fi <- matrix(NA, ncol=4, nrow=iters)
speed <- matrix(NA, ncol=4, nrow=iters)
mat <- vector(length=iters, mode='list')
set.seed(123456790)
for (i in 1:iters) {
  if (i%%10==0) print(i)
  x <- get.mat(500, .3, .1)
  mat[[i]] <- x
  speed[i,1] <- system.time({
    fi[i,1] <- bin.fi.walsh(x, get.p)$FI
  })[1]
  speed[i,2] <- system.time({
    fi[i,2] <- bin.fi.exact2(x, get.p)$FI
  })[1]
  speed[i,3] <- system.time({
    fi[i,3] <- bin.fi.greedy(x, get.p)$FI
  })[1]
  speed[i,4] <- system.time({
    fi[i,4] <- bin.fi(x, test='fisher', alg='hybrid')$FI # backwards search alg, a little bit slower due to other checks
  })[1]
}

# get output for paper
ggplot2::qplot(fi[,1]-fi[,2])
apply(speed, 2, median)
apply(fi, 2, median)

mean(fi[,1]!=fi[,2])
mean(fi[,3]!=fi[,2])

max(fi[,1]-fi[,2])
max(fi[,3]-fi[,2])

mean(fi[,4]!=fi[,2])

# run code for no signal setting
fi2 <- matrix(NA, ncol=4, nrow=iters)
speed2 <- matrix(NA, ncol=4, nrow=iters)
mat2 <- vector(length=iters, mode='list')
for (i in 1:iters) {
  if (i%%10==0) print(i)
  x <- get.mat(500, .4, .4)
  mat2[[i]] <- x
  speed2[i,1] <- system.time({
    fi2[i,1] <- bin.fi.walsh(x, get.p)$FI
  })[1]
  speed2[i,2] <- system.time({
    fi2[i,2] <- bin.fi.exact2(x, get.p)$FI
  })[1]
  speed2[i,3] <- system.time({
    fi2[i,3] <- bin.fi.greedy(x, get.p)$FI
  })[1]
  speed2[i,4] <- system.time({
    fi2[i,4] <- bin.fi(x, test='fisher', alg='exact')$FI # backwards search alg, a little bit slower due to other checks
  })[1]
}

# get ouput for paper
ggplot2::qplot(fi2[,1]-fi2[,2])
apply(speed2, 2, median)
apply(fi2, 2, median)

mean(fi2[,1]!=fi2[,2])
mean(fi2[,3]!=fi2[,2])

max(fi2[,1]-fi2[,2])
max(fi2[,3]-fi2[,2])

mean(fi2[,4]!=fi2[,2])














#### binary: do avandia meta analysis ####
avandia <- readr::read_csv("exec/generalizedfi_avandia.csv")
avandia <- as.matrix(avandia)#[,-1]

outmeta <- meta::metabin(avandia[,2], avandia[,1], avandia[,4], avandia[,3], sm='OR', method='Peto') # reproducing their output
#outmeta


# get vector of q's to consider
avandia2 <- avandia[,c(2,1,4,3)]
get_phatsl <- function(dat.row) c(rep(1-dat.row[1]/dat.row[2], dat.row[1]), rep(dat.row[1]/dat.row[2], dat.row[2]-dat.row[1]),
                                  rep(1-dat.row[3]/dat.row[4], dat.row[3]), rep(dat.row[3]/dat.row[4], dat.row[4]-dat.row[3]))
phat_sl <- apply(avandia2, 1, get_phatsl)
phat_sl <- Reduce(c, phat_sl)
phat_sl <- sort(unique(phat_sl))

q_vals <- vector(length=length(phat_sl)-1)
for (i in 1:length(q_vals)) {
  q_vals[i] <- (phat_sl[i]+phat_sl[i+1])/2
}
q_vals <- c(0, q_vals)

## get GFI output for all qs
metaout <- vector(mode='list', length=length(q_vals))
library(parallel)
cl <- parallel::makeCluster(7)
avandia_loop <- cbind(id=1:nrow(avandia), avandia2)
metalooptime <- system.time({
  for (i in 1:length(q_vals)) {
    print(i/length(q_vals))
    metaout[[i]] <- binmeta.fi(data=avandia_loop, method='Peto', sm='OR', verbose=TRUE, q=q_vals[i])
  }
}) # user = 183.911 s
parallel::stopCluster(cl)

# arrange into "out" for plotting the per-patient
metaoutmat <- cbind(q_vals, plyr::laply(metaout, function(ll) ll[[1]])) # get q, fi extracted
metaoutmat2 <- metaoutmat

metaoutmat2[-1,1] <- phat_sl[-length(phat_sl)] # get phat_sl back in
metaoutmat2 <- rbind(metaoutmat2, c(phat_sl[length(phat_sl)], Inf))
colnames(metaoutmat2) <- c('min.p', 'fi')

metaoutmat2 <- metaoutmat2[!duplicated(metaoutmat2[,2]),]# remove dupes

ip <- incidence.plot(metaoutmat2, ylab=latex2exp::TeX("Generalized fragility index $GFI_q$"))+theme_bw()
ggsave('meta_sl.eps', ip, device='eps', width = 35, height = 32, units = "mm", dpi=300, scale=2.5)

## get overall likelihood
metaout_zmod <- plyr::llply(metaout, function(ll) ll[[8]]) # get zmod extracted

# get estimated likelihood function
likelihoodratio <- function(Z, Zmod, estOR) {
  loglike <- function(ZZ) {
    # get matrix with each row being a trial's info
    numTrials <- length(unique(ZZ[[1]]))
    trialInfo <- matrix(nrow=numTrials, ncol=4)
    for (t in 1:numTrials) {
      ZZcur <- ZZ[ZZ[[1]]==t,2:3]
      ak <- sum(ZZcur[[1]]==1 & ZZcur[[2]]==1)
      apbk <- sum(ZZcur[[1]]==1)
      cpdk <- sum(ZZcur[[1]]==0)
      apck <- sum(ZZcur[[2]]==1)
      trialInfo[t,] <- c(ak, apbk, cpdk, apck)
    }

    pr <- apply(trialInfo, 1, function(tIR) BiasedUrn::dFNCHypergeo(tIR[1], tIR[2], tIR[3], tIR[4], estOR)) # apply biasedurn to each row
    pr <- log(pr)
    return(sum(pr)) # add lls and return
  }

  Z_ll <- loglike(Z) # get loglike for each Z and Zmod
  Zmod_ll <- loglike(Zmod)

  # return exp diff
  return(exp(Z_ll - Zmod_ll))
}

estOR <- 1.4283 # don't know how to get from meta::metabin

row.to.df <- function(dat.row) { # ripped from front.R
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
all.dfs <- plyr::alply(as.data.frame(avandia_loop), 1, row.to.df)
Z <- dplyr::bind_rows(all.dfs)

overall_lratios <- plyr::laply(metaout_zmod, function(Zmod) likelihoodratio(Z, Zmod, estOR)) # get the overall likelihoods
metaoutmat_overall <- cbind(qq=overall_lratios, FI=plyr::laply(metaout, function(ll) ll[[1]])) # get q, fi extracted

pltol <- ggplot(data.frame(metaoutmat_overall), aes(x=qq, y=FI), log10='x')+
  scale_x_continuous(trans='log10')+annotation_logticks(sides='b')+
  labs(x='Overall likelihood (ratio) of modified data', y=latex2exp::TeX("Generalized fragility index"))+
  geom_point()+theme_bw()
ggsave('meta_sl_overall.eps', pltol, device='eps', width = 35, height = 32, units = "mm", dpi=300, scale=2.5) #width=75












#### binary: e value example ####

library(survival)
pelotas <- haven::read_dta("~/Box Sync/school/research/marty/fragility_index/pelotas/pelotas.dta") # not public, please contact if needed

# clean data
dat <- subset(pelotas, (cdeath==1 & case==1) | case==0) # pretend like other cases dont exist
dat$y <- (dat$case==1)*1
#dat <- subset(dat, milkill %in% c(1,5))
dat$milktype <- as.factor(dat$milkill)

dat <- dat[dat$stratum %in% as.numeric(names(table(dat$stratum)[table(dat$stratum)>1])),]
dat$stratum <- as.factor(dat$stratum) # remove nonrepeated then make factor

# initial p value
m.e.init <- survival::clogit(y~milktype+age+bw+educmoth+strata(stratum),
                      data=dat,
                      method='approximate',
                      init=c(0.98697, 0.33889, 1.26744, 1.12445, -0.03217, -0.00078, 0.00613)
)
summary(m.e.init)

# set up gfi
get.replacements <- function(y, x, rn, Y, X) data.frame("y" = setdiff(unique(Y[[1]]), y)) # recall Y is left unchanged in the algorithm
get.p.val <- function(X,Y) { # returns e value
  # initial p value
  m <- survival::clogit(y~milktype+age+bw+educmoth+strata(stratum),
                        data=cbind(X,Y),
                        method='approximate',
                        init=coef(m.e.init)
  )
  m.lowerlim <- summary(m)$conf.int[3,3]
  ee <- function(r) {
    if (r >= 1) {
      r + sqrt(r*(r-1))
    } else {
      1/2 # arbitrary number < 1
    }
  }
  return(ee(m.lowerlim))
}

X.e <- as.data.frame(dat[,c('milktype', 'age', 'bw', 'educmoth', 'stratum')])
Y.e <- as.data.frame(dat[,'y'])
out.e <- greedy.fi(X.e, Y.e, get.replacements=get.replacements, get.p.val=get.p.val, alpha=1, verbose=TRUE) # non traditional use of inputs

#### gaussian: one sample t test ####
## looking at fake data of focus
set.seed(1234567890)
y <- rnorm(500, mean = 0.1)
p.grid <- c(seq(.01, .5, by = .02), seq(.5, .99, by = .02))
timeee <- system.time({
  #fi.grid <- sapply(p.grid, function(p) unlist(ttest.fi(y, q=p)[c('FI', 'sl')]))
  out.grid.t <- sapply(p.grid, function(p) unlist(ttest.fi(y, q=p)[c('FI', 'sl', 'max.likelihood')]))
})
out.grid.t <- t(out.grid.t)

out.grid.t[,3] <- rev(out.grid.t[,3]) ### hacky fix for now
plt1.t <- ggplot2::qplot(out.grid.t[,3], out.grid.t[,1],
               geom = "line", xlim = c(min(p.grid), max(p.grid)),
               xlab = latex2exp::TeX("Modification likelihood threshold $q$"), ylab = latex2exp::TeX("Generalized fragility index $GFI_q$")
)+theme_bw()
plt2.t <- ggplot(as.data.frame(out.grid.t), aes(x=sl, y=FI), log10='x')+
  scale_x_continuous(trans='log10')+annotation_logticks(sides='b')+
  labs(x='Overall likelihood (ratio) of modified data', y=latex2exp::TeX("Generalized fragility index"))+
  geom_point()+theme_bw()
ggplot2::qplot(out.grid.t[,3], out.grid.t[,2],
               #xlim = c(min(p.grid), max(p.grid)),
               xlab = "Within-patient Likelihood bound", ylab = "Overall modified data likelihood ratio"
)+theme_bw()

ggsave('t_sl.eps', plt1.t, device='eps', width = 35, height = 32, units = "mm", dpi=300, scale=2.5) #width=75
ggsave('t_sl_overall.eps', plt2.t, device='eps', width = 35, height = 32, units = "mm", dpi=300, scale=2.5) #width=75



## simulating for sampling distribution
niteriter <- 100
gfi_t3 <- matrix(nrow=niteriter, ncol=3)
set.seed(1234567489)
cl <- parallel::makeCluster(7)
## look at 1000 samples of GFI calculations
for (iteriter in 1:niteriter) {
  print(paste0(iteriter, 'a'))
  yy <- rnorm(500, mean = 0.1)
  gfi_t3[iteriter,1] <- ttest.fi(yy, q=.05, cl=cl, verbose=TRUE)$FI
  print(paste0(iteriter, 'b'))
  gfi_t3[iteriter,2] <- ttest.fi(yy, q=.5, cl=cl)$FI
  print(paste0(iteriter, 'c'))
  gfi_t3[iteriter,3] <- ttest.fi(yy, q=.95, cl=cl)$FI
}
parallel::stopCluster(cl)

dat.t3 <- data.frame(xx = c(gfi_t3), yy = rep(as.factor(c(.05, .5, .95)) ,each = niteriter))
t3.plots <- ggplot(dat.t3, aes(x=xx, fill=yy)) +
  geom_histogram(alpha=.4, position="identity", bins=50) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2", name='Likelihood threshold q') + #scale_fill_manual(values=c('blue', 'red', 'green'))
  labs(x=latex2exp::TeX("Generalized fragility index $GFI_q$"), y='Count') +
  theme(legend.position=c(.71,.75))
t3.plots
ggsave('gfi_t_samplingdistn.eps', t3.plots, device=cairo_ps, fallback_resolution = 300,
       width = 35, height = 32, units = "mm", dpi=300, scale=2.5) #width=75

apply(gfi_t3, 2, summary)





#### gaussian: modifying covariate example ####
## focusing on fake data
set.seed(123456789)
n.age <- 200
beta.age <- rep(.2, 3)
x.age <- rnorm(n.age)
z.age <- rbinom(n.age, 1, 1/2)
eta.age <- apply(t(beta.age*t(cbind(1,x.age,z.age))),1,sum)
y.age <- rbinom(n.age, 1, binomial()$linkinv(eta.age))

library(parallel)
cl <- parallel::makeCluster(7)

FIrun.out.age <- vector(mode='list', length=10)
time.out.age <- rep(NA, 10)
counter.age <- 0
for (q.maxl in seq(.01, .99, length.out=10)[-c(1:3)]) { # removing the first 3 because they're infinite (faster)
  print(q.maxl)
  counter.age <- counter.age + 1; print(counter.age)
  time.out.age[counter.age] <- system.time({
    FIrun.out.age[[counter.age]] <- glm.gaussian.covariate.fi(cbind(x.age, z.age), y.age, q = q.maxl, cl=cl)
  })[1]
}
parallel::stopCluster(cl)

overall.sl.age <- plyr::laply(FIrun.out.age, function(listt) normal.sl(x.age, listt[[8]][[3]]))
fi.age <- plyr::laply(FIrun.out.age, function(listt) listt$FI)
out.age.allinfo <- data.frame(FI=fi.age, overall=overall.sl.age, q=seq(.01, .99, length.out=10))[-c(1:3),]

plt.age.q <- ggplot(out.age.allinfo, aes(x=q, y=FI))+
  labs(x=latex2exp::TeX("Modification likelihood threshold $q$"),
       y=latex2exp::TeX("Generalized fragility index $GFI_q$"))+
  geom_line()+theme_bw()
plt.age.lr <- ggplot(out.age.allinfo, aes(x=overall, y=FI), log10='x')+
  scale_x_continuous(trans='log10')+annotation_logticks(sides='b')+
  labs(x='Overall likelihood (ratio) of modified data', y=latex2exp::TeX("Generalized fragility index"))+
  geom_point()+theme_bw()

ggsave('age_q.eps', plt.age.q, device='eps', width = 35, height = 32, units = "mm", dpi=300, scale=2.5)
ggsave('age_lr.eps', plt.age.lr, device='eps', width = 35, height = 32, units = "mm", dpi=300, scale=2.5)

## imagine the age is \mu = 50 and \sigma = 10.. then how far do these values actually move?
with(FIrun.out.age[[4]],
     median(abs(new_responses$unique.c.y..out.optim.par.. - old_responses$X.regr...1.)))



## simulating from sampling distribution
niteriter <- 100
gfi_age3 <- matrix(nrow=niteriter, ncol=2)
time_age3 <- matrix(nrow=niteriter, ncol=2)
set.seed(1234567489)
cl <- parallel::makeCluster(7)
## look at 1000 samples of GFI calculations
for (iteriter in 1:niteriter) {
  print(paste0(iteriter, 'a'))

  x.age.s <- rnorm(n.age)
  z.age.s <- rbinom(n.age, 1, 1/2)
  eta.age.s <- apply(t(beta.age*t(cbind(1,x.age.s,z.age.s))),1,sum)
  y.age.s <- rbinom(n.age, 1, binomial()$linkinv(eta.age.s))

  time_age3[iteriter,1] <- system.time({
    gfi_age3[iteriter,1] <- glm.gaussian.covariate.fi(cbind(x.age.s, z.age.s), y.age.s, q = .5, cl=cl, verbose=TRUE)$FI
  })[1]
  print(paste0(iteriter, 'b'))
  time_age3[iteriter,2] <- system.time({
    gfi_age3[iteriter,2] <- glm.gaussian.covariate.fi(cbind(x.age.s, z.age.s), y.age.s, q = .95, cl=cl, verbose=TRUE)$FI
  })[1]
}
parallel::stopCluster(cl)

# dat.t3 <- data.frame(xx = c(gfi_t3), yy = rep(as.factor(c(.05, .5, .95)), each = niteriter))
# t3.plots <- ggplot(dat.t3, aes(x=xx, fill=yy)) +
#   geom_histogram(alpha=.4, position="identity", bins=50) +
#   theme_bw() +
#   scale_fill_brewer(palette = "Dark2", name='Likelihood threshold q') + #scale_fill_manual(values=c('blue', 'red', 'green'))
#   labs(x=latex2exp::TeX("Generalized fragility index $GFI_q$"), y='Count') +
#   theme(legend.position=c(.65,.75))
# t3.plots
# ggsave('gfi_t_samplingdistn.eps', t3.plots, device=cairo_ps, fallback_resolution = 300, width = 35, height = 32, units = "mm", dpi=300, scale=2.5) #width=75
#
# apply(gfi_t3, 2, min)

apply(time_age3, 2, summary)
apply(gfi_age3, 2, summary)











#### poisson: outcome w/ covariates (incomplete) ####
## simulate data
set.seed(123456789)
n.pois <- 500
x.pois <- exp(rnorm(n.pois))
z.pois <- rbinom(n.pois,1,1/2)
beta.pois <- rep(1,3)
eta.pois <- apply(t(beta.pois*t(cbind(1,x.pois,z.pois))),1,sum)
y.pois <- rpois(n.pois, poisson()$linkinv(eta.pois))#eta.pois)#

# get point estimates
library(quantreg)
(out.pois.full <- quantreg::rq(y.pois~x.pois+z.pois, method='fn'))
out.pois.reduced <- quantreg::rq(y.pois~z.pois, method='fn')
anova(out.pois.full, out.pois.reduced, test='rank')

## get gfi

#### time to event ####

library(flexsurv)
plot(flexsurvreg(Surv(time, status) ~ as.factor(lung$age>60), data = lung, dist='weibull')) # event time
plot(flexsurvreg(Surv(time, 3-status) ~ as.factor(lung$age>60), data = lung, dist='weibull')) # censor time

time = lung$time; status = lung$status-1; group = as.factor((lung$age>60)+0);
test = "logrank"; max.likelihood = 0.01; max.step = 0; cl = NULL; alpha = .05;
min.time = 0; max.time = NULL; tau = NULL; verbose = TRUE
only.consider = NULL; dont.consider = c()

# i <- 27; y=YY.red[i, ]; x=XX.red[i, ]; rn=row.names(YY.red)[i]
# get.replacements(y,x,rn,Y,X)

#devtools::load_all()
system.time({
  out.surv1 <- surv.fi(time, status, group, verbose=TRUE, q=0.001) # -44
})
system.time({
  out.surv2 <- surv.fi(time, status, group, verbose=TRUE, q=0.9)  # -3
})
