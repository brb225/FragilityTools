library(FragilityTools)

#### motivating example ####
mat <- matrix(nrow=2,ncol=2,byrow=TRUE,c(100,300,15,15))
fisher.test(mat)$p.value
bin.fi(mat)
bin.fi(mat, r=1/2, nsim=100, verbose=TRUE)$FI # 24




#### voter example ####
fn <- function(n) pbinom(200, n, 1/1000)

fl_vote_totals <- c(2912790,	2912253,	16415,	97488,	562,	2281,	17484,	622,	1371,	1804,	34,	6) # numbers from florida website
turnout <- .7 # from FL 2000 presidential vote page; https://www.dos.myflorida.com/elections/data-statistics/elections-data/voter-turnout/

fl_nonvoters <- sum(fl_vote_totals)/turnout - sum(fl_vote_totals[c(1,2)])

us_vote_totals <- c(50999987, 50456002, 92875537) # gore, bush, nonvoter
us_eligible_voters <- sum(us_vote_totals)

p <- fl_nonvoters/us_eligible_voters
GFInon <- 538

MGFI <- ceiling(GFInon/p)

phyper(GFInon, fl_nonvoters, us_eligible_voters-fl_nonvoters, MGFI) # around 1/2
#pbinom(GFInon, MGFI, p) # around 1/2




#### intro to NHEFS ####
library(ggplot2)
set.seed(0123456789)
theme_set(theme_bw())

nhefs <- readr::read_csv("exec/nhefs.csv")
nhefs <- as.data.frame(nhefs)
dat <- nhefs[,c('death', 'qsmk', 'smokeyrs', 'age', 'wt71', 'smokeintensity', 'sex', 'race', 'education', 'exercise', 'active')]

mod.form <- death ~ qsmk + smokeyrs +
  age + wt71 + I(wt71*wt71) + smokeintensity +
  sex + race + as.factor(education) +
  as.factor(exercise) + as.factor(active)
mod <- glm(mod.form, dat = dat, family=binomial())
summary(mod)
start.coef <- coef(mod)

with(dat, ggplot2::qplot(smokeyrs, death, data=dat, col=as.factor(qsmk))+geom_smooth())






#### traditional FI ####
crosstab <- table(dat[,c('qsmk', 'death')])[2:1,2:1]
crosstab

fisher.test(crosstab)
bin.fi(crosstab, test='fisher', alg='exact')

incidence.plot(bin.fi.incidence(crosstab, test='fisher', alg='exact'))
crosstab[1,1]/sum(crosstab)

bin.fi(crosstab, test='fisher', verbose=TRUE, r=.5)

#### generalized FIs ####
get_viz <- function(out, qpr) {
  var.to.viz <- predict(mod)

  selected.pats <-  (1:nrow(dat)) %in% as.numeric(out$patients)
  cell00 <- dat$death==0 & dat$qsmk==0
  cell01 <- dat$death==0 & dat$qsmk==1
  cell10 <- dat$death==1 & dat$qsmk==0
  cell11 <- dat$death==1 & dat$qsmk==1

  ul_txt <- paste0('q=',qpr)

  datt <- dat
  datt$var.to.viz <- var.to.viz
  datt <- rbind(cbind(datt, 'type'='All'),
                cbind(datt[selected.pats & cell00, ], 'type'=rep('Modified 00', sum(selected.pats & cell00))),
                #cbind(datt[selected.pats & cell01, ], 'type'=rep('Modified; quit smoking, died', sum(selected.pats & cell01))),
                #cbind(datt[selected.pats & cell10, ], 'type'=rep('Modified; continued smoking, survived', sum(selected.pats & cell10)))#,
                cbind(datt[selected.pats & cell11, ], 'type'=rep('Modified 11', sum(selected.pats & cell11)))
  )
  datt$type <- factor(datt$type, levels=c('All',
                                          'Modified; quit smoking, died',
                                          'Modified; continued smoking, survived'))#'Modified 00', 'Modified 01', 'Modified 10', 'Modified 11'))
  ggplot(datt, aes(x=var.to.viz, fill=type)) +
    geom_histogram(position="identity") +
    scale_fill_manual(name = "Patients",
                      values=c("grey", "red", "blue"), #"black", "red", "purple", "blue"),
                      drop=F
    ) +
    xlab('')+#xlab(paste0('Number of years smoked (q=',qpr,')')) +
    ylab('')+#ylab('Count')+
    theme_void()+ #theme_bw()
    annotate("text", label= ul_txt, size=5, x=50, y=120)#x=5, y=70)
}
get.p <- function(X, Y) {
  modd <- glm(mod.form, data=cbind(X,Y), family=binomial(), start=start.coef)
  return(summary(modd)$coefficients[2,4]) # p for qsmk
}
get.replacements <- function(y, x, rn, Y, X) data.frame(Y=setdiff(unique(Y[[1]]), y))

out <- greedy.fi(within(dat, rm('death')),
                 dat[,'death',drop=FALSE],
                 get.replacements, get.p, verbose=TRUE)
get_viz(out, 0)


# get material for q>0
mod.repl <- glm(death~qsmk+smokeyrs, family=binomial(), data=dat)
get.replacements2 <- function(y, x, rn, Y, X, q) {
  datmod <- data.frame(Y=setdiff(unique(Y[[1]]), y))

  # remove if not sufficiently likely
  if (y==1 & q > 1-predict(mod.repl, newdata=x, type='response')) datmod <- data.frame(Y=numeric())
  if (y==0 & q > predict(mod.repl, newdata=x, type='response')) datmod <- data.frame(Y=numeric())
  return(datmod)
}
get.repl2 <- function(q) return(function(y, x, rn, Y, X) get.replacements2(y, x, rn, Y, X, q))

q.grid <- c(0, .25, .57, .6, .8, .9)
plt.grid <- vector('list', length=length(q.grid))
gfi.grid <- rep(NA_real_, length(q.grid))
for (ind in 1:length(q.grid)) {
  q <- q.grid[ind]
  out <- greedy.fi(within(dat, rm('death')),
                   dat[,'death',drop=FALSE],
                   get.repl2(q), get.p, verbose=TRUE)
  gfi.grid[ind] <- out$FI
  plt.grid[[ind]] <- get_viz(out, q=q)
}

library(patchwork)
plt2by2 <- plt.grid[[1]] + plt.grid[[2]] + plt.grid[[3]] +
  plt.grid[[4]] + plt.grid[[5]] + plt.grid[[6]] +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
#gridExtra::grid.arrange(plt.grid[[1]], plt.grid[[3]],
#                        plt.grid[[4]], plt.grid[[6]], ncol=2,
#                        common.legend=TRUE)

ggsave('hnefs_2by2.eps', plt2by2, device='eps', width = 300, height = 200, units = "mm", dpi=300, scale=.6) #width=75


plt <- ggplot(data.frame(q=q.grid, FI=gfi.grid), aes(x=q, y=FI))+
  labs(x='Modification likelihood threshold q', y=latex2exp::TeX("Generalized fragility index $GFI_q$"))+
  geom_line()+geom_point()+theme_bw()
plt

#### Stochastic FIs ####
q.grid <- c(0, .25, .57, .6, .8, .9)
r.grid <- c(0, 1/4, 1/2, 3/4)
FI.mat <- matrix(nrow=length(q.grid), ncol=length(r.grid))
rownames(FI.mat) <- q.grid; colnames(FI.mat) <- r.grid

alg.init <- matrix(10, nrow=length(q.grid), ncol=length(r.grid))
alg.init[,2] <- c(174, 316, 530, 615, 990, 1470) # inital guesses from rough previous run
alg.init[,3] <- c(205, 340, 610, 690, 1060, 1510)
alg.init[,4] <- c(235, 430, 750, 760, 1160, 1600)

set.seed(1234567890)
for (ind.r in 1:length(r.grid)) {
  for (ind.q in 1:length(q.grid)) {
    print(c(ind.r, ind.q))
    q <- q.grid[ind.q]
    r <- r.grid[ind.r]
    FI.mat[ind.q, ind.r] <- stochastic.fi(dat[,c('qsmk', 'smokeyrs')],
                                          dat[,'death',drop=FALSE],
                                          get.repl2(q), get.p,
                                          r = r, qfi.init=alg.init[ind.q, ind.r],
                                          init.step = alg.init[ind.q, ind.r]==10,
                                          verbose=TRUE, D=150)$FI
  }
}

# R crashed, so i'm manually pulling in the output of the above (it's the same)
FI.mat[1:24] <- -c(10,11,11,15,25,30,166,303,540,604,975,1458,209,368,612,695,1076,1517,
                   249,420,706,755,1154,1569)
FIlong <- reshape2::melt(FI.mat)
colnames(FIlong) <- c('q', 'r', 'FI')
qr.plot <- ggplot(FIlong, aes(x=q, y=FI, group=as.factor(r))) +
  geom_line(aes(color=as.factor(r)))+
  geom_point(aes(color=as.factor(r)))+
  theme_bw()+
  labs(x = latex2exp::TeX("Modification likelihood threshold $q$"),
       y = latex2exp::TeX("Stochastic generalized fragility index $SGFI_{q,r}$"),
       colour = latex2exp::TeX("Stochastic threshold $r$"))+
  theme(legend.position=c(.2,.3))
ggsave('qr_plot.eps', qr.plot, device='eps', width = 300, height = 200, units = "mm", dpi=300, scale=.6) #width=75

