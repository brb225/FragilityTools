set.seed(1234567890)
#### do calculation ####
# power = .8, alpha=.05
print('get sample sizes for plot')
tictoc::tic()
out_pow.8_tau.5 <- general.fi.sscurve(phi.grid, min.power=.8, alpha=0.05, tau=.5,
                                get.p.val=get.p.val, get.replacements=get.replacements, get.sample=get.sample.alt,
                                gamma=0.3, nsim=nsim_low, eps=eps, niters=niters,
                                verbose=FALSE, cl=cl, algorithm='walsh')
tictoc::toc()

print(out_pow.8_tau.5)

# power = .7, alpha=.05
tictoc::tic()
out_pow.9_tau.5 <- general.fi.samplesize(min.fi=NULL, min.power=.9, alpha=0.05, tau=.5,
                                         get.p.val=get.p.val, get.replacements=get.replacements, get.sample=get.sample.alt,
                                         gamma=.3, nsim=nsim_low, eps=eps, niters=niters,
                                         verbose=FALSE, cl=cl, algorithm='walsh')
tictoc::toc()

# power = .8, alpha=.01
tictoc::tic()
out_pow.8_tau.8 <- general.fi.sscurve(phi.grid, min.power=.8, alpha=0.05, tau=.2,
                                get.p.val=get.p.val, get.replacements=get.replacements, get.sample=get.sample.alt,
                                gamma=.3, nsim=nsim_high, eps=eps, niters=niters,
                                verbose=FALSE, cl=cl, algorithm='walsh')
tictoc::toc()

# power = .7, alpha=.01
tictoc::tic()
out_pow.9_tau.8 <- general.fi.samplesize(min.fi=NULL, min.power=.9, alpha=0.05, tau=.2,
                                         get.p.val=get.p.val, get.replacements=get.replacements, get.sample=get.sample.alt,
                                         gamma=.3, nsim=nsim_high, eps=eps, niters=niters,
                                         verbose=FALSE, cl=cl, algorithm='walsh')
tictoc::toc()


#### plot #####
print('make plot')
out <- rbind(cbind(min.power=.8,tau=0.5, out_pow.8_tau.5),
             cbind(min.power=.9,tau=0.5, cbind(out_pow.8_tau.5[1,1],out_pow.9_tau.5,out_pow.8_tau.5[1,3])),
             cbind(min.power=.8,tau=0.2, out_pow.8_tau.8),
             cbind(min.power=.9,tau=0.2, cbind(out_pow.8_tau.8[1,1],out_pow.9_tau.8,out_pow.8_tau.8[1,3]))
)
out <- as.data.frame(out)
out2 <- out

out2$tau <- 1-out2$tau
out2$tau <- as.factor(out2$tau)

df <- out2[,c('min.power','tau', 'power_ss')]
rownames(df) <- NULL
df <- df[!duplicated(df),]

df$min.power <- as.factor(df$min.power)

ggplot(out2, aes(x=min.fi, y=fi_ss, shape=tau))+
  geom_point()+
  xlab(TeX("$\\varphi'$, minimum tolerable fragility index"))+
  ylab('Required Sample size')+
  geom_hline(aes(yintercept=power_ss, linetype=min.power), df[1:2,])+
  #ggtitle('FAME: Sample size calculations accounting for the fragility index')+
  labs(shape = TeX("Power $1-\\tau'$ (FI)"), linetype=TeX("Power $\\pi'$ ($p$ value)"))+
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.14,.68)
  )
# save as 600x380 #564x387 #590x370
# AAAAAA_ss_plot.png

#summary(lm(fi_ss~min.fi, data=as.data.frame(out_pow.8_tau.5)))

#### get error rates of tests at both extremes ####
print('get error rates for varying phi')
get_rates <- function(phi, n) as.data.frame(t(get.rejection.rates(get.p.val=get.p.val,
                                                                  get.sample.null=get.sample.null, get.sample.alt=get.sample.alt,
                                                                  phi=phi, n=n, alpha=0.05, cl=cl, algorithm='walsh', nsim=1*10^k)))
tic()
if (chosen.pi==.8) {
  my_rates_phivary <- as.data.frame(out_pow.8_tau.5) %>% rowwise() %>%
    do(get_rates(.$min.fi, max(.$power_ss, .$fi_ss)))
} else if (chosen.pi==.9) {
  my_rates_phivary <- as.data.frame(out_pow.8_tau.5) %>% rowwise() %>%
    do(get_rates(.$min.fi, max(out_pow.9_tau.5, .$fi_ss)))
}
toc()

my_rates_phivary <- as.data.frame(my_rates_phivary)
rownames(my_rates_phivary) <- NULL
if (chosen.pi==.8) {
  my_rates_phivary <- cbind(tau=0.50, phi=out_pow.8_tau.5[,1], n=pmax(out_pow.8_tau.5[,2],out_pow.8_tau.5[,3]), my_rates_phivary)
} else if (chosen.pi==.9) {
  my_rates_phivary <- cbind(tau=0.50, phi=out_pow.8_tau.5[,1], n=pmax(out_pow.9_tau.5,out_pow.8_tau.5[,3]), my_rates_phivary)
}

#### look at what happens when varying tau ####
tic()
print('getting n for varying tau')
tau.ss <- sapply(tau.grid.long, function(tau)
                      general.fi.samplesize(min.fi=chosen.phi, min.power=chosen.pi, alpha=0.05, tau=tau,
                              get.p.val=get.p.val, get.replacements=get.replacements, get.sample=get.sample.alt,
                              gamma=.3, nsim=nsim_high, eps=eps, niters=niters,
                              verbose=FALSE, cl=cl, algorithm='walsh'))
tau.calc.df <- data.frame('tau'=tau.grid.long, t(tau.ss))

print('getting error rates for varying tau')
my_rates_tauvary <- tau.calc.df %>% rowwise() %>% do(get_rates(chosen.phi, max(.$power_ss, .$fi_ss)))
toc()

my_rates_tauvary <- as.data.frame(my_rates_tauvary)
rownames(my_rates_tauvary) <- NULL
my_rates_tauvary <- cbind(tau=tau.calc.df[,1], phi=chosen.phi, n=pmax(tau.calc.df[,2],tau.calc.df[,3]), my_rates_tauvary)

#### getting traditional ss param table ####
print('get traditional param table')
if (chosen.pi==.8) {
  inp <- pmax(out_pow.8_tau.5[,2], out_pow.8_tau.5[,3])
} else if (chosen.pi==.9) {
  inp <- pmax(out_pow.9_tau.5[1], out_pow.8_tau.5[,3])
}

trad.out <- sapply(inp, FUN=function(n) get.traditional.ss.params(0.05, chosen.pi, my.delta,
                                                               n, get.sample.alt.f, get.p.val=get.p.val, nsim=1*10^k
                   )
)
trad.out <- t(trad.out)
rownames(trad.out) <- NULL

trad.out <- cbind(tau=chosen.tau, phi=phi.grid, n=pmax(out_pow.8_tau.5[,2], out_pow.8_tau.5[,3]), trad.out)
