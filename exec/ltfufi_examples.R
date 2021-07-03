library(FragilityTools) #devtools::load_all() # load FragilityTools package
library(ggplot2)
set.seed(1234567890)

#### get beta plot ####
p.grid <- seq(0, 1, by=.0001)
s <- 2000
po <- .131 # gopcabe off pump

mult.grid <- c(1.1, 1.3, 1.5, 2)
d.out <- c()
s.out <- c()
for (j in 1:length(mult.grid)) {
  s <- get_beta_parameters(po, mult = mult.grid[j])
  s.out <- c(s.out, round(s, 1))
  d.out <- cbind(d.out, sapply(p.grid, function(p) dbeta(p, s*po+1, s*(1-po)+1)))
}
colnames(d.out) <- s.out
out <- cbind(p.grid, d.out)
out <- reshape2::melt(as.data.frame(out), id='p.grid')
beta.plt <- ggplot(data=out, aes(x=p.grid, y=value, col=variable))+
  geom_line()+
  xlim(c(0,.4))+
  labs(x='LTFU incidence', y='Density', color='s')+
  theme(legend.position=c(.75,.7))#85,.7))
beta.plt

#### look at a bunch of LTFU-aware fragility indices ####
mat_gopcabe <- matrix(nrow=2,byrow=TRUE,c(154, 1179-154, 167, 1191-167))
ltfu_gopcabe <- c(4+8, 5+16)
out_gopcabe <- ltfu.fi(mat_gopcabe, ltfu_gopcabe, function(m)fisher.test(m)$p.value,
        sampling.control = list(ndraws=5000000, mult=1.3),
        xlab='Off-pump LTFU Incidence', ylab='On-pump LTFU Incidence')

mat_excel <- matrix(byrow=TRUE, nrow=2, c(203,681,176,686))
ltfu_excel <- c(64,95)
out_excel <- ltfu.fi(mat_excel, ltfu_excel, function(m)fisher.test(m)$p.value,
        sampling.control = list(ndraws=5000000, mult=1.3),
        fig.size=.001,#.3,
        xlab='PCI LTFU Incidence', ylab='CABG LTFU Incidence')

mat_peterson <- matrix(byrow=TRUE, nrow=2, c(32, 101-32, 18, 91-18))
ltfu_peterson <- c(124-101, 124-91)
out_peterson <- ltfu.fi(mat_peterson, ltfu_peterson, function(m)fisher.test(m)$p.value,
        sampling.control = list(ndraws=5000000, mult=1.3),
        gradientn.scale = .95,
        xlab='Low MAP LTFU Incidence', ylab='High MAP LTFU Incidence')

#### save output ####
mydpi <- 300#400
myscale1 <- 1.2#1.2
myscale2 <- 2.1#2
ggsave('beta_curves.eps', beta.plt, device='eps', width = 65, height = 60, units = "mm", dpi=mydpi, scale=myscale1) #width=75
ggsave('gopcabe_ltfuplot.eps', out_gopcabe$plot, device='eps', width = 70, height = 60, units = "mm", dpi=mydpi, scale=myscale2)
ggsave('excel_ltfuplot.eps', out_excel$plot, device='eps', width = 70, height = 60, units = "mm", dpi=mydpi, scale=myscale2)
ggsave('peterson_ltfuplot.eps', out_peterson$plot, device='eps', width = 70, height = 60, units = "mm", dpi=mydpi, scale=myscale2)
