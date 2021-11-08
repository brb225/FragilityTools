library(ggplot2)
set.seed(1234567890)

#### get normal distribution plots ####
mu.true <- 1
y <- rnorm(150, mean = mu.true)
qplot(y, bins=13,xlim=c(-3,3)+mu.true, main='Histogram of the data')

x.grid <- seq(-3+mu.true,3+mu.true,by=.01)
y.grid <- sapply(x.grid, function(x) dnorm(x, mean=mu.true))
dd <- data.frame(x=x.grid, y=y.grid)
get_plt <- function(l,u) {
  conff <- round(pnorm(u, mean=mu.true)-pnorm(l, mean=mu.true),3)
  print(paste0(c(conff, c(l,u)), sep='; '))
  qplot(x,y,data=dd,geom="line",ylab='Density',xlab='Outcome')+
  geom_ribbon(data=subset(dd,x>l & x<u),
              aes(ymax=y,ymin=0),
              fill="grey10",colour=NA,alpha=0.5)+
  #annotate('text', x=2.5,y=.3,label=conff)+
  geom_segment(aes(x=obs, xend=obs, y=0, yend=dnorm(obs, mean=mu.true)), col='black', size=.7)+
  theme_bw()
}

obs <- -1.2
#l<-obs;u<-obs+.05;plt1 <- get_plt(l,u)
l<-obs;u<-0;plt2 <- get_plt(l,u)
l<-obs;u<-1.5;plt3 <- get_plt(l,u)
l<-obs;u<-abs(obs-mu.true)+mu.true;plt4 <- get_plt(l,u)
l<--1.7;u<-abs(l-mu.true)+mu.true;plt5 <- get_plt(l,u)
#l<--4;u<-abs(l-mu.true)+mu.true; plt6 <- get_plt(l,u)
plt <- gridExtra::grid.arrange(#plt1,
                        plt2, plt3, plt4, plt5,
                        #plt6,
                        nrow=2)

ggsave('normal_hdr_grid.eps', plt, device=cairo_ps, fallback_resolution = 300, width = 65, height = 60, units = "mm", dpi=300, scale=2) #width=75




#### get poisson distribution plots ####
set.seed(1234567890)
mu.true <- 9
y <- rpois(150, mu.true)
qplot(y, bins=13, main='Histogram of the data')

x.grid <- 0:20
y.grid <- sapply(x.grid, function(x) dpois(x, mu.true))
dd <- data.frame(x=x.grid, y=y.grid)
get_plt_pois <- function(l,u) {
  n <- nrow(subset(dd, x>=l & x<=u))

  conff <- round(dpois(u, mu.true)-dpois(l, mu.true),3)
  print(paste0(c(conff, c(l,u)), sep='; '))
  qplot(x,y,data=dd,geom=c("point"),ylab='density',xlab='outcome')+
    geom_ribbon(data=rbind(subset(dd, x>=l & x<=u),
                           subset(with(dd, data.frame(x=x+1, y)), x>=l+1 & x<=u+1))[
                             rep(c(0,n), n)+rep(1:n, each=2),
                           ],
                aes(ymax=y,ymin=0),
                fill="grey10",colour=NA,alpha=0.5)+
    #annotate('text', x=2.5,y=.3,label=conff)+
    geom_segment(aes(x=obs, xend=obs, y=0, yend=dpois(obs, mu.true)), col='black', size=1)+
    geom_segment(data=dd, aes(x=x, xend=x+1, y=dpois(x, mu.true), yend=dpois(x, mu.true)), col='black', size=.7)+
    theme_bw()
}

obs <- 13
u<-obs;l<-12;plt2 <- get_plt_pois(l,u)
u<-obs;l<-8;plt3 <- get_plt_pois(l,u)
u<-obs;l<-4;plt4 <- get_plt_pois(l,u)
u<-16;l<-3;plt5 <- get_plt_pois(l,u)
plt_pois <- gridExtra::grid.arrange(plt2, plt3, plt4, plt5, nrow=2)

ggsave('poisson_hdr_grid.eps', plt_pois, device=cairo_ps, fallback_resolution = 300, width = 65, height = 60, units = "mm", dpi=300, scale=2) #width=75
