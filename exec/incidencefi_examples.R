library(FragilityTools) #devtools::load_all()
library(tidyverse)
theme_set(theme_bw())

# example given in Walter: The fragility of trial results involves more than statistical significance alone
# example 4a in paper
mat_walter <- matrix(nrow=2, byrow=TRUE, c(5, 90, 0, 96))
rownames(mat_walter) <- c('Treatment A', 'Treatment B')
colnames(mat_walter) <- c('Died', 'Alive')
mat_walter

out1 <- bin.fi.incidence(mat_walter)
plt1 <- incidence.plot(out1)
#ggsave('tab1_fi_plot.eps', plt1, device='eps',  width = 75, height = 60, units = "mm", dpi=320)
# 457 x 316

# LIMIT-2 TRIAL
# given in walsh paper as example of
# example 4b in paper
mat_limit <- matrix(nrow=2,byrow=TRUE,c(90, 1150-90, 118, 1150-118))
rownames(mat_limit) <- c('Magnesium', 'Placebo')
colnames(mat_limit) <- c('death', 'survived')
mat_limit

out_limit <- bin.fi.incidence(mat_limit)
plt_limit <- incidence.plot(out_limit)

# table 2
# example 4c in paper
mat2 <- matrix(nrow=2, byrow=TRUE, c(2,1,1,2))#c(6,191,19,0))
rownames(mat2) <- c('Group 1', 'Group 2')
colnames(mat2) <- c('Event', 'Nonevent')

out2 <- bin.fi.incidence(mat2)
#plt2 <- incidence.plot(out2) # not meaningful since FI is always -Inf

# table 3
# example 4d in paper
mat3 <- matrix(nrow=2, byrow=TRUE, c(24, 126, 13, 67))#c(29,12,46,101))
rownames(mat3) <- c('Treatment A', 'Treatment B')
colnames(mat3) <- c('Event', 'Nonevent')

out3 <- bin.fi.incidence(mat3)
plt3 <- incidence.plot(out3)
#ggsave('tab3_fi_plot.eps', plt3, device='eps',  width = 75, height = 60, units = "mm", dpi=320)

# table 4
# example 5a in paper
mat4 <- matrix(nrow=2, byrow=TRUE, c(75,75,5,75)) #c(2,73,47,217))
rownames(mat4) <- c('Treatment A', 'Treatment B')
colnames(mat4) <- c('Event', 'Nonevent')

out4 <- bin.fi.incidence(mat4)
plt4 <- incidence.plot(out4)
#ggsave('tab4_fi_plot.eps', plt4, device='eps',  width = 75, height = 60, units = "mm", dpi=320)

# table 6
# example 5b in the paper
mat6 <- matrix(nrow=2, byrow=TRUE, c(10, 17, 27, 65))#c(38,226,5,23))
rownames(mat6) <- c('Treatment A', 'Treatment B')
colnames(mat6) <- c('Event', 'Nonevent')

out6 <- bin.fi.incidence(mat6)
plt6 <- incidence.plot(out6)
#ggsave('tab6_fi_plot.eps', plt6, device='eps',  width = 75, height = 60, units = "mm", dpi=320)

# table 7
# "table 6" in the paper
mat7 <- matrix(byrow=TRUE, nrow=2, c(23,87,44,46))
rownames(mat6) <- c('Treatment A', 'Treatment B')
colnames(mat6) <- c('Event', 'Nonevent')

out7 <- bin.fi.incidence(mat7)
plt7 <- incidence.plot(out7)
ggsave('tab7_fi_plot.eps', plt7, device='eps',  width = 75, height = 60, units = "mm", dpi=320)

bin.fi(mat7, alg='original')$FI
