devtools::load_all() # library(FragilityTools)
set.seed(12336475)

#### binary response ####
# look at usual FI calcations for 2x2 tables
mat <- matrix(nrow = 2, byrow = TRUE, c(1, 8, 7, 1))
colnames(mat) <- c("event", "nonevent")
rownames(mat) <- c("treatment", "control")
mat # different format than usual

out <- bin.fi(mat, alg = "walsh")
out$FI # add 3 events in the control group

# look at corresponding patient-level (X,Y) data
dat <- dplyr::bind_cols(mat_to_xy(mat))
dat

# compute greedy algorithm
X <- subset(dat, select = 1)
Y <- subset(dat, select = 2)
out <- bin.fi(X = X, Y = Y, alg = "greedy", verbose = TRUE)
out$FI

#### survival response ####
# look at survival data
dat <- get.survival.data(25, 6) # status=1 (event)
dat

(X <- subset(dat, select = 1))
(Y <- subset(dat, select = -1))

# compute FI for 0 time step
dat
out <- surv.fi(dat$time, dat$status, dat$group, max.step = 0, verbose = TRUE)
out$FI

# compute FI for any time step
dat
out <- surv.fi(dat$time, dat$status, dat$group, max.step = Inf, verbose = TRUE)
out$FI
