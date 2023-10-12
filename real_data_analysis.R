# data cleaning
library(survival)
first <- subset(rhDNase, !duplicated(id)) # first row for each subject
dnase <- tmerge(first, first, id=id, tstop=as.numeric(end.dt -entry.dt))

# Subjects whose fu ended during the 6 day window are the reason for
#  this next line
temp.end <- with(rhDNase, pmin(ivstop+6, end.dt-entry.dt))
dnase <- tmerge(dnase, rhDNase, id=id,
                infect=event(ivstart),
                end=  event(temp.end))
# toss out the non-at-risk intervals, and extra variables
#  3 subjects had an event on their last day of fu, infect=1 and end=1
dnase <- subset(dnase, (infect==1 | end==0), c(id:trt, fev:infect))
dnase <- subset(dnase, !duplicated(id))

my.rhDNase <- dnase[, c(3, 4, 7, 8, 9)]
my.rhDNase$time <- my.rhDNase$tstop - my.rhDNase$tstart
my.rhDNase <- my.rhDNase[, -c(3, 4)]

# data preparation
Y <- log(my.rhDNase$time)
n <- nrow(my.rhDNase)
X <- matrix(c(rep(1, n), my.rhDNase$trt, my.rhDNase$fev), nrow = n)
delta <- my.rhDNase$infect

# priors
mu_0 <- c(4.4, 0.25, 0.04)
v_0 <- 1
alpha_0 <- 501
omega_0 <- 500

start.time <- Sys.time()
rhDNase.vb.fit <- vb_aft_final(Y, X, delta, mu_0, v_0, alpha_0, omega_0, max_iteration = 10000, threshold = 0.0005)
end.time  <- Sys.time()
end.time - start.time # take 0.8793399 secs to get the fitted model

round(rhDNase.vb.fit$mu_vector, 3)
round(rhDNase.vb.fit$omega / (rhDNase.vb.fit$alpha - 1), 3)

inv.gam(rhDNase.vb.fit$alpha, rhDNase.vb.fit$omega)
round(hdi(qinvgamma, 0.95, shape=rhDNase.vb.fit$alpha, rate=rhDNase.vb.fit$omega), 3)

round(sqrt(diag(rhDNase.vb.fit$Sigma)), 3)
rhDNase.vb.fit$converged_ELBO






