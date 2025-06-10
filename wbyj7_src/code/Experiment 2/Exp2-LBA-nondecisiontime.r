# hierarchical LBA fit with PMwG to consumer SAT Experiment 2 phase 1 and 2 data
rm(list=ls())

# load required packages
require(rtdists)
require(tidyverse)
require(mvtnorm)
require(pmwg)

Expt <- "Exp2"

# load source functions and data
# generates four data sets: FS, SS
source(paste0("prep-", Expt, "-data.r"))
source("pmwg-DIC.r")

# estimation settings
n.cores <- 8
n.posterior <- 25
burn.sets <- c(n=200, particles=200)
adapt.sets <- c(n=5000, particles=200)
sample.sets <- c(n=5000, particles=100)



### Name for current model
model.par <- "NonDecisionTime"

## Name the parameters:
par_names <- c(paste(c("A", "s", "utility", "threshold", "t0P1", "t0P2"), "store", sep="_"),
               paste(c("A", "s", "drift_mean", "drift_sd", "threshold", "t0"), "dots", sep="_"))

## Some good starts for a chain.
start_points <- list(
  mu = log(c(2, .5, .9, 2, .5, .5,
             2, .5, 50, 5, 2, .5)),
  sig2 = diag(rep(.01, length(par_names)))
)

## Some priors.
theta_mu = rep(0, length(par_names))
# centre drift_mean to 4 (log(50)) for unbiased responding...
theta_mu[par_names == "drift_mean_dots"] <- log(50)

priors <- list(
  theta_mu_mean= theta_mu,
  theta_mu_var= diag(rep(1, length(par_names)))
)


# log density for the data (RT & response) under subject-level parameters
ll <- function(x, data, sample=FALSE) {
  names(x) <- par_names
  # probit transformation to place u1/u2 on raw scale [0,1]
  up1 <- pnorm(x["utility_store"])
  up2 <- pnorm(x["utility_store"])
  # log transform for all other parameters
  x[!(par_names %in% c("utility_store"))] <- exp(x[!(par_names %in% c("utility_store"))])

  # independent estimation for store and dots tasks
  store.data <- subset(data, task == "store")
  dots.data <- subset(data, task == "dots")

  ps.store <- diff(range(store.data$RT))
  ps.dots <- diff(range(dots.data$RT))
  # number of trials in the store and dots tasks
  n.store <- nrow(store.data)
  n.dots <- nrow(dots.data)
  nds <- 1:4       # number of drift rates for store task
  phase1 <- store.data$phase == "1"

  # vectorise model parameters first for the store task
  thresholds.store <- rep(x["threshold_store"], n.store)
  thresholds.store <- thresholds.store + x["A_store"]
  t0s.store <- rep(x["t0P1_store"], n.store)
  t0s.store[!phase1] <- x["t0P2_store"]

  # transform utility parameters to drift rates for the three options
  # drift rates = exp(-utility(option))/sum(exp(-utility(option)))
  store.sums <- as.matrix(store.data[,paste0("S",nds)])
  drifts.store <- rbind(t(apply(store.sums[phase1,], 1,function(z) { tmp=exp(-z^up1) ; tmp/sum(tmp) })),
                        t(apply(store.sums[!phase1,],1,function(z) { tmp=exp(-z^up2) ; tmp/sum(tmp) })))

  # parameters for the dots task
  thresholds.dots <- rep(x["threshold_dots"], n.dots)
  thresholds.dots <- thresholds.dots + x["A_dots"]
  t0s.dots <- rep(x["t0_dots"], n.dots)
  # cumulative normal link function with freely estimated mean and SD. generates drift rate for option 2
  tmp <- pnorm(dots.data$numdot, mean=x["drift_mean_dots"], sd=x["drift_sd_dots"])
  # calculate option 1 drift rate with sum-to-1 constraint on drift rates (like drifts.store)
  drifts.dots <- cbind(1-tmp, tmp)

  if(sample == FALSE) {
    drifts.store <- as.list(data.frame(drifts.store))
    drifts.dots <- as.list(data.frame(drifts.dots))

    # call LBA likelihood
    like.store <- dLBA(rt=store.data$RT, response=store.data$R,
                       A=x["A_store"], b=thresholds.store, t0=t0s.store, mean_v=drifts.store, sd_v=x["s_store"],
                       dist="norm",silent=TRUE)
    like.dots <-  dLBA(rt=dots.data$RT, response=dots.data$R,
                       A=x["A_dots"], b=thresholds.dots, t0=t0s.dots, mean_v=drifts.dots, sd_v=x["s_dots"],
                       dist="norm",silent=TRUE)
    # contaminant mixture process
    like.store <- (1-p.contaminant)*like.store + (p.contaminant/max(nds))*(1/ps.store)
    like.dots <- (1-p.contaminant)*like.dots + (p.contaminant/2)*(1/ps.dots)
    like <- c(like.store, like.dots)
    # Return sum log likelihood. Include protection against log(0) problems
    return(sum(log(pmax(like, 1e-10))))

  } else {
    # data is subjects data - modified with model predictions and returned
    #    remove responses/RTs
    store.data$R <- store.data$correct <- store.data$RT <- NA
    dots.data$R <- dots.data$correct <- dots.data$RT <- NA
    # first simulate store task
    for(i in 1:n.store) {
      tmp <- rLBA(n=1, A=x["A_store"], b=thresholds.store[i], t0=t0s.store[i],
                  mean_v=drifts.store[i,], sd_v=x["s_store"], dist="norm",silent=TRUE)
      store.data$RT[i] <- tmp$rt
      store.data$R[i] <- tmp$response
    }
    store.data$correct <- as.numeric(apply(store.data[,paste0("S",nds)], 1, which.min) == store.data$R)
    # second simulate dots task
    for(i in 1:n.dots) {
      tmp <- rLBA(n=1, A=x["A_dots"], b=thresholds.dots[i], t0=t0s.dots[i],
                  mean_v=drifts.dots[i,], sd_v=x["s_dots"], dist="norm",silent=TRUE)
      dots.data$RT[i] <- tmp$rt
      dots.data$R[i] <- tmp$response
    }
    dots.data$correct <- 0
    dots.data$correct[dots.data$numdot < 50 & dots.data$R == 1] <- 1
    dots.data$correct[dots.data$numdot > 50 & dots.data$R == 2] <- 1
    data <- bind_rows(store.data, dots.data)
    return(data)
  }
}



#####################################################################
# estimate model independently for each condition in the experiment #
#####################################################################

for(condition in names(all.data)) {
  cat("\n\n\n\nEstimating model for: ", condition, "\n\n")
  fnam <- paste0(Expt, "-LBA-", model.par, "-", condition, ".RData")

  ## Initialise things.
  sampler <- pmwgs(
    data = all.data[[condition]],
    pars = par_names,
    prior = priors,
    ll_func = ll
  )

  ## Initialse more stuff, with starts above
  sampled <- init(sampler, start_mu=start_points$mu, start_sig=start_points$sig2)
  save.image(paste0("modelFits/",fnam))

  # Run some burn-in samples
  sampled <- run_stage(sampled, stage="burn", iter=burn.sets[1], particles=burn.sets[2], n_cores=n.cores)
  save.image(paste0("modelFits/",fnam))

  # adaptation phase
  sampled <- run_stage(sampled, stage="adapt", iter=adapt.sets[1], particles=adapt.sets[2], n_cores=n.cores)
  save.image(paste0("modelFits/",fnam))

  # sampling phase
  sampled <- run_stage(sampled, stage="sample", iter=sample.sets[1], particles=sample.sets[2], n_cores=n.cores)
  save.image(paste0("modelFits/",fnam))

  # calculate DIC and then perform final save
  DIC <- pmwg_DIC(sampled)
  save.image(paste0("modelFits/",fnam))
}
