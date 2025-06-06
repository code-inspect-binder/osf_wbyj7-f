# hierarchical LBA fit with PMwG to consumer SAT Experiment 1 phase 1 and 2 data
rm(list=ls())

# load required packages
require(rtdists)
require(tidyverse)
require(mvtnorm)
require(pmwg)

# load source functions and data
# generates four data sets: FF, FS, SF, SS
source("prep-Exp1-data.r")
source("pmwg-DIC.r")

# estimation settings
n.cores <- 10
n.posterior <- 25
burn.sets <- c(n=200, particles=100)
adapt.sets <- c(n=5000, particles=100)
sample.sets <- c(n=5000, particles=100)
epsilon <- .4


### Name for current model
model.par <- "Null"

## Name the parameters:
par_names <- c("A", "s", "utility", "threshold", "t0")

## Some good starts for a chain.
start_points <- list(
  mu = log(c(2, .5, .9, 2, .5)),
  sig2 = diag(rep(.01, length(par_names)))
)

## Some priors.
priors <- list(
  theta_mu = rep(0, length(par_names)),
  theta_sig = diag(rep(1, length(par_names)))
)


# log density for the data (RT & response) under subject-level parameters
ll <- function(x, data, sample=FALSE) {
  names(x) <- par_names
  # probit transformation to place u1/u2 on raw scale [0,1]
  up1 <- pnorm(x["utility"])
  up2 <- pnorm(x["utility"])
  # log transform for all other parameters
  x[!(par_names %in% c("utility"))] <- exp(x[!(par_names %in% c("utility"))])

  ps <- diff(range(data$RT))
  n <- nrow(data)  # number of trials
  nds <- 1:3       # number of drift rates
  phase1 <- data$phase == "1"

  # vectorise model parameters
  thresholds <- rep(x["threshold"], n)
  thresholds <- thresholds + x["A"]
  t0s <- rep(x["t0"], n)

  # transform utility parameters to drift rates for the three options
  # drift rates = exp(-utility(option))/sum(exp(-utility(option)))
  store.sums <- as.matrix(data[,c("S1","S2","S3")])
  drifts <- rbind(t(apply(store.sums[phase1,], 1,function(z) { tmp=exp(-z^up1) ; tmp/sum(tmp) })),
                  t(apply(store.sums[!phase1,],1,function(z) { tmp=exp(-z^up2) ; tmp/sum(tmp) })))

  if(sample == FALSE) {
    drifts <- as.list(data.frame(drifts))
    # call LBA likelihood
    like <- dLBA(rt=data$RT,response=data$R,
                 A=x["A"], b=thresholds, t0=t0s, mean_v=drifts, sd_v=x["s"],
                 dist="norm",silent=TRUE)
    # contaminant mixture process
    like <- (1-p.contaminant)*like + (p.contaminant/max(nds))*(1/ps)
    # Return sum log likelihood. Include protection against log(0) problems
    return(sum(log(pmax(like,1e-10))))

  } else {
    # data is subjects data - modified with model predictions and returned
    #    remove responses/RTs
    data$R <- data$correct <- data$RT <- NA
    for(i in 1:n) {
      tmp <- rLBA(n=1, A=x["A"], b=thresholds[i], t0=t0s[i], mean_v=drifts[i,], sd_v=x["s"],
                   dist="norm",silent=TRUE)
      data$RT[i]=tmp$rt
      data$R[i]=tmp$response
    }
    data$correct <- as.numeric(apply(data[,c("S1","S2","S3")], 1, which.min) == data$R)
    return(data)
  }
}



#####################################################################
# estimate model independently for each condition in the experiment #
#####################################################################

for(condition in names(all.data)) {
  cat("\n\n\n\nEstimating model for: ", condition, "\n\n")
  fnam <- paste0("Exp1-LBA-", model.par, "-", condition, ".RData")

  ## Initialise things.
  sampler <- pmwgs(
    data = all.data[[condition]],
    pars = par_names,
    prior = priors,
    ll_func = ll
  )

  ## Initialse more stuff, with starts above
  sampled <- init(sampler, theta_mu=start_points$mu, theta_sig=start_points$sig2)
  save.image(paste0("modelFits/",fnam))

  # Run some burn-in samples
  sampled <- run_stage(sampled, stage="burn", iter=burn.sets[1], particles=burn.sets[2],
                       epsilon=epsilons, n_cores=n.cores)
  save.image(paste0("modelFits/",fnam))

  # adaptation phase
  sampled <- run_stage(sampled, stage="adapt", iter=adapt.sets[1], particles=adapt.sets[2],
                       epsilon=epsilon, n_cores=n.cores)
  save.image(paste0("modelFits/",fnam))

  # sampling phase
  sampled <- run_stage(sampled, stage="sample", iter=sample.sets[1], particles=sample.sets[2],
                       epsilon=.3, n_cores=n.cores)
  save.image(paste0("modelFits/",fnam))

  # calculate DIC and then perform final save
  DIC <- pmwg_DIC(sampled)
  save.image(paste0("modelFits/",fnam))
}
