# takes PMwG sampled object and calculates DIC.
# subj_ll is a 2D array of subject-level log-likelihoods, with structure [subjects,iterations].

pmwg_DIC <- function(sampled){
  nsubj <- length(unique(sampled$data$subject))
  subject.names <- as.character(unique(sampled$data$subject))
  keep <- x$samples$stage == "sample"

  # the mean likelihood of the overall (sampled-stage) model, separately for each subject
  mean.like <- apply(sampled$samples$subj_ll[,keep],1,mean)

  # the mean of each parameter across iterations. Keep dimensions for parameters and subjects
  mean.params <- t(apply(sampled$samples$alpha[,,keep],1:2,mean))

  # name mean.params here so it can be used by the log_like function
  dimnames(mean.params) <- list(subject.names, sampled$par_names)

  # log-likelihood for each subject using their mean parameter vector
  mean.params.like <- numeric(nrow(mean.params))
  names(mean.params.like) <- subject.names

  for(j in subject.names) {
    mean.params.like[j] <- sampled$ll_func(mean.params[j,], data=sampled$data[sampled$data$subject==j,], sample=FALSE)
  }

  # Effective number of parameters
  pD <- sum(-2*mean.like + 2*mean.params.like)

  # Deviance Information Criterion
  DIC <- sum(-4*mean.like + 2*mean.params.like)

  return(c("DIC"=DIC, "effective parameters"=pD))
}
