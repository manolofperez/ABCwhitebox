## number of simulations
numsim <- 500

###Pop params
##number of pops
Npop = 2
## sample sizes (number of alleles). For diploid data, multiply by 2.
Nsam = 93*2
## sample size of Pop1.
Pop1 = 73*2
## sample size of Pop2.
Pop2 = 20*2

###Loci params
##number of nuclear loci
Nnloci = 10
## average number of base pairs for nuclear loci
Len_nloci <- 400
## mutation rate for nuclear loci
nmutrate <- 7*10^(-9)

##number of mitochondrial/chloroplast loci
Nctloci = 10
## average number of base pairs for nuclear loci
Len_ctloci <- 400
## mutation rate for nuclear loci
ctmutrate <- 7*10^(-9)

##number of SSR loci
NSSRloci = 10
## mutation rate for nuclear loci
SSRmutrate <- 7*10^(-6)

## years per generation
genlen <- 30
##Upper and lower bounds for effective pop size
Ne_upper = 100
Ne_lower = 30000

coalescent_theta = F

#' Sample the theta value for MS's simulations
#'
#' Takes the theta values informed by the user or calculates them from Ne
#' and mutation rates.
#' @param coalescent_theta A boolean, if TRUE the theta is informed by the user,
#' if FALSE should be calculated from Ne and mutation rates.
#' @return Theta values for each partition (nuclear, citoplasmatic and SSR).
#' @export
sample_theta <- function(coalescent_theta){
  if (coalescent_theta == T){
    theta_nloci = runif(1, theta_nloci_lower, theta_nloci_upper)
    theta_ctloci = runif(1, theta_ctloci_lower, theta_ctloci_upper)
    theta_SSRloci = runif(1, theta_SSRloci_lower, theta_SSRloci_upper)
  }
  else {
    Ne = runif(1, Ne_lower, Ne_upper)
    theta_nloci = Ne*4*nmutrate*Len_nloci
    theta_ctloci = Ne*ctmutrate*Len_ctloci
    theta_SSRloci = tNe*4*SSRmutrate
  }
}

##number of models
Nmodels=3

model1 = "-ej %f 1 2"

model2 = "-n 1 %f -n 2 %f -ej %f 1 2"

model3 = "-m 1 2 %f -m 2 1 %f -ej %f 1 2"

model_opt = list(model1, model2, model3)

###TO DO: use a conditional to sample times in years or in Coalescent_times prior
#coalDivTime <- DivTime/(genlen*4*Ne)

### Define parameters
prior1 = list(0.1, 4)
prior2 = list(0.001, 1, 0.001, 1, 0.1, 4)
prior3 = list(0.1, 5, 0.1, 5, 0.1, 4)

prior_list = list(prior1, prior2, prior3)

#' Sample priors for MS's simulations
#'
#' Sample priors from a uniform distribution using the lower and upper bounds
#' informed by the user.
#' @param prior_list A list containing lower and upper bounds for all the priors required
#' by the model
#' @return sampled priors from uniform distributions
#' @export
sample_priors <- function(model_number){
  parameter_list = c()
  for (bound in seq(1, length(prior_list[[model_number]]), 2)){
    prior <- runif(1,  prior_list[[model_number]][[bound]], prior_list[[model_number]][[(bound+1)]])
    parameter_list <- c(parameter_list, prior)
    return(parameter_list)
  }
}

for(model in Nmodels){
#  for (simulation in 1:numsim) {
    sample_priors(model)
#  }
}
#    ## use the DivTime in years to calculte divergence time in colaescent units (required by ms)
#
#    ## ms's command
#    system(sprintf("./ms %d %d -t %f -I 2 %f %f -ej %f 1 2 | ./sample_stats >> const.txt", Nsam, numloc, theta, Pop1, Pop2, coalDivTime))
#
#    ## save parameter values
#    parameters <- rbind(parameters, data.frame(Ne, DivTime, RatioPop1, RatioPop2, Migration12, Migration21))
#  }
#}
#
### bottleneck in both pops model
#for (i in 1:numsim) {
#
#  ### Define parameters
#  Ne <- runif(1, 100, 30000)
#  theta <- Ne*4*mutrate*L
#  DivTime <- runif(1, 300000, 4000000)
#  coalDivTime <- DivTime/(genlen*4*Ne)
#  ##Bottleneck intensity in each population, measured as a rate of the Ne in the original population. Prior sampled from an uniform distribution from 0.01 to 1.
#  RatioPop1 <- runif(1, 0.01, 1)
#  RatioPop2 <- runif(1, 0.01, 1)
#
#  ## ms's command
#  system(sprintf("./ms %d %d -t %f -I 2 %f %f -n 1 %f -n 2 %f -ej %f 1 2 | ./sample_stats >> bott.txt", Nsam, numloc, theta, Pop1, Pop2, RatioPop1, RatioPop2, coalDivTime))
#
#  ## save parameter values
#  parameters <- rbind(parameters, data.frame(Ne, DivTime, RatioPop1, RatioPop2, Migration12, Migration21))
#}
#
### Bidirectional migration model
#for (i in 1:numsim) {
#
#  ### Define parameters
#  Ne <- runif(1, 100, 30000)
#  theta <- Ne*4*mutrate*L
#  DivTime <- runif(1, 300000, 4000000)
#  coalDivTime <- DivTime/(genlen*4*Ne)
#  ## migration prior sampled sampled from an uniform distribution from 0.1 to 5 migrants per generation.
#  Migration12 <- runif(1, 0.1, 5)
#  Migration21 <- runif(1, 0.1, 5)
#
#  ## ms's command
#  system(sprintf("./ms %d %d -t %f -I 2 %f %f -m 1 2 %f -m 2 1 %f -ej %f 1 2 | ./sample_stats >> migr.txt", Nsam, numloc, theta, Pop1, Pop2, Migration12, Migration21, coalDivTime))
#
#  ## save parameter values
#  parameters <- rbind(parameters, data.frame(Ne, DivTime, RatioPop1, RatioPop2, Migration12, Migration21))
#}
#
###sample_stats exports data in the following format: pi:	0.404999	ss:	1	D:	1.381008	thetaH:	0.157164	H:	0.247835
### calculate summary statistics mean and variance over the 25 loci
#const <- read.table("const.txt")
#const <- data.frame(pi.m=tapply(const[,2]/L, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                    ss.m=tapply(const[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                    D.m=tapply(const[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                    T.m=tapply(const[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                    TH.m=tapply(const[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                    pi.v=tapply(const[,2]/L, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
#                    ss.v=tapply(const[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
#                    D.v=tapply(const[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T))
#bott <- read.table("bott.txt")
#bott <- data.frame(pi.m=tapply(bott[,2]/L, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   ss.m=tapply(bott[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   D.m=tapply(bott[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   T.m=tapply(bott[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   TH.m=tapply(bott[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   pi.v=tapply(bott[,2]/L, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
#                   ss.v=tapply(bott[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
#                   D.v=tapply(bott[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T))
#migr <- read.table("migr.txt")
#migr <- data.frame(pi.m=tapply(migr[,2]/L, factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   ss.m=tapply(migr[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   D.m=tapply(migr[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   T.m=tapply(migr[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   TH.m=tapply(migr[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   pi.v=tapply(migr[,2]/L, factor(rep(1:numsim, each=numloc)), var,na.rm=T),
#                   ss.v=tapply(migr[,4], factor(rep(1:numsim, each=numloc)), var,na.rm=T),
#                   D.v=tapply(migr[,6], factor(rep(1:numsim, each=numloc)), var,na.rm=T))
#
#models <- rep(c("const", "bott", "migr"), each=numsim)
#
### join data
#sust <- rbind(const, bott, migr)
#names(sust) <- c("pi.m", "ss.m", "TajD.m", "T.m", "TH.m","pi.v", "ss.v","TajD.v")
#
###optional lines to save the files in the folder, not only in the R environment
###write.table(models, file="models.txt", quote=F, row.names=F, col.names=F)
###write.table(parameters, file="parameters.txt", quote=F, row.names=F, col.names=F)
###write.table(sust, file="sust.txt", quote=F, row.names=F, col.names=F)
#
### erase temporary files
###system(sprintf("rm -rf const.txt bott.txt migr.txt"))
#