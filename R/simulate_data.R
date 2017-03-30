library(phyclust)
library(PopGenome)
library(doParallel)
dir.create("temp")

n.cores <- parallel::detectCores() #multicore:::detectCores()
n.cores <- n.cores - 1 
registerDoParallel(n.cores)
## number of simulations
numsim <- 2

###Pop params
##number of pops
Npop = 2
## sample sizes (number of alleles). For diploid data, multiply by 2.
Nsam = 93*2
## sample size of Pop1.
Pop1 = 73*2
## sample size of Pop2.
Pop2 = 20*2

ncss_table=c()
ctss_table=c()
popNlist = paste(Pop1,Pop2,sep=" ")

###Loci params
##number of nuclear loci
Nnloci = 2
## average number of base pairs for nuclear loci
Len_nloci <- 400
## mutation rate for nuclear loci
nmutrate <- 7*10^(-9)

##number of mitochondrial/chloroplast loci
Nctloci = 2
## average number of base pairs for nuclear loci
Len_ctloci <- 400
## mutation rate for nuclear loci
ctmutrate <- 7*10^(-9)

##number of SSR loci
#NSSRloci = 10
## mutation rate for nuclear loci
#SSRmutrate <- 7*10^(-6)

## years per generation
genlen <- 30
##Upper and lower bounds for effective pop size
Ne_upper = 300000
Ne_lower = 10000

theta_nloci = c()
theta_ctloci = c()
#theta_SSRloci = c()
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
    theta_nloci <<- cbind(runif(numsim, theta_nloci_lower, theta_nloci_upper))
    theta_ctloci <<- cbind(runif(numsim, theta_ctloci_lower, theta_ctloci_upper))
    #    theta_SSRloci <<- cbind(runif(numsim, theta_SSRloci_lower, theta_SSRloci_upper))
  }
  else {
    Ne = cbind(runif(numsim, Ne_lower, Ne_upper))
    theta_nloci <<- Ne*4*nmutrate*Len_nloci
    theta_ctloci <<- Ne*ctmutrate*Len_ctloci
    #    theta_SSRloci <<- Ne*4*SSRmutrate
  }
}

##number of models
Nmodels=3

nc_modelstr=sprintf("-t tbs -I %d %s",Npop, popNlist)

ct_modelstr=sprintf("-t tbs -I %d %s",Npop,popNlist)

seq_model1 = "-ej tbs 1 2"

seq_model2 = "-n 1 tbs -n 2 tbs -ej tbs 1 2"

seq_model3 = "-m 1 2 tbs -m 2 1 tbs -ej tbs 1 2"

seq_model_opt = list(seq_model1, seq_model2, seq_model3)

#SSR_model1 = "-ej %f 1 2"
#
#SSR_model2 = "-n 1 %f -n 2 %f -ej %f 1 2"
#
#SSR_model3 = "-m 1 2 %f -m 2 1 %f -ej %f 1 2"
#
#SSR_model_opt = list(SSR_model1, SSR_model2, SSR_model3)


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
    prior <- cbind(runif(numsim,  prior_list[[model_number]][[bound]], prior_list[[model_number]][[(bound+1)]]))
    parameter_list<-cbind(parameter_list, prior)
  }
  return(parameter_list)
}

#' Simulates genetic datasets using Hudson's MS and calculates summary statistics
#' from the generated data
#'
#' Uses the priors informed to simulate genetic datasets for each model 
#' using Hudson's MS, the calculates summary statistics with the 
#' PopGenome package.
#' @param A Boolean informing the genomic partition of the markers 
#' (nuclear and/or cytoplasmatic)
#' @return A table with the calculated summary statistics
#' @export
simulate_seqs <- function(seq_genome){
  if (nDNA == T){
    write(ms(nsam=Nsam, nreps=1, opts= paste(nc_modelstr, seq_model_opt[[model]]), tbs.matrix=rbind2(nc_tbs_values[rep,])),file=file.path("temp",sprintf("ncsim_%d_%09.f.txt", model, rep)))
    nc_ms<-readMS(file.path("temp",sprintf("ncsim_%d_%09.f.txt", model, rep)))
    nc_ms<-set.populations(nc_ms,list(c(1:Pop1),c((Pop1+1):(Pop1+Pop2))))
    nc_ms<-F_ST.stats(nc_ms)
    nc_ms<-diversity.stats(nc_ms)
    nc_ms<-neutrality.stats(nc_ms)
    nss<-c(nc_ms@n.segregating.sites,nc_ms@nuc.diversity.within, nc_ms@nuc.diversity.between, nc_ms@Tajima.D, nc_ms@theta_Watterson, nc_ms@Hudson.G_ST,nc_ms@hap.F_ST.pairwise)
    ncss_table<<-rbind(ncss_table, nss)
  }
  if (ctDNA == T){
    write(ms(nsam=Nsam, nreps=1, opts= paste(ct_modelstr, seq_model_opt[[model]]), tbs.matrix=rbind2(ct_tbs_values[rep,])),file=file.path("temp",sprintf("ctsim_%d_%09.f.txt", model, rep)))
    ct_ms<-readMS(file.path("temp",sprintf("ctsim_%d_%09.f.txt", model, rep)))
    ct_ms<-set.populations(ct_ms,list(c(1:Pop1),c((Pop1+1):(Pop1+Pop2))))
    ct_ms<-F_ST.stats(ct_ms)
    ct_ms<-diversity.stats(ct_ms)
    ct_ms<-neutrality.stats(ct_ms)
    css<-c(ct_ms@n.segregating.sites,ct_ms@nuc.diversity.within, ct_ms@nuc.diversity.between, ct_ms@Tajima.D, ct_ms@theta_Watterson, ct_ms@Hudson.G_ST,ct_ms@hap.F_ST.pairwise)
    ctss_table<<-rbind(ctss_table, css)
  }
}

nDNA=T
ctDNA=T
for(model in 1:Nmodels) {
#foreach(model=1:Nmodels, .packages='phyclust') %dopar%{
  tbs_values<-sample_priors(model)
  sample_theta(coalescent_theta)
  ct_tbs_values<-cbind(theta_ctloci,tbs_values)
  nc_tbs_values<-cbind(theta_nloci,tbs_values)
  #foreach (rep=1:numsim, .packages=c('phyclust','PopGenome')) %dopar%{
  for (rep in 1:numsim){
    for (nloc in 1:Nnloci){
      simulate_seqs(nDNA)
    }
    for (ctloc in 1:Nctloci){
      simulate_seqs(ctDNA)
    }
  }
  
}

##to runs microsat data (TO DO
#  for(simulation in 1:numsim){
#  simd_data <- sim_microsats(theta = 1,
#                             n_ind = c(73,20),
#                             n_loc = 10,
#                             n_pop = 2,
#                             ms_options = "-I 2 146 40 -n 1 tbs -n 2 tbs -ej tbs 1 2",
#                             tbs_matrix = cbind(0.5,0.5,2.3))
#  }
#unlink("temp",recursive=TRUE)

#nc_ms.txt<-readms(nc_ms)
#lapply(tbs_values[simulation,], function(tbs){seq_model_opt[[1]]})

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
