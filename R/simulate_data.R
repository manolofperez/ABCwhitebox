library(phyclust)
library(PopGenome)
library(foreach)
dir.create("step1")

dir.create("temp")

n.cores <- detectCores() #multicore:::detectCores()
n.cores <- n.cores - 7 
cl<-makeCluster(n.cores)
registerDoParallel(cl)
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
simulate_seqs <- function(nDNA, cpDNA){
  if (nDNA == T){
    for (nloc in 1:Nnloci){
      write(ms(nsam=Nsam, nreps=1, opts= paste(nc_modelstr, seq_model_opt[[model]]), tbs.matrix=rbind2(nc_tbs_values[rept,])),file=file.path("temp","ncsim.txt"))
      nc_ms<-readMS(file.path("temp","ncsim.txt"))
      nc_ms<-set.populations(nc_ms,list(c(1:Pop1),c((Pop1+1):(Pop1+Pop2))))
      nc_ms<-F_ST.stats(nc_ms)
      nc_ms<-diversity.stats(nc_ms)
      nc_ms<-neutrality.stats(nc_ms)
      nss<-c(nc_ms@n.segregating.sites,nc_ms@nuc.diversity.within, nc_ms@nuc.diversity.between, nc_ms@Tajima.D, nc_ms@theta_Watterson, nc_ms@Hudson.G_ST,nc_ms@hap.F_ST.pairwise)
      ncss_table<<-rbind(ncss_table, nss)
    }
    
  }
  if (ctDNA == T){
    for (ctloc in 1:Nctloci){
      write(ms(nsam=Nsam, nreps=1, opts= paste(ct_modelstr, seq_model_opt[[model]]), tbs.matrix=rbind2(ct_tbs_values[rept,])),file=file.path("temp","ctsim.txt"))
      ct_ms<-readMS(file.path("temp","ctsim.txt"))
      ct_ms<-set.populations(ct_ms,list(c(1:Pop1),c((Pop1+1):(Pop1+Pop2))))
      ct_ms<-F_ST.stats(ct_ms)
      ct_ms<-diversity.stats(ct_ms)
      ct_ms<-neutrality.stats(ct_ms)
      css<-c(ct_ms@n.segregating.sites,ct_ms@nuc.diversity.within, ct_ms@nuc.diversity.between, ct_ms@Tajima.D, ct_ms@theta_Watterson, ct_ms@Hudson.G_ST,ct_ms@hap.F_ST.pairwise)
      ctss_table<<-rbind(ctss_table, css)
    }
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
  for (rept in 1:numsim){
      simulate_seqs(nDNA, ctDNA)
  }
}
unlink("temp",recursive=TRUE)

##to runs microsat data (TO DO)
#  for(simulation in 1:numsim){
#  simd_data <- sim_microsats(theta = 1,
#                             n_ind = c(73,20),
#                             n_loc = 10,
#                             n_pop = 2,
#                             ms_options = "-I 2 146 40 -n 1 tbs -n 2 tbs -ej tbs 1 2",
#                             tbs_matrix = cbind(0.5,0.5,2.3))
#  }
#

### calculate summary statistics mean the simulated loci
ncss_avg <- function(ncss_table){
  ncavgss_table = c()
  for (sust in 1:length(ncss_table[1,])){
    ncavgss_table<-cbind(ncavgss_table, tapply(ncss_table[,sust], factor(rep(1:(numsim*Nmodels), each=Nnloci)), mean,na.rm=T))
  }
  return(ncavgss_table)
}

ctss_avg <- function(ctss_table){
  ctavgss_table = c()
  for (sust in 1:length(ctss_table[1,])){
    ctavgss_table<-cbind(ctavgss_table, tapply(ctss_table[,sust], factor(rep(1:(numsim*Nmodels), each=Nnloci)), mean,na.rm=T))
  }
  return(ctavgss_table)
}

### join data

ncss_table = ncss_avg(ncss_table)
ctss_table = ctss_avg(ctss_table)
ss_table = cbind(ncss_table, ctss_table)

##enter the models as a vector
models=c()
for (m in 1:Nmodels){
  models<-c(models, rep(m, each=numsim))
}
#
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
