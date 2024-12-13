###Title: Fitting mixed effects model nimble
###Author: Abbey Feuka
###Date: 03042023
###Notes: all counties
library(nimble)
library(snow)
library(tidyverse)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")

source("./crop_planting_anom/Functions/clean_crop_dat.R")

#all counties-------------
load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)
# commod_idx <-1
for(commod_idx in 1:length(commod_names_c)){
  dat <- clean_crop_dat(dat_orig=dat_orig,commod_name = commod_names[commod_idx], only_pigs = F)
  
  dat_clean <- dat$dat_clean
  covsx  <- dat$covsx
  xmat <- dat$xmat
  reg_count_idx <- dat$reg_count_idx

  if("take.hog.intens.5yeartrend.sc"%in%colnames(dat_clean)){
    subfolder <- "Take Trend"
  } else if("take.hog.intens.prev.sc"%in%colnames(dat_clean)) {
    subfolder <- "Take Previous" 
  } else {
    stop("Must include take covariate. Check clean_crop_dat.R")
  }
  # #correlation check 
  # x_cor<- cor(xmat)
  # x_cor_idx <- as.data.frame(which(abs(x_cor)>0.7 &x_cor!=1,arr.ind = T))
  # x_cor_idx[,1] <- covsx$cov[x_cor_idx[,1]]
  # x_cor_idx[,2] <- covsx$cov[x_cor_idx[,2]]
  # x_cor_idx
  
  # hist(1/rgamma(10000,1,1))#sd
  ##model setup ----------------------------
  nimbleMod <- nimbleCode({
    tau ~ dgamma(1,1)
    sd <- 1/tau
    
    tau_region ~ dgamma(1,1)
    sd_region <- 1/tau_region
    
    tau_county ~ dgamma(1,1)
    sd_county <- 1/tau_county
    
    lscale ~ dunif(0.1,1000)
    
    # beta[1:nbeta] ~ dmnorm(mu_beta[1:nbeta],Sig_beta[1:nbeta,1:nbeta])
    for(i in 1:nbeta){
      beta[i] ~ dlaplace(0,lscale)
    }
    
    for(i in 1:nregion){
      beta_region[i] ~ dnorm(beta[1],sd=sd_region)
    }
    
    for(i in 1:ncounty){
      beta_county[i] ~ dnorm(beta_region[region_idx_county[i]],sd_county)
      # alpha[i] <- beta_county[i] - beta[1]
    }
    
    for(i in 1:nsamp){
      mu[i] <- beta_county[county_idx[i]] + inprod(x[i,1:(nbeta-1)], beta[2:nbeta]) 
      y[i] ~ dnorm(mu[i],sd=sd)
      ypred[i] ~ dnorm(mu[i],sd=sd)
    }
    
  })#end of nimblecode
  
  ##fit model --------------------
  ###multi chain --------------------
  nChains <- 3
  niter <- 50000
  burnProp <- 0.5
  thin <- 5
  
  mcmcPar <- function(j){
    require(nimble)
    modDat <- list(y = dat_clean$plant.anom,
                   x = xmat)
    
    nbeta<-ncol(modDat$x)+1
    
    const <- list(nbeta=nbeta,
                  ncounty=length(unique(dat_clean$county_idx)),
                  nsamp=nrow(dat_clean),
                  county_idx=dat_clean$county_idx,
                  mu_beta=rep(0,nbeta),
                  Sig_beta=diag(0.0001,nbeta)
    )
    
    inits <- list(beta = rnorm(nbeta,0,10),
                  tau = rgamma(1,0.1,0.1),
                  beta_county = rnorm(const$ncounty,0,10))
    names(inits) <- c("beta","tau","beta_county")
    
    ##set up mcmc -------------------
    modDat <- list(y = dat_clean$plant.anom,
                   x = xmat)
    
    #constants (not dat_cleana, not estimated)
    const <- list(nbeta=nbeta,
                  ncounty=length(unique(dat_clean$county_idx)),
                  nregion=length(unique(dat_clean$region_idx)),
                  nsamp=nrow(dat_clean),
                  county_idx=dat_clean$county_idx,
                  # region_idx=dat_clean$region_idx,
                  region_idx_county=reg_count_idx$reg_count_idx
                  # mu_beta=rep(0,nbeta),
                  # Sig_beta=diag(0.0001,nbeta)
    )
    
    #initial values (all estimated parameters)
    inits <- list(beta = rnorm(nbeta,0,10),
                  tau = rgamma(1,1,1),
                  beta_county = rnorm(const$ncounty,0,10),
                  tau_region=rgamma(1,1,1),
                  tau_county=rgamma(1,1,1))
    names(inits) <- c("beta","tau","beta_county","tau_region","tau_county")
    
    ##set up mcmc -------------------
    mod <- nimbleModel(code = nimbleMod,
                       data = modDat,
                       constants = const,
                       inits = inits)
    
    mcmc.conf <- configureMCMC(mod, enableWAIC=F) #default MCMC configuration
    monitors <- c("beta","tau",#"beta_county",
                  "beta_region","ypred",
                  "lscale")
    mcmc.conf$setMonitors(monitors)
    
    Cmcmc <- buildMCMC(mcmc.conf,resetFuncitions=T) #uncompiled MCMC
    Cmod <- compileNimble(mod,Cmcmc) #compiled model
    
    ##fit model --------------------
    mcmc.out <- runMCMC(Cmod$Cmcmc, 
                        niter = niter, 
                        nburn = niter*burnProp,
                        thin=thin,
                        setSeed = 1, 
                        nchains = 1, 
                        WAIC = F)
    return(mcmc.out)
  }
  
  cl <- makeCluster(nChains, "SOCK")
  clusterExport(cl, list("nimbleMod","dat_clean","xmat","reg_count_idx","niter","burnProp","thin"))
  
  system.time(
    parSamples<- clusterApply(cl, 1:nChains, mcmcPar)
  )
  stopCluster(cl)
  
  samp1 <- cbind(parSamples[[1]],chain_idx=1,samp_idx=1:nrow(parSamples[[1]]))
  samp2 <- cbind(parSamples[[2]],chain_idx=2,samp_idx=1:nrow(parSamples[[1]]))
  samp3 <- cbind(parSamples[[3]],chain_idx=3,samp_idx=1:nrow(parSamples[[1]]))
  samples_all <- rbind.data.frame(samp1,samp2,samp3)
  save(samples_all,file=paste0("./Model outputs/",subfolder,"/",
                               commod_names[commod_idx],"_region_multi.RData"))
  rm(samp1,samp2,samp3,parSamples)
}


