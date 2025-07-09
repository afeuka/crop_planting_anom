###Title: Fitting mixed effects model nimble
###Author: Abbey Feuka
###Date: 03042023
###Notes: only counties with pigs

library(nimble)
library(snow)
library(tidyverse)
library(tigris)
library(sf)
options(tigris_use_cache = TRUE)

# setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")
setwd("/Users/afeuka/Desktop/Crops")

source("./crop_planting_anom/Functions/clean_crop_dat.R")

nb=FALSE
temporal=TRUE

counties <- counties()
counties <- counties %>% select(GEOID,NAME)

#all counties-------------
load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

# subfolder <- "Spatial"
subfolder <- "CRPxPigs"

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)

commod_idx <-5

for(commod_idx in 1:length(commod_names_c)){
  # only counties with pigs ----------------------------
  dat_op  <- clean_crop_dat(dat_orig=dat_orig,
                            nb=nb,
                            temporal=temporal,
                            commod_name = commod_names[commod_idx], 
                            only_pigs = T)

  dat_clean_op <- dat_op$dat_clean
  covsx_op  <- dat_op$covsx
  xmat <- dat_op$xmat
  reg_count_idx <- dat_op$reg_count_idx

  counties_samp <- counties %>% filter(GEOID%in%dat_clean_op$GEOID)
  counties_mat <- st_intersects(counties_samp,remove_self=T)
  counties_car <- as.carAdjacency(counties_mat)
  
  # check
  # sum(!dat_clean_op$GEOID%in%counties_samp$GEOID)
  
  # #correlation check 
  # x_cor<- cor(xmat)
  # x_cor_idx <- as.data.frame(which(abs(x_cor)>0.7 &x_cor!=1,arr.ind = T))
  # x_cor_idx[,1] <- covsx_op$cov[x_cor_idx[,1]]
  # x_cor_idx[,2] <- covsx_op$cov[x_cor_idx[,2]]
  # x_cor_idx
 
  nimbleMod <- nimbleCode({
    tau ~ dgamma(1,1)
    sd <- 1/tau
    
    lscale ~ dunif(0.01,10)
   
    for(i in 1:nbeta){
      beta[i] ~ dlaplace(0,lscale)
    }
  
    tau_s ~ dgamma(0.001, 0.001)
    s[1:ncounty] ~ dcar_normal(adj=adj[1:L], 
                               weights=weights[1:L], 
                               num=num[1:ncounty], 
                               tau=tau_s, 
                               zero_mean = 0)
    
    for(i in 1:nsamp){
      mu[i] <- inprod(x[i,1:(nbeta)], beta[1:nbeta]) + s[county_idx[i]]
      y[i] ~ dnorm(mu[i],sd=sd)
      ypred[i] ~ dnorm(mu[i],sd=sd)
    }
    
  })#end of nimblecode

  ##fit model --------------------
  # nChains <- 1
  # niter <- 10000
  # burnProp <- 0.2
  # thin <- 1
  
  ###multi chain --------------------
  nChains <- 3
  niter <- 100000
  burnProp <- 0.75
  thin <- 5
  
  mcmcPar <- function(j){
    require(nimble)
    
    modDat <- list(y = dat_clean_op$plant.anom,
                   x = xmat)
    
    nbeta<-ncol(modDat$x)
    
    ##set up mcmc -------------------
    #constants (not dat_cleana, not estimated)
    const <- list(nbeta=nbeta,
                  ncounty=length(unique(dat_clean_op$county_idx)),
                  nsamp=nrow(dat_clean_op),
                  county_idx=dat_clean_op$county_idx,
                  adj=counties_car$adj,
                  weights=as.numeric(counties_car$weights),
                  num=counties_car$num,
                  L=length(counties_car$adj)
    )

    #initial values (all estimated parameters)
    inits <- list(beta = rnorm(nbeta,0,10),
                  tau = rgamma(1,0.1,0.1),
                  s=rep(0,const$ncounty),
                  tau_s = 1)
    
    names(inits) <- c("beta","tau",
                      "s","tau_s")
    
    ##set up mcmc -------------------
    mod <- nimbleModel(code = nimbleMod,
                       data = modDat,
                       constants = const,
                       inits = inits)
    
    mcmc.conf <- configureMCMC(mod, enableWAIC=F) #default MCMC configuration
    monitors <- c("beta",
                  "tau",
                  "tau_s",
                  "s",
                  "mu",
                  "ypred",
                  "lscale"
                  )
    mcmc.conf$setMonitors(monitors)
    
    Cmcmc <- buildMCMC(mcmc.conf,resetFuncitions=T) #uncompiled MCMC
    Cmod <- compileNimble(mod,Cmcmc) #compiled model
    
    ##fit model --------------------
    # mcmc.out <- runMCMC(Cmod$Cmcmc,
    #                     niter = 0,
    #                     thin= thin,
    #                     nburn = 0,
    #                     # setSeed = 1:nChains,
    #                     nchains = nChains)
    
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
  clusterExport(cl, list("nimbleMod","dat_clean_op","xmat",
                         "counties_car","niter","burnProp","thin"))
  
  system.time(
    parSamples<- clusterApply(cl, 1:nChains, mcmcPar)
  )
  stopCluster(cl)
  
  samp1 <- cbind(parSamples[[1]],chain_idx=1,samp_idx=1:nrow(parSamples[[1]]))
  samp2 <- cbind(parSamples[[2]],chain_idx=2,samp_idx=1:nrow(parSamples[[1]]))
  samp3 <- cbind(parSamples[[3]],chain_idx=3,samp_idx=1:nrow(parSamples[[1]]))
  samples_all <- rbind.data.frame(samp1,samp2,samp3)
  save(samples_all,file=paste0("./Model outputs/",subfolder,"/",
                               commod_names[commod_idx],"_op_spatial.RData"))
  rm(samp1,samp2,samp3,parSamples)
}
