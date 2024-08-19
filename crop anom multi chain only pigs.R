###Title: Fitting mixed effects model nimble
###Author: Abbey Feuka
###Date: 03042023
###Notes: only counties with pigs
library(nimble)
library(snow)
library(tidyverse)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")

#all counties-------------
load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)

for(commod_idx in 1:length(commod_names_c)){
  dat <- dat_orig %>% 
    filter(commodity_desc==commod_names_c[commod_idx]) %>% 
    filter(!is.nan(plant.anom) & 
             !is.na(plant.anom) & 
             !is.na(GEOID) & 
             year>=2009 
    )
  # only counties with pigs ----------------------------
  dat_op <- dat %>% filter(year<2023) %>% 
    filter(!is.nan(plant.anom) & !is.na(plant.anom) & 
             !is.na(GEOID) & year>=2009 &
             (pig.last.year==1 | pig.in.year==1 )#& long>=-112
    )%>% 
    distinct() 
  
  pig_cov <- "take.hog.intens" 
  
  covs_op <- data.frame(
    cov=c("plant.anom",
          "GEOID",
          "division_grp",
          # "pig.last.year.sc",
          "temp5trend.sc",
          "precip5trend.sc",
          "reg.roi5trend.sc",
          "plant.anom.nb.sc",
          "plant.anom.prev.sc",
          # "take5trend.sc",
          "take.hog.intens.sc",
          "prop.nfsp.sc",
          "crp.prop.sc"
    ),
    name=c("Planting anomaly",
           "County",
           "Ecoregion",
           # "Pig presence previous year",
           "Temperature 5 yr trend",
           "Precipitation 5 yr trend",
           "ROI 5 yr trend",
           "Neighboring planting anomaly",
           "Previous year's planting anomaly",
           # "Take 5 yr trend",
           "Take per hog intensity",
           "Prop. of county with pigs",
           "Prop. CRP land"
    ))
  pig_cov_name <- covs_op$name[which(covs_op$cov==paste0(pig_cov,".sc"))]
  
  covsx_op <- covs_op[-(1:3),]
  
  dat_clean_op <- na.omit(dat_op[,covs_op$cov])
  dat_clean_op$county_idx <- as.numeric(factor(dat_clean_op$GEOID))
  dat_clean_op$region_idx <- as.numeric(factor(dat_clean_op$division_grp))
  
  reg_count_idx <- dat_clean_op %>% group_by(county_idx) %>% 
    summarise(reg_count_idx=unique(region_idx)) %>% 
    arrange(county_idx)
  
  xmat <- dat_clean_op[,covsx_op$cov]
  xmat <- sapply(xmat,as.numeric)
  
  if(pig_cov=="ever.pigs" | pig_cov=="nfsp.level" | pig_cov=="hog.intensity"){
    tidx <- which(covsx_op$cov=="take5trend.sc")
    xmat <- xmat[,-tidx]
    covsx_op <- covsx_op[-tidx,]
  }
 
  #correlation check 
  x_cor<- cor(xmat)
  x_cor_idx <- as.data.frame(which(abs(x_cor)>0.7 &x_cor!=1,arr.ind = T))
  x_cor_idx[,1] <- covsx_op$cov[x_cor_idx[,1]]
  x_cor_idx[,2] <- covsx_op$cov[x_cor_idx[,2]]
  x_cor_idx
 
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
  niter <- 100000
  burnProp <- 0.5
  thin <- 5
  
  mcmcPar <- function(j){
    require(nimble)
    modDat <- list(y = dat_clean_op$plant.anom,
                   x = xmat)
    
    nbeta<-ncol(modDat$x)+1
    
    const <- list(nbeta=nbeta,
                  ncounty=length(unique(dat_clean_op$county_idx)),
                  nsamp=nrow(dat_clean_op),
                  county_idx=dat_clean_op$county_idx,
                  mu_beta=rep(0,nbeta),
                  Sig_beta=diag(0.0001,nbeta)
    )
    
    inits <- list(beta = rnorm(nbeta,0,10),
                  tau = rgamma(1,0.1,0.1),
                  beta_county = rnorm(const$ncounty,0,10))
    names(inits) <- c("beta","tau","beta_county")
    
    ##set up mcmc -------------------
    modDat <- list(y = dat_clean_op$plant.anom,
                   x = xmat)
    
    #constants (not dat_cleana, not estimated)
    const <- list(nbeta=nbeta,
                  ncounty=length(unique(dat_clean_op$county_idx)),
                  nregion=length(unique(dat_clean_op$region_idx)),
                  nsamp=nrow(dat_clean_op),
                  county_idx=dat_clean_op$county_idx,
                  # region_idx=dat_clean_op$region_idx,
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
  clusterExport(cl, list("nimbleMod","dat_clean_op","xmat","reg_count_idx","niter","burnProp","thin"))
  
  system.time(
    parSamples<- clusterApply(cl, 1:nChains, mcmcPar)
  )
  stopCluster(cl)
  
  samp1 <- cbind(parSamples[[1]],chain_idx=1,samp_idx=1:nrow(parSamples[[1]]))
  samp2 <- cbind(parSamples[[2]],chain_idx=2,samp_idx=1:nrow(parSamples[[1]]))
  samp3 <- cbind(parSamples[[3]],chain_idx=3,samp_idx=1:nrow(parSamples[[1]]))
  samples_all <- rbind.data.frame(samp1,samp2,samp3)
  save(samples_all,file=paste0("./Model outputs/",commod_names[commod_idx],"_op_region_rs_multi.RData"))
  rm(samp1,samp2,samp3,parSamples)
}


##plots --------------------
###trace plots -------------------------
plot(samples_op[,"tau"],typ="l",main="trace")

beta <- samples_op[,grepl("beta",colnames(samples_op)) & !grepl("beta_county",colnames(samples_op))]
beta <- as.data.frame(beta)

beta_cor <- cor(beta)
beta_cor_idx <- as.data.frame(which(abs(beta_cor)>0.7 &beta_cor!=1,arr.ind = T))
beta_cor_idx[,1] <- c("int",covsx_op$cov)[beta_cor_idx[,1]]
beta_cor_idx[,2] <- c("int",covsx_op$cov)[beta_cor_idx[,2]]
beta_cor_idx

beta$idx <- 1:nrow(beta)
beta_long <- beta %>%  
  pivot_longer(grep("beta",colnames(beta)),
               values_to="value",names_to="beta")
beta_long$cov <- rep(c("Intercept",covsx_op$name),nrow(beta))

ggplot(beta_long %>% filter(beta!="beta[1]"))+
  geom_line(aes(x=idx,y=value))+
  geom_hline(yintercept=0,col="blue",lty=2)+
  facet_wrap(.~cov)

# plot(samples_op[,paste0("beta[",which(c("int",covsx_op$cov)=="lat.sc"),"]")],typ="l")

##model checks ---------------------
ypred <- samples_op[,grep("ypred",colnames(samples_op))]
ypred_sum <- data.frame(mn=colMeans(ypred),
                        lci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.025)),
                        uci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.975)),
                        obs = dat_clean_op$plant.anom)

ypred_sum$resid <- ypred_sum$obs - ypred_sum$mn
ypred_sum$std_resid <- scale(ypred_sum$resid)

ypred_sum = cbind(ypred_sum,dat_clean_op)
ypred_sum <- ypred_sum %>% 
  left_join(dat %>% dplyr::select(GEOID,lat,long) %>% distinct())

###posterior predictive distribution ----------------
ggplot(ypred_sum %>% pivot_longer(cols=c("mn","obs"),
                                  values_to="value",
                                  names_to="typ"))+
  geom_density(aes(x=value,fill=typ),alpha=0.5,col="transparent")+
  scale_fill_discrete(name="",labels=c("Post pred",
                                       "Data"))

##simulated model stats --------------------
ypred_stat <- data.frame(mn=rowMeans(ypred),
                         sd=sapply(1:nrow(ypred),function(i)sqrt(var(ypred[i,]))),
                         idx=1:nrow(ypred))
ggplot()+
  geom_density(data=ypred_stat,aes(x=mn))+
  geom_vline(xintercept=mean(dat_clean_op$plant.anom),col="red",lty=2)

ggplot()+
  geom_density(data=ypred_stat,aes(x=sd))+
  geom_vline(xintercept=sqrt(var(dat_clean_op$plant.anom)),col="red",lty=2)

### standardized residuals--------------------
ggplot(ypred_sum)+
  geom_point(aes(x=mn,y=std_resid,size=long),alpha=0.7)+
  geom_hline(yintercept=0,col="blue",lty=2)

covs_outliers <- ypred_sum %>% filter(std_resid<=-5) %>% 
  pivot_longer(cols=c(grep("sc",colnames(ypred_sum)),
                      grep("lat",colnames(ypred_sum)),
                      grep("long",colnames(ypred_sum))),
               names_to="cov",values_to="value") %>% 
  group_by(cov) %>% 
  summarise(min=min(value),
            max=max(value),
            typ="outlier")

covs_op <- ypred_sum %>% 
  pivot_longer(cols=c(grep("sc",colnames(ypred_sum)),
                      grep("lat",colnames(ypred_sum)),
                      grep("long",colnames(ypred_sum))),
               names_to="cov",values_to="value") %>% 
  group_by(cov) %>% 
  summarise(min=min(value),
            max=max(value),
            typ="all")

covs_resid <- full_join(covs_outliers,covs_op) 

ggplot(covs_resid)+
  geom_errorbar(aes(x=cov,ymin=min,ymax=max,col=typ),
                position = position_dodge(width=0.5))+
  theme(axis.text.x=element_text(angle=90,hjust=1))

##fixed effects --------------------------
beta_sum <- beta_long %>% 
  group_by(cov) %>% 
  summarise(mn=mean(value),
            lci=quantile(value,probs=0.025),
            uci=quantile(value,probs=0.975))

beta_sum$signif <- ifelse(sign(beta_sum$lci)==sign(beta_sum$uci),TRUE,FALSE)

ggplot(beta_sum %>% filter(cov!="Intercept"))+
  geom_point(aes(y=cov,x=mn,alpha=signif),size=2.5)+
  geom_errorbar(aes(y=cov,xmin=lci,xmax=uci,alpha=signif),width=0,lwd=2)+
  geom_vline(xintercept = 0,col="blue",lty=2)+
  scale_alpha_manual(values=c(0.3,1),name="Significant effect")+
  ylab("Covariate")+xlab("Coefficient estimate")+
  ggtitle("Counties with pigs")+
  theme(text=element_text(size=15))
ggsave(filename = "./Model outputs/Plots/op_counties_betas.jpeg",
       device = "jpeg",
       width=8,height=8,units="in")

##random effects ---------------------------------
beta_county <- samples_op[,grep("beta_county",colnames(samples_op))]
beta_county_sum <-data.frame(mn=colMeans(beta_county),
                             lci=sapply(1:ncol(beta_county),function(i)quantile(beta_county[,i],probs=0.025)),
                             uci=sapply(1:ncol(beta_county),function(i)quantile(beta_county[,i],probs=0.975)),
                             GEOID=unique(dat_clean_op$GEOID))

ggplot(beta_county_sum)+
  geom_point(aes(y=reorder(GEOID,mn),x=mn),alpha=0.5)+
  geom_errorbar(aes(y=reorder(GEOID,mn),xmin=lci,xmax=uci),width=0,alpha=0.3,col="grey44")+
  geom_vline(xintercept = 0,col="blue",lty=2)+
  ylab("County")+xlab("Random effect estimate")+
  theme(text=element_text(size=15),
        axis.text.y=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "./Model outputs/Plots/op_counties_re.jpeg",
       device = "jpeg",
       width=8,height=8,units="in")

##covariate relationships --------------------------------
### neighboring ----------------------
anom_nb <- data.frame(sc=seq(range(dat_clean_op$plant.anom.nb.sc)[1],
                             range(dat_clean_op$plant.anom.nb.sc)[2],by=0.5))
anom_nb$bt <- anom_nb$sc *attr(dat_clean_op$plant.anom.nb.sc, 'scaled:scale') + 
  attr(dat_clean_op$plant.anom.nb.sc, 'scaled:center')

beta_int <- beta_long %>% filter(cov=="Intercept") %>% select(value)
beta_int <- unlist(c(beta_int),use.names = F)
beta_anom <- beta_long %>% filter(cov=="Neighboring planting anomaly") %>% select(value)
beta_anom <- unlist(c(beta_anom),use.names = F)
est_nb_sc <- sapply(1:length(beta_anom),function(i)beta_int[i] + beta_anom[i]*anom_nb$sc)
est_nb_sum <- data.frame(mn=rowMeans(est_nb_sc),
                         lci=sapply(1:nrow(est_nb_sc),function(i)quantile(est_nb_sc[i,],probs=0.025)),
                         uci=sapply(1:nrow(est_nb_sc),function(i)quantile(est_nb_sc[i,],probs=0.975)))
est_nb_sum <- cbind.data.frame(anom_nb,est_nb_sum)

ggplot(est_nb_sum)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Neighboring corn planting anomaly") +
  ylab("County-level corn planting anomaly")+
  geom_hline(yintercept=0,lty=2)+
  ggtitle("Counties with pigs")+
  theme(text=element_text(size=15))+
  ylim(c(-150,300))
ggsave("./Model outputs/Plots/op_counties_nb.jpeg",
       width=5,height=5,units="in")

### previous ----------------------
anom_pr <- data.frame(sc=seq(range(dat_clean_op$plant.anom.prev.sc)[1],
                             range(dat_clean_op$plant.anom.prev.sc)[2],by=0.5))
anom_pr$bt <- anom_pr$sc *attr(dat_clean_op$plant.anom.prev.sc, 'scaled:scale') + 
  attr(dat_clean_op$plant.anom.prev.sc, 'scaled:center')

beta_pr <- beta_long %>% filter(cov=="Previous year's planting anomaly") %>% select(value)
beta_pr <- unlist(c(beta_pr),use.names = F)
est_pr_sc <- sapply(1:length(beta_pr),function(i)beta_int[i] + beta_pr[i]*anom_pr$sc)
est_pr_sum <- data.frame(mn=rowMeans(est_pr_sc),
                         lci=sapply(1:nrow(est_pr_sc),function(i)quantile(est_pr_sc[i,],probs=0.025)),
                         uci=sapply(1:nrow(est_pr_sc),function(i)quantile(est_pr_sc[i,],probs=0.975)))
est_pr_sum <- cbind.data.frame(anom_pr,est_pr_sum)

ggplot(est_pr_sum)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Previous year's corn planting anomaly") +
  ylab("County-level corn planting anomaly")+
  geom_hline(yintercept=0,lty=2)+
  ggtitle("Counties with pigs")+
  theme(text=element_text(size=15))+
  ylim(c(-150,300))
ggsave("./Model outputs/Plots/op_counties_prev.jpeg",
       width=5,height=5,units="in")

### years post nfsp -------------------------
anom_yrs <- data.frame(sc=unique(dat_clean_op$yrs.post.nfsp.sc))

beta_yrs <- beta_long %>% filter(cov=="Years post NFSP") %>% select(value)
beta_yrs <- unlist(c(beta_yrs),use.names = F)

est_yrs_sc <- sapply(1:length(beta_yrs),function(i)beta_int[i] + beta_yrs[i]*anom_yrs$sc)
est_yrs_sum <- data.frame(mn=rowMeans(est_yrs_sc),
                          lci=sapply(1:nrow(est_yrs_sc),function(i)quantile(est_yrs_sc[i,],probs=0.025)),
                          uci=sapply(1:nrow(est_yrs_sc),function(i)quantile(est_yrs_sc[i,],probs=0.975)))
est_yrs_sum <- cbind.data.frame(anom_yrs,est_yrs_sum)
est_yrs_sum$bt <- anom_yrs$sc *attr(dat_clean_op$yrs.post.nfsp.sc, 'scaled:scale') + 
  attr(dat_clean_op$yrs.post.nfsp.sc, 'scaled:center')
est_yrs_sum$bt <- factor(est_yrs_sum$bt)

ggplot(est_yrs_sum)+
  geom_errorbar(aes(x=bt,ymax=uci,ymin=lci))+
  geom_point(aes(x=bt,y=mn))+
  geom_hline(yintercept = 0,lty=2)+
  xlab("Years after NFSP") +
  ylab("County-level corn planting anomaly")+
  ggtitle("Counties with pigs")+
  theme(text=element_text(size=15))
ggsave("./Model outputs/Plots/op_counties_yrs_nfsp.jpeg",
       width=8,height=5,units="in")


### years post nfsp x longitude  ----------------------
anom_long <- data.frame(sc=seq(range(dat_clean_op$long.sc)[1],
                               range(dat_clean_op$long.sc)[2],by=0.5))
anom_yrs <- data.frame(sc=unique(dat_clean_op$yrs.post.nfsp.sc))
anom_long_yrs <- expand.grid(long.sc=anom_long$sc,yrs.sc=anom_yrs$sc)

beta_l <- beta_long %>% filter(cov=="Longitude") %>% select(value)
beta_yrs <- beta_long %>% filter(cov=="Years post NFSP") %>% select(value)
beta_lyrs <- beta_long %>% filter(cov=="Longitude x Years post NFSP") %>% select(value)
beta_l <- unlist(c(beta_l),use.names = F)
beta_yrs <- unlist(c(beta_yrs),use.names = F)
beta_lyrs <- unlist(c(beta_lyrs),use.names = F)

est_lyrs_sc <- sapply(1:length(beta_long),function(i)beta_int[i] + 
                        beta_l[i]*anom_long_yrs$long.sc + 
                        beta_yrs[i]*anom_long_yrs$yrs.sc +
                        beta_lyrs[i]*anom_long_yrs$yrs.sc*anom_long_yrs$long.sc)
est_lyrs_sum <- data.frame(mn=rowMeans(est_lyrs_sc),
                           lci=sapply(1:nrow(est_lyrs_sc),function(i)quantile(est_lyrs_sc[i,],probs=0.025)),
                           uci=sapply(1:nrow(est_lyrs_sc),function(i)quantile(est_lyrs_sc[i,],probs=0.975)))
est_lyrs_sum <- cbind.data.frame(anom_long_yrs,est_lyrs_sum)
est_lyrs_sum$long.bt <- anom_long_yrs$long.sc *attr(dat_clean_op$long.sc, 'scaled:scale') + 
  attr(dat_clean_op$long.sc, 'scaled:center')
est_lyrs_sum$yrs.bt <- anom_long_yrs$yrs.sc *attr(dat_clean_op$yrs.post.nfsp.sc, 'scaled:scale') + 
  attr(dat_clean_op$yrs.post.nfsp.sc, 'scaled:center')
est_lyrs_sum$yrs.bt <- factor(est_lyrs_sum$yrs.bt)

ggplot(est_lyrs_sum)+
  geom_ribbon(aes(x=long.bt,ymax=uci,ymin=lci,
                  fill=yrs.bt,group=yrs.bt),alpha=0.2)+
  geom_line(aes(x=long.bt,y=mn,col=yrs.bt,group=yrs.bt),lwd=1.3)+
  geom_hline(yintercept=0,lty=2)+
  xlab("Longitude") +
  ylab("County-level corn planting anomaly")+
  ggtitle("Counties with pigs")+
  theme(text=element_text(size=15))+
  scale_color_discrete(guide="legend",name="Years of \nNFSP")+
  scale_fill_discrete(name="Years of \nNFSP")
ggsave("./Model outputs/Plots/op_counties_long_yrs_nfsp.jpeg",
       width=8,height=5,units="in")

