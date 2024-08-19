###Title: Fitting mixed effects model nimble
###Author: Abbey Feuka
###Date: 03042023
###Notes: all counties
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
  dat <- clean_crop_dat(dat_orig=dat_orig,commod_name = commod_names[commod_idx], only_pigs = F)
  
  dat_clean <- dat$dat_clean
  covsx  <- dat$covsx
  xmat <- dat$xmat
  reg_count_idx <- dat$reg_count_idx
  
  # #correlation check 
  # x_cor<- cor(xmat)
  # x_cor_idx <- as.data.frame(which(abs(x_cor)>0.7 &x_cor!=1,arr.ind = T))
  # x_cor_idx[,1] <- covsx_all$cov[x_cor_idx[,1]]
  # x_cor_idx[,2] <- covsx_all$cov[x_cor_idx[,2]]
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
  niter <- 100000
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
  save(samples_all,file=paste0("./Model outputs/",commod_names[commod_idx],"_region_rs_multi.RData"))
  rm(samp1,samp2,samp3,parSamples)
}


##plots ----------------------
# load("corn_planting_09_22_normal_pig_last_year.RData")

###trace plots -------------------------
ggplot(samples_all) + geom_line(aes(x=samp_idx,y=tau,col=idx))

beta <- samples_all[,(grepl("beta",colnames(samples_all)) & !grepl("beta_county",colnames(samples_all))) |
                      grepl("idx",colnames(samples_all))]
beta <- as.data.frame(beta)
# beta$idx <- samples_all[,"idx"]

beta_cor<- cor(beta[,!grepl("idx",colnames(beta))])
# beta_cor <- cor(beta %>% select(-idx))
beta_cor_idx <- as.data.frame(which(abs(beta_cor)>0.7 &beta_cor!=1,arr.ind = T))
beta_cor_idx[,1] <- c("int",covsx_all$cov)[beta_cor_idx[,1]]
beta_cor_idx[,2] <- c("int",covsx_all$cov)[beta_cor_idx[,2]]
beta_cor_idx

# beta$idx_samp <- 1:nrow(beta)
# beta$idx_samp <- rep(1:(nrow(beta)/3),3)
beta_long <- beta %>%  
  pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
beta_long$cov <- rep(c("Intercept",covsx_all$name),nrow(beta))

ggplot(beta_long %>% filter(beta!="beta[1]"))+
  geom_line(aes(x=samp_idx,y=value,col=idx))+
  geom_hline(yintercept=0,col="blue",lty=2)+
  facet_wrap(.~cov)

# plot(samples_all[,paste0("beta[",which(c("int",covsx_all)=="ever.pigs"),"]")],typ="l")

##model checks ---------------------
ypred <- samples_all[,grep("ypred",colnames(samples_all))]
ypred_sum <- data.frame(mn=colMeans(ypred),
                        lci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.025)),
                        uci=sapply(1:ncol(ypred),function(i)quantile(ypred[,i],probs=0.975)),
                        obs = dat_clean$plant.anom)

ypred_sum$resid <- ypred_sum$obs - ypred_sum$mn
ypred_sum$std_resid <- scale(ypred_sum$resid)

ypred_sum = cbind(ypred_sum,dat_clean)
ypred_sum <- ypred_sum %>% 
  left_join(dat %>% dplyr::select(GEOID,lat,long) %>% distinct())

###posterior predictive distribution ----------------
ggplot(ypred_sum %>% pivot_longer(cols=c("mn","obs"),
                                  values_to="value",
                                  names_to="typ"))+
  geom_density(aes(x=value,fill=typ),alpha=0.5)+
  scale_fill_discrete(name="",labels=c("Post pred",
                                       "Data"))

##simulated model stats --------------------
ypred_stat <- data.frame(mn=rowMeans(ypred),
                         sd=sapply(1:nrow(ypred),function(i)sqrt(var(ypred[i,]))),
                         idx=1:nrow(ypred))
ggplot()+
  geom_density(data=ypred_stat,aes(x=mn))+
  geom_vline(xintercept=mean(dat_clean$plant.anom),col="red",lty=2)

ggplot()+
  geom_density(data=ypred_stat,aes(x=sd))+
  geom_vline(xintercept=sqrt(var(dat_clean$plant.anom)),col="red",lty=2)


### standardized residuals--------------------
ggplot(ypred_sum)+
  geom_point(aes(x=mn,y=std_resid,size=long),alpha=0.7)+
  geom_hline(yintercept=0,col="blue",lty=2)

covs_outliers <- ypred_sum %>% filter(abs(std_resid)<=5) %>% 
  pivot_longer(cols=c(grep("ever.pigs",colnames(ypred_sum)),
                      grep("sc",colnames(ypred_sum)),
                      grep("lat",colnames(ypred_sum)),
                      grep("long",colnames(ypred_sum))),
               names_to="cov",values_to="value") %>% 
  group_by(cov) %>% 
  summarise(min=min(value),
            max=max(value),
            typ="outlier")

covs_all <- ypred_sum %>% 
  pivot_longer(cols=c(grep("ever.pigs",colnames(ypred_sum)),
                      grep("sc",colnames(ypred_sum)),
                      grep("lat",colnames(ypred_sum)),
                      grep("long",colnames(ypred_sum))),
               names_to="cov",values_to="value") %>% 
  group_by(cov) %>% 
  summarise(min=min(value),
            max=max(value),
            typ="all")

covs_resid <- full_join(covs_outliers,covs_all) 

ggplot(covs_resid)+
  geom_errorbar(aes(x=cov,ymin=min,ymax=max,col=typ),
                position = position_dodge(width=0.5))+
  theme(axis.text.x=element_text(angle=90,hjust=1))

###fixed effects --------------------------
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
  ggtitle("All counties")+
  theme(text=element_text(size=15))

ggsave(filename = "./Model outputs/Plots/all_counties_betas.jpeg",
       device = "jpeg",
       width=8,height=8,units="in")

###random effects ---------------------------------
beta_county <- samples_all[,grep("beta_county",colnames(samples_all))]
beta_county_sum <-data.frame(mn=colMeans(beta_county),
                             lci=sapply(1:ncol(beta_county),function(i)quantile(beta_county[,i],probs=0.025)),
                             uci=sapply(1:ncol(beta_county),function(i)quantile(beta_county[,i],probs=0.975)),
                             GEOID=unique(dat_clean$GEOID))

ggplot(beta_county_sum)+
  geom_point(aes(y=reorder(GEOID,mn),x=mn),alpha=0.5)+
  geom_errorbar(aes(y=reorder(GEOID,mn),xmin=lci,xmax=uci),width=0,alpha=0.3,col="grey44")+
  geom_vline(xintercept = 0,col="blue",lty=2)+
  ylab("County")+xlab("Random effect estimate")+
  theme(text=element_text(size=15),
        axis.text.y=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "./Model outputs/Plots/all_counties_re.jpeg",
       device = "jpeg",
       width=8,height=8,units="in")

##covariate relationships --------------------------------
### neighboring ----------------------
anom_nb <- data.frame(sc=seq(range(dat_clean$plant.anom.nb.sc)[1],
                             range(dat_clean$plant.anom.nb.sc)[2],by=0.5))
anom_nb$bt <- anom_nb$sc *attr(dat_clean$plant.anom.nb.sc, 'scaled:scale') + 
  attr(dat_clean$plant.anom.nb.sc, 'scaled:center')

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
  ggtitle("All counties")+
  geom_hline(yintercept=0,lty=2)+
  theme(text=element_text(size=15))+
  ylim(c(-250,450))
ggsave("./Model outputs/Plots/all_counties_nb.jpeg",
       width=5,height=5,units="in")

### previous ----------------------
anom_pr <- data.frame(sc=seq(range(dat_clean$plant.anom.prev.sc)[1],
                             range(dat_clean$plant.anom.prev.sc)[2],by=0.5))
anom_pr$bt <- anom_pr$sc *attr(dat_clean$plant.anom.prev.sc, 'scaled:scale') + 
  attr(dat_clean$plant.anom.prev.sc, 'scaled:center')

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
  ggtitle("All counties")+
  theme(text=element_text(size=15))+
  ylim(c(-250,450))
ggsave("./Model outputs/Plots/all_counties_prev.jpeg",
       width=5,height=5,units="in")

### crp ----------------------
anom_crp <- data.frame(sc=seq(range(dat_clean$crp.prop.sc)[1],
                              range(dat_clean$crp.prop.sc)[2],by=0.5))
anom_crp$bt <- anom_crp$sc *attr(dat_clean$crp.prop.sc, 'scaled:scale') + 
  attr(dat_clean$crp.prop.sc, 'scaled:center')

beta_crp <- beta_long %>% filter(cov=="Prop. CRP land") %>% select(value)
beta_crp <- unlist(c(beta_pr),use.names = F)
est_crp_sc <- sapply(1:length(beta_pr),function(i)beta_int[i] + beta_pr[i]*anom_crp$sc)
est_crp_sum <- data.frame(mn=rowMeans(est_crp_sc),
                          lci=sapply(1:nrow(est_crp_sc),function(i)quantile(est_crp_sc[i,],probs=0.025)),
                          uci=sapply(1:nrow(est_crp_sc),function(i)quantile(est_crp_sc[i,],probs=0.975)))
est_crp_sum <- cbind.data.frame(anom_crp,est_crp_sum)

ggplot(est_crp_sum)+
  geom_line(aes(x=bt,y=mn))+
  geom_ribbon(aes(x=bt,ymax=uci,ymin=lci),alpha=0.5)+
  xlab("Land in Conservation Reserve Program") +
  geom_hline(yintercept=0,lty=2)+
  ylab("County-level corn planting anomaly")+
  ggtitle("All counties")+
  theme(text=element_text(size=15))+
  ylim(c(-250,450))
ggsave("./Model outputs/Plots/all_counties_crp.jpeg",
       width=5,height=5,units="in")

