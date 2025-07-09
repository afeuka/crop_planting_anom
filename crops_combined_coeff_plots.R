### Title: Make combined betas plot
### Author: Abbey Feuka
### Date: 04NOV24

library(tidyverse)

# setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")

source("./crop_planting_anom/Functions/clean_crop_dat.R")

load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

nb=FALSE
temporal=TRUE
mod_typ <- "spatial"

if("take.hog.intens.5yeartrend.sc"%in%colnames(dat_clean) & mod_typ!="spatial"){
  subfolder <- "Take Trend"
} else if("take.hog.intens.prev.sc"%in%colnames(dat_clean) & mod_typ!="spatial") {
  subfolder <- "Take Previous" 
} else if("crp.nfsp.sc"%in%colnames(dat_clean) & mod_typ=="spatial") {
  subfolder <- "CRPxPigs"
} else if(mod_typ=="spatial") {
  subfolder <- "Spatial" 
} else {
  stop("Must include take covariate. Check clean_crop_dat.R")
}

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)

beta_list <- re_list<- tau_long_trace <- lscale_long_trace <-
  beta_long_trace <- list()
if(mod_typ=="spatial"){tau_s_long_trace<-list()}

commod_idx <-1

for(commod_idx in 1:length(commod_names_c)){
  
  ## load samples (only pigs)-----------------
  load(paste0("./Model outputs/",subfolder,"/",commod_names[commod_idx],"_op_",mod_typ,".RData"))
  ## load data------------
  dat_op <- clean_crop_dat(dat_orig=dat_orig,
                           nb=nb,
                           temporal=temporal,
                           commod_name = commod_names[commod_idx],
                           only_pigs = T)
  
  dat_clean_op <- dat_op$dat_clean
  covsx_op <- dat_op$covsx

  
  beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                        !grepl("beta_county",colnames(samples_all)) &
                        !grepl("beta_region",colnames(samples_all))]
  beta <- cbind(beta,
                chain_idx=samples_all$chain_idx,
                samp_idx=samples_all$samp_idx)
  
  beta_long_all <- beta %>%
    pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
  beta_long_all$cov <- rep(c(#"Intercept",
                             covsx_op$name),nrow(beta))
  
  # ##fixed effects --------------------------
  beta_sum <- beta_long_all %>%
    group_by(cov) %>%
    summarise(mn=mean(value),
              md=median(value),
              lci=quantile(value,probs=0.025),
              uci=quantile(value,probs=0.975))
  
  beta_sum$signif <- ifelse(sign(beta_sum$lci)==sign(beta_sum$uci),TRUE,FALSE)
  beta_sum$crop <- commod_names_t[commod_idx]
  beta_list[[commod_idx]] <- beta_sum
  
  if(mod_typ!="spatial"){
    ## random effects ---------------------------
    beta_region <- samples_all[,grep("beta_region",colnames(samples_all))]
    beta_region_sum <-data.frame(mn=colMeans(beta_region),
                                 md=apply(beta_region,2,median),
                                 lci=sapply(1:ncol(beta_region),function(i)quantile(beta_region[,i],probs=0.025)),
                                 uci=sapply(1:ncol(beta_region),function(i)quantile(beta_region[,i],probs=0.975)),
                                 division_grp=unique(dat_clean_op$division_grp))
    
    beta_region_sum$signif <- ifelse(sign(beta_region_sum$lci)==sign(beta_region_sum$uci),TRUE,FALSE)
    beta_region_sum$crop <- commod_names_t[commod_idx]
    re_list[[commod_idx]] <- beta_region_sum
  }
  
  # trace plots -----------------
  tau <- cbind.data.frame(samples_all[,"tau"],
                          chain_idx=samples_all$chain_idx,
                          samp_idx=samples_all$samp_idx)
  lscale <- cbind.data.frame(samples_all[,"lscale"],
                             chain_idx=samples_all$chain_idx,
                             samp_idx=samples_all$samp_idx)
  
  tau_long <- tau %>%
    pivot_longer(grep("tau",colnames(tau)),values_to="value",names_to="beta")
  tau_long$crop <- commod_names_t[commod_idx]
  tau_long_trace[[commod_idx]] <- tau_long
  
  lscale_long <- lscale %>%
    pivot_longer(grep("lscale",colnames(lscale)),values_to="value",names_to="beta")
  lscale_long$crop <- commod_names_t[commod_idx]
  lscale_long_trace[[commod_idx]] <- lscale_long
  
  if(mod_typ=="spatial"){
    tau_s <- cbind.data.frame(samples_all[,"tau_s"],
                              chain_idx=samples_all$chain_idx,
                              samp_idx=samples_all$samp_idx)
    tau_s_long <- tau %>%
      pivot_longer(grep("tau_s",colnames(tau_s)),values_to="value",names_to="beta")
    tau_s_long$crop <- commod_names_t[commod_idx]
    tau_s_long_trace[[commod_idx]] <- tau_s_long
  }

  beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                        !grepl("beta_county",colnames(samples_all)) &
                        !grepl("beta_region",colnames(samples_all))]
  beta <- cbind(beta,
                chain_idx=samples_all$chain_idx,
                samp_idx=samples_all$samp_idx)
  
  beta_long_all <- beta %>%
    pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
  beta_long_all$cov <- rep(c(#"Intercept",
    covsx_op$name),nrow(beta))
  beta_long_all$crop <- commod_names_t[commod_idx]
  beta_long_trace[[commod_idx]] <- beta_long_all

}

beta_all <- do.call("rbind",beta_list)
beta_all <- beta_all %>% mutate(mn=round(mn,2),
                                md=round(md,2),
                                lci=round(lci,2),
                                uci=round(uci,2))
write.csv(beta_all,paste0("./Model outputs/",subfolder,"/Plots/Combined Figures/op_betas_all_table.csv"))

if(mod_typ!='spatial'){
  re_all <- do.call("rbind",re_list)
  write.csv(re_all,paste0("./Model outputs/",subfolder,"/Plots/Combined Figures/op_re_all_table.csv"))
}

beta_long_trace_all <- do.call("rbind",beta_long_trace)
tau_long_trace_all <- do.call("rbind",tau_long_trace)
lscale_long_trace_all <- do.call("rbind",lscale_long_trace)

if(mod_typ=="spatial"){
  tau_s_long_trace_all <- do.call("rbind",tau_s_long_trace)
}

# beta_all$cov[beta_all$cov=="Take per hog intensity 5 yr trend"] <- "Take per wild pig intensity 5 yr trend"
# beta_all$cov[beta_all$cov=="Prop. of county with pigs"] <- "Prop. of county with wild pigs"

beta_all$cov <- factor(beta_all$cov,
                       levels=rev(c(#"Intercept",
                                    "Take per wild pig intensity 5 yr trend",
                                    "Prop. of county with wild pigs",
                                    "Prop. CRP land",
                                    "CRP x pigs",
                                    "ROI 5 yr trend",
                                    "Temperature 5 yr trend",
                                    "Precipitation 5 yr trend",
                                    if(temporal)"Previous year's planting anomaly",
                                    if(nb)"Neighboring planting anomaly"))
                       )



ggplot(beta_all %>% filter(cov!="Intercept"))+
  geom_point(aes(y=cov,x=md,alpha=signif),size=2.5)+
  geom_errorbar(aes(y=cov,xmin=lci,xmax=uci,alpha=signif),width=0,lwd=1)+
  geom_vline(xintercept = 0,col="blue",lty=2)+
  scale_alpha_manual(values=c(0.3,1),name="Significant effect")+
  ylab("")+
  xlab("")+
  facet_wrap(.~crop,nrow=1)+
  guides(alpha="none")+
  theme(text=element_text(size=20),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
ggsave(filename = paste0("./Model outputs/",subfolder,"/Plots/Combined Figures/op_betas_all.jpeg"),
       device = "jpeg",
       width=14,height=8,units="in",dpi=800)

if(mod_typ!="spatial"){
  ggplot(re_all)+
    geom_errorbar(aes(y=division_grp,col=division_grp,
                      xmin=lci,xmax=uci,alpha=signif),width=0,lwd=1)+
    geom_point(aes(y=division_grp,x=md,alpha=signif,col=division_grp),size=2.5)+
    geom_vline(xintercept = 0,col="blue",lty=2)+
    scale_color_manual(values=reg_col)+
    scale_alpha_manual(values=c(0.4,1),name="Significant effect")+
    facet_wrap(.~crop,nrow=1)+
    guides(col="none",alpha="none")+
    ylab("Ecoregion")+xlab("Random intercept estimate")+
    theme(text=element_text(size=20),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
  ggsave(filename = paste0("./Model outputs/",subfolder,"/Plots/Combined Figures/re_all.jpeg"),
         device = "jpeg",
         width=14,height=8,units="in",dpi=800)
}


ggplot(tau_long_trace_all)+
  geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
  ylab("tau")+
  facet_wrap(.~crop,nrow=1)+
  scale_color_discrete(name="chain")+
  theme(text=element_text(size=15))
ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/Combined Figures/tau_trace_op_all.jpeg"),
       device="jpeg",width=14,height=8,units="in",dpi=800)

if(mod_typ=="spatial"){
  ggplot(tau_s_long_trace_all)+
    geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
    ylab("tau_spatial")+
    facet_wrap(.~crop,nrow=1)+
    scale_color_discrete(name="chain")+
    theme(text=element_text(size=15))
  ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/Combined Figures/tau_s_trace_op_all.jpeg"),
         device="jpeg",width=14,height=8,units="in",dpi=800)
}

ggplot(lscale_long_trace_all)+
  geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
  ylab("Laplace scale")+
  facet_wrap(.~crop,nrow=1)+
  scale_color_discrete(name="chain")+
  theme(text=element_text(size=15))
ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/Combined Figures/lscale_trace_op_all.jpeg"),
       device="jpeg",width=14,height=8,units="in",dpi=800)

ggplot(beta_long_trace_all %>% filter(beta!="beta[1]"))+
  geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
  geom_hline(yintercept=0,col="blue",lty=2)+
  facet_grid(crop~cov)+
  scale_color_discrete(name="chain")
ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/Combined Figures/beta_trace_op_all.jpeg"),
       device="jpeg",width=17,height=10,units="in",dpi=800)
