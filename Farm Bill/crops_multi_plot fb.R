### Title: Plotting multi-chain crop analysis
### Author: Abbey Feuka
### Date: 15AUG24
### Notes: 

library(nimble)
library(tidyverse)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops")

source("./crop_planting_anom/Functions/clean_crop_dat.R")

nb=FALSE
temporal=FALSE

if(nb==FALSE & temporal==FALSE){
  subfolder <- "Farm Bill Sites No NB No Temp"
} else if (nb==FALSE & temporal==TRUE){
  subfolder <- "Farm Bill Sites No NB"
} else if (nb==TRUE & temporal==FALSE){
  subfolder <- "Farm Bill Sites No Temp"
} else {
  subfolder <- "Farm Bill Sites"
}


load("./Data/all_crops_anom_scaled_2023_2009_anom.RData")
dat_orig <- dat

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)

count_all <- read.csv("./Data/county_state_names.csv")
count_fb <- read.csv("./Data/FB_Projects_Counties_long.csv")
count_fb <- count_fb %>% 
  rename(ST_ABBR=State,
         CNTY_NAME=County) %>% 
  select(-X)
count_fb$FB <-TRUE
count_all <- count_all %>% 
  mutate(CNTY_NAME=toupper(CNTY_NAME))
count_all <- count_all %>% left_join(count_fb)
count_all$FB[is.na(count_all$FB)] <- FALSE

count_all %>% filter(grepl("CLAIR",CNTY_NAME))
count_fb %>% filter(grepl("CLAIR",CNTY_NAME))
dat_orig %>% filter(grepl("CLAIR",county_name) & state_name=="MISSOURI") %>% 
  summarise(unique(FB))

dat_orig_fb <- dat_orig %>% left_join(count_all %>% 
                                        select(-CNTY_NAME) %>% #use geo id to join, not names (. and spaces cause differences)
                                        mutate(state_name=toupper(ST_NAME),
                                               GEOID=factor(GEOID))) %>% 
  filter(FB==TRUE)

# commod_idx <- 5

for(commod_idx in 1:length(commod_names_c)){
  if(commod_idx!=4){
    
    # if("take.hog.intens.5yeartrend.sc"%in%colnames(dat_clean)){
    #   subfolder <- "Take Trend"
    # } else if("take.hog.intens.prev.sc"%in%colnames(dat_clean)) {
    #   subfolder <- "Take Previous" 
    # } else {
    #   stop("Must include take covariate. Check clean_crop_dat.R")
    # }
    
    
    ## random effects ---------------------------------
    ### regions--------------------
    reg_col<- c("Hot Continental"="#F8766D",
                "Warm Continental"="#DD9600",
                "Subtropical"="springgreen4",
                "Tropical/Subtropical Steppe"="#39B600",
                "Prairie"="darksalmon",
                "Marine"="#00B0F6",
                "Temperate Desert"="goldenrod4",
                "Temperate Steppe"="hotpink4",
                "Tropical/Subtropical Desert"="chocolate1",
                "Tropical/Subtropical"="#A3A500",
                "Mediterranean"="paleturquoise4")
    
    # only counties with pigs -------------
    ## load samples -----------------
    load(paste0("./Model outputs/",subfolder,"/",commod_names[commod_idx],"_op_region_multi.RData"))
    
    # create directories 
    if(!file.exists(paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx]))){
      dir.create(paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx]))
    }
    
    ## load data------------
    dat_op <- clean_crop_dat(dat_orig=dat_orig_fb,
                             nb=nb,
                             temporal=temporal,
                             commod_name = commod_names[commod_idx], only_pigs = T)
    
    dat_clean_op <- dat_op$dat_clean
    covsx_op <- dat_op$covsx
    
    ## trace plots -------------------------
    tau <- cbind.data.frame(samples_all[,"tau"],
                            chain_idx=samples_all$chain_idx,
                            samp_idx=samples_all$samp_idx)
    lscale <- cbind.data.frame(samples_all[,"lscale"],
                               chain_idx=samples_all$chain_idx,
                               samp_idx=samples_all$samp_idx)
    
    tau_long <- tau %>%
      pivot_longer(grep("tau",colnames(tau)),values_to="value",names_to="beta")
    lscale_long <- lscale %>%
      pivot_longer(grep("lscale",colnames(lscale)),values_to="value",names_to="beta")
    
    ggplot(tau_long)+
      geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
      ylab("tau")+
      scale_color_discrete(name="chain")+
      theme(text=element_text(size=15))
    ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/tau_trace_op.jpeg"),
           device="jpeg",width=7,height=5,units="in")
    
    ggplot(lscale_long)+
      geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
      ylab("logistic scale")+
      scale_color_discrete(name="chain")+
      theme(text=element_text(size=15))
    ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/lscale_trace_op.jpeg"),
           device="jpeg",width=7,height=5,units="in")
    
    beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                          !grepl("beta_county",colnames(samples_all)) &
                          !grepl("beta_region",colnames(samples_all))]
    beta <- cbind(beta,
                  chain_idx=samples_all$chain_idx,
                  samp_idx=samples_all$samp_idx)
    
    beta_long_all <- beta %>%
      pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
    beta_long_all$cov <- rep(c("Intercept",covsx_op$name),nrow(beta))
    
    ggplot(beta_long_all %>% filter(beta!="beta[1]"))+
      geom_line(aes(x=samp_idx,y=value,col=factor(chain_idx)))+
      geom_hline(yintercept=0,col="blue",lty=2)+
      facet_wrap(.~cov)+
      ggtitle(paste(commod_names_t[commod_idx],"- Counties with pigs"))+
      scale_color_discrete(name="chain")+
      theme(text=element_text(size=15))
    ggsave(filename=paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/beta_trace_op.jpeg"),
           device="jpeg",width=7,height=5,units="in")
    
    # ##fixed effects --------------------------
    beta_sum <- beta_long_all %>%
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
      ggtitle(paste(commod_names_t[commod_idx],"- Counties with pigs"))+
      theme(text=element_text(size=15))
    ggsave(filename = paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/op_betas.jpeg"),
           device = "jpeg",
           width=8,height=8,units="in")
    
    # ## random effects ---------------------------
    beta_region <- samples_all[,grep("beta_region",colnames(samples_all))]
    beta_region_sum <-data.frame(mn=colMeans(beta_region),
                                 lci=sapply(1:ncol(beta_region),function(i)quantile(beta_region[,i],probs=0.025)),
                                 uci=sapply(1:ncol(beta_region),function(i)quantile(beta_region[,i],probs=0.975)),
                                 division_grp=unique(dat_clean_op$division_grp))
    
    beta_region_sum$signif <- ifelse(sign(beta_region_sum$lci)==sign(beta_region_sum$uci),TRUE,FALSE)
    
    ggplot(beta_region_sum)+
      geom_errorbar(aes(y=division_grp,col=division_grp,
                        xmin=lci,xmax=uci,alpha=signif),width=0,lwd=2)+
      geom_point(aes(y=division_grp,x=mn,alpha=signif,col=division_grp),size=3)+
      geom_vline(xintercept = 0,col="blue",lty=2)+
      scale_color_manual(values=reg_col)+
      scale_alpha_manual(values=c(0.3,1),name="Significant effect")+
      guides(col="none")+
      ylab("Ecoregion")+xlab("Random intercept estimate")+
      theme(text=element_text(size=15),
            axis.line = element_blank(),
            axis.ticks = element_blank())+
      ggtitle(paste(commod_names_t[commod_idx],"- Counties with pigs"))
    
    ggsave(filename = paste0("./Model outputs/",subfolder,"/Plots/",commod_names_t[commod_idx],"/op_regions_re.jpeg"),
           device = "jpeg",
           width=8,height=8,units="in")
  }
}

