### Title: Functional response plots for crop planting analysis
### Author: Abbey Feuka
### Date: 19AUG24
### Notes:

# setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Crops/crop_planting_anom")
setwd("~/Desktop/Crops/crop_planting_anom")
library(tidyverse)

source("./Functions/fun_response.R")
source("./Functions/clean_crop_dat.R")

mod_typ="spatial"
subfolder="CRPxPigs"#'Spatial'
model_dir <- "~/Desktop/Crops/"

#original data --------------------
load(paste0(model_dir,"/Data/all_crops_anom_scaled_2023_2009_anom.RData"))
dat_orig <- dat 

commod_names_c <- unique(dat_orig$commodity_desc)
commod_names_t <- str_to_title(commod_names_c)
commod_names <- tolower(commod_names_c)


#load beta table --------------
beta_tab_orig <- read.csv(paste0(model_dir,"/Model outputs/",
                           subfolder,"/Plots/Combined Figures/op_betas_all_table.csv"))
beta_tab_orig <- beta_tab_orig %>% select(-X) 
beta_tab <- beta_tab_orig %>% 
  filter((grepl("pig",cov) & signif==TRUE))
interaction<- beta_tab_orig %>% 
  filter(crop%in%c("Soybeans","Peanuts") & (cov=="Prop. CRP land" | cov=="Prop. of county with wild pigs"))
beta_tab<- rbind.data.frame(beta_tab,interaction)

fun_take_op<-list()
for(i in 1:nrow(beta_tab)){
  dat_op <- clean_crop_dat(dat_orig=dat,commod_name = tolower(beta_tab$crop[i]),
                           nb=FALSE,
                           temporal=TRUE,
                           only_pigs = T)
  
  if(!grepl("x",beta_tab$cov[i])){ #no interaction
    load(paste0(model_dir,"/Model outputs/",subfolder,"/",
                tolower(beta_tab$crop[i]),"_op_",mod_typ,".RData"))
    
    beta <- samples_all[,grepl("beta",colnames(samples_all)) &
                          !grepl("beta_county",colnames(samples_all)) &
                          !grepl("beta_region",colnames(samples_all))]
    beta <- cbind(beta,
                  chain_idx=samples_all$chain_idx,
                  samp_idx=samples_all$samp_idx)
    
    beta_long_all <- beta %>%
      pivot_longer(grep("beta",colnames(beta)),values_to="value",names_to="beta")
    beta_long_all$cov <- rep(c(#"Intercept",
      dat_op$covsx$name),nrow(beta))
    
    s <- samples_all[,(grepl("s",colnames(samples_all)) | grepl("chain_idx",colnames(samples_all))) &
                       !grepl("lscale",colnames(samples_all)) &
                       !grepl("tau_s",colnames(samples_all))]
    # tau <- samples_all[,"tau"]
    
    fun_take_op[[i]] <- fun_response(dat_df=dat_op$dat_clean,
                                     scale_cov_name=dat_op$covsx$cov[dat_op$covsx$name==beta_tab$cov[i]],
                                     covsx = dat_op$covsx,
                                     beta_long=beta_long_all,
                                     s=s,
                                     # tau=tau,
                                     x_incr=0.5,
                                     spatial=T) 
  }
}

g<-list()
for(i in 1:nrow(beta_tab)){
  g[[i]] <- ggplot(fun_take_op[[i]])+
    geom_ribbon(aes(x=cov_bt,ymin=lci,ymax=uci,fill=county_idx),alpha=0.1)+
    geom_line(aes(x=cov_bt,y=mn,col=county_idx),lwd=1)+
    guides(fill="none",color="none")+
    # xlab(dat_op$covsx$name[dat_op$covsx$name==beta_tab$cov[i]]) +
    xlab("Proportion of county with wild pigs")+
    geom_hline(yintercept=0,lty=2)+
    ylab(expression(paste("County-level planting anomaly (km"^"2",")")))+
    ggtitle(beta_tab$crop[i])+
    theme(text=element_text(size=15),
          panel.background = element_blank(),
          axis.line=element_line(color="grey55"))
}

cowplot::plot_grid(g[[1]],g[[2]])
g[[3]]

ggsave(filename=paste0(model_dir,"/Model outputs/",subfolder,"/Plots/Combined Figures/func_response.jpeg"),
       width=5,height=5,units="in",dpi=800)

#effect size------------
beta_tab<- read.csv(paste0(model_dir,"/Model outputs/",
                           subfolder,"/Plots/Combined Figures/op_betas_all_table.csv"))
beta_tab <- beta_tab %>% select(-X) 
beta_tab$cov_sd <- NA
# beta_tab$cov[beta_tab$cov=="Prop. of county with wild pigs"] <-"Prop. of county with pigs"
# beta_tab$cov[beta_tab$cov=="Take per wild pig intensity 5 yr trend"] <-"Take per hog intensity 5 yr trend"

for(i in 1:nrow(beta_tab)){
  dat_op <- clean_crop_dat(dat_orig=dat,commod_name = tolower(beta_tab$crop[i]),
                           nb=FALSE,
                           temporal=TRUE,
                           only_pigs = T)
  
  dat_df <- as.data.frame(dat_op$dat_clean)
  
  if(beta_tab$cov[i]!="Intercept"){
    
    scale_cov_name <- dat_op$covsx$cov[dat_op$covsx$name==beta_tab$cov[i]]
    beta_tab$cov_sd[i] <- attr(dat_df[,scale_cov_name], 'scaled:scale')
    
  }

}

## proportion of pigs in county
beta_tab %>% filter(cov=="Prop. of county with pigs") %>% 
  filter(signif==TRUE) %>% 
  mutate(cov_eff=md/cov_sd*0.1) #for every 10% change in wild pigs

#temperature 
beta_tab %>% filter(cov=="Temperature 5 yr trend") %>% 
  filter(signif==TRUE) %>% 
  mutate(cov_eff=md/cov_sd) #for every 1 degree increase in 5 year trend

#precip 
beta_tab %>% filter(cov=="Precipitation 5 yr trend") %>% 
  filter(signif==TRUE) %>% 
  mutate(cov_eff=md/cov_sd*10) #for every 1 cm increase in 5 year trend

#crp
beta_tab %>% filter(cov=="Prop. CRP land") %>% 
  filter(signif==TRUE) %>% 
  mutate(cov_eff=md/cov_sd*0.1) #for every 10% change in crp

#crp interaction
beta_tab %>% filter(cov=="CRP x pigs") %>% 
  filter(signif==TRUE) %>% 
  mutate(cov_eff=md/cov_sd*0.1) #for every 10% change in crp

#roi
beta_tab %>% filter(cov=="ROI 5 yr trend") %>% 
  filter(signif==TRUE) %>% 
  mutate(cov_eff=md/cov_sd) #for every 1 USD change in r yr trend

#prev plant anom
beta_tab %>% filter(cov=="Previous year's planting anomaly") %>% 
  filter(signif==TRUE) %>% 
  mutate(cov_eff=md/cov_sd) #for every 1 USD change in r yr trend

