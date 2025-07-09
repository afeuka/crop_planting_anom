library(tidyverse)
library(rnassqs)

#api key
nassqs_auth('7838DEE7-D776-3645-BA19-2FA3803005E0')

#sales
sales <- nassqs(year='2022',agg_level_desc="NATIONAL",statisticcat_desc="SALES",
                group_desc="FIELD CROPS")

total_sales <- sum(sales$Value,na.rm=T)

sales_sum <- sales %>% 
  group_by(commodity_desc) %>% 
  summarise(sales=sum(Value,na.rm = T),
            prop_sales=sales/total_sales) %>% 
  arrange(desc(sales)) %>% 
  mutate(rank=1:length(unique(commodity_desc)))

sales_sum %>% filter(commodity_desc=="PEANUTS")

#acres
acres <- nassqs(year='2022',agg_level_desc="NATIONAL",statisticcat_desc="AREA HARVESTED",
                group_desc="FIELD CROPS")
total_km2 <- acres %>% 
  filter(unit_desc=="ACRES") %>% 
  summarise(total_acres=sum(Value,na.rm = T)/247.1)

km2_sum <- acres %>% 
  filter(unit_desc=="ACRES") %>% 
  group_by(commodity_desc) %>% 
  summarize(km2=sum(Value,na.rm = T)/247.1,
            prop_acres=km2/total_km2) %>% 
  arrange(desc(km2)) %>% 
  mutate(rank=1:length(unique(commodity_desc)))

km2_sum %>% 
  filter(commodity_desc=="PEANUTS") %>% 
  select(km2,rank,prop_acres)

#yield
yield <- nassqs(year='2022',agg_level_desc="NATIONAL",statisticcat_desc="YIELD",
                group_desc="FIELD CROPS")

yield_sum <- yield %>% 
  group_by(commodity_desc,unit_desc) %>% 
  summarize(total_yield=sum(Value,na.rm = T)) %>% 
  arrange(desc(total_yield)) 

yield_sum %>% 
  filter(commodity_desc%in%c("CORN","WHEAT","SOYBEANS","COTTON","PEANUTS")) %>% 
  as.data.frame()

#lb/acre conversion
yield_sum %>% filter(commodity_desc%in%c("PEANUTS", "COTTON")) %>% 
  mutate(kg_acre=total_yield/2.205) %>% 
  mutate(kg_km2=kg_acre/(1/247.1))

#bu/acre conversion
yield_sum %>% filter(commodity_desc%in%c("CORN", "WHEAT","SOYBEANS")) %>% 
  mutate(bu_km2 =total_yield/2.205)

