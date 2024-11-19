#Written by Cindy Norton 2022
#install.packages('car')
library('ggplot2')
library(dplyr)
library(car)
library(reshape2)
library(gridExtra)
library(tidyverse)


thermal_data_dry <- readxl::read_xlsx('//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/Figures_Metadata_sapwood_sampling_080222.xlsx', sheet = 1) %>% 
  filter(weather == 'dry') %>%
  mutate(
    Heartwood_radius_thermal = as.numeric(Core_radius_cm) - as.numeric(Sapwood_width_thermal_cm),
    Heartwood_radius_visual = as.numeric(Core_radius_cm) - as.numeric(Sapwood_width_visual_cm),
    SH_ratio_thermal = as.numeric(Sapwood_width_thermal_cm) / Heartwood_radius_thermal,
    SH_ratio_visual = as.numeric(Sapwood_width_visual_cm) / Heartwood_radius_visual,
    Sapwood_area_thermal = pi * as.numeric(Core_radius_cm)^2 - pi * Heartwood_radius_thermal^2,
    Sapwood_area_visual = pi * as.numeric(Core_radius_cm)^2 - pi * Heartwood_radius_visual^2,
    Heartwood_area_thermal = pi * Heartwood_radius_thermal^2,
    Heartwood_area_visual = pi * Heartwood_radius_visual^2
  ) %>%
  group_by(Species) %>%
  na.omit() %>%
  mutate(cooksd = cooks.distance(lm(as.numeric(Sapwood_width_thermal_cm) ~ as.numeric(Sapwood_width_visual_cm))))

thermal_data_dry <- thermal_data_dry %>%
  subset(as.numeric(cooksd) <= quantile(as.numeric(cooksd), 0.90,na.rm=TRUE), 
         as.numeric(cooksd) >= quantile(as.numeric(cooksd), 0.10,na.rm=TRUE))

thermal_data_dry_outliers <- thermal_data_dry %>%
  subset(as.numeric(cooksd) <= quantile(as.numeric(cooksd), 0.10,na.rm=TRUE),
         as.numeric(cooksd) >= quantile(as.numeric(cooksd), 0.90,na.rm=TRUE))

core1 <- thermal_data_dry %>% filter(Sample =="1") %>% 
  dplyr::select(Tree_ID, Species,Sapwood_width_thermal_cm,Sapwood_width_visual_cm,DBH_cm) %>% 
  rename(Sapwood_width_thermal_1 = Sapwood_width_thermal_cm ,Sapwood_width_visual_1 = Sapwood_width_visual_cm)%>% select(-c(DBH_cm))
core2 <- thermal_data_dry %>% filter(Sample =="2") %>%
  dplyr::select(Tree_ID, Species,Sapwood_width_thermal_cm,Sapwood_width_visual_cm,DBH_cm) %>% 
  rename(Sapwood_width_thermal_2 = Sapwood_width_thermal_cm ,Sapwood_width_visual_2 = Sapwood_width_visual_cm)

core_difference <- core1%>%left_join(core2, by = c("Tree_ID","Species")) %>%
  transform(thermal_core_diff = as.numeric(Sapwood_width_thermal_1) - as.numeric(Sapwood_width_thermal_2))%>%
  transform(visual_core_diff = as.numeric(Sapwood_width_visual_1) - as.numeric(Sapwood_width_visual_2))

mean_core_difference<-core_difference %>%                                # Specify data frame
  group_by(Species) %>%                         # Specify group indicator
  na.omit()%>%     
  summarise_at(vars(thermal_core_diff),              # Specify column
               list(name = mean))

mean_species_measurements<-thermal_data_dry %>%                                # Specify data frame
  group_by(Species) %>%                         # Specify group indicator
  mutate(Core_radius_cm = as.numeric(Core_radius_cm), 
            Sapwood_width_visual_cm = as.numeric(Sapwood_width_visual_cm),
            Sapwood_width_thermal_cm = as.numeric(Sapwood_width_thermal_cm),
            Sapwood_area_thermal = as.numeric(Sapwood_area_thermal),
            Sapwood_area_visual = as.numeric(Sapwood_area_visual))%>%
  na.omit()%>%     
  summarise_all(mean)

#write.csv(mean_species_measurements,'mean_species_measurements.csv')
#getwd()

mean_difference_visual_thermal<-thermal_data_dry %>%                                # Specify data frame
  group_by(Species) %>%                         # Specify group indicator
  na.omit()%>%     
  transform(sapwood_diff = as.numeric(Sapwood_width_visual_cm) - as.numeric(Sapwood_width_thermal_cm))%>%
  transform(heartwood_diff = as.numeric(Heartwood_radius_visual) - as.numeric(Heartwood_radius_thermal))%>%
  transform(ratio_diff = as.numeric(SH_ratio_visual) - as.numeric(SH_ratio_thermal))%>%
  transform(sArea_diff = as.numeric(Sapwood_area_visual) - as.numeric(Sapwood_area_thermal))%>%
  group_by(Species) %>%                         # Specify group indicator
  summarise_at(vars(c(sapwood_diff,heartwood_diff,ratio_diff,sArea_diff)),              # Specify column
               list(mean=mean,sd= sd))

total_diff_boundary<- thermal_data_dry %>%                                # Specify data frame
  na.omit()%>%     
  transform(core_diff = as.numeric(Sapwood_width_visual_cm) - as.numeric(Sapwood_width_thermal_cm))%>%
  summarise_at(vars(core_diff),              # Specify column
               list(mean=mean,sd= sd,mad= mad))

perc_diff_area_visual_thermal<-thermal_data_dry %>%                                # Specify data frame
  group_by(Species) %>%                         # Specify group indicator
  na.omit()%>%     
  transform(core_diff = 100*(1-(as.numeric(Sapwood_area_thermal)/as.numeric(Sapwood_area_visual))))%>%
  group_by(Species) %>%                         # Specify group indicator
  summarise_at(vars(core_diff),              # Specify column
               list(mean=mean,sd= sd,mad= mad))

percent_diff_area<- thermal_data_dry %>%                                # Specify data frame
  na.omit()%>%     
  transform(core_diff = 100*(1-(as.numeric(Sapwood_area_thermal)/as.numeric(Sapwood_area_visual))))%>%
  summarise_at(vars(core_diff),              # Specify column
               list(mean=mean,sd= sd,mad= mad))

thermal_dry_species_mean <- thermal_data_dry %>%  group_by(Species, Tree_ID)%>%
  mutate_at(c("DBH_cm", "Core_radius_cm","Sapwood_width_thermal_cm","Heartwood_radius_thermal","SH_ratio_thermal", "Sapwood_area_thermal","Heartwood_area_thermal"), as.numeric)%>%
  summarise_at(vars(DBH_cm, Core_radius_cm,Sapwood_width_thermal_cm,Heartwood_radius_thermal,SH_ratio_thermal,Sapwood_area_thermal,Heartwood_area_thermal), mean, na.rm = TRUE)

visual_dry_species_mean <- thermal_data_dry %>%  group_by(Species, Tree_ID) %>%
  mutate_at(c("DBH_cm", "Core_radius_cm","Sapwood_width_visual_cm","Heartwood_radius_visual","SH_ratio_visual","Sapwood_area_visual","Heartwood_area_visual"), as.numeric)%>%
  summarise_at(vars(DBH_cm, Core_radius_cm,Sapwood_width_visual_cm,Heartwood_radius_visual,SH_ratio_visual,Sapwood_area_visual,Heartwood_area_visual), mean, na.rm = TRUE)


species_mean_1 <- left_join(thermal_dry_species_mean,visual_dry_species_mean, by = c("Species", "Tree_ID","DBH_cm"))





thermal_dry <- thermal_data_dry %>% dplyr::select(-c(Sapwood_width_visual_cm,Heartwood_radius_visual,SH_ratio_visual, Sapwood_area_visual,Heartwood_area_visual)) %>% 
  rename(Sapwood_width_cm = Sapwood_width_thermal_cm, Heartwood_width_cm = Heartwood_radius_thermal,SH_ratio = SH_ratio_thermal, Sapwood_area=Sapwood_area_thermal, Heartwood_area = Heartwood_area_thermal)%>%
  mutate(Category = "thermal")

visual_dry <-thermal_data_dry %>% dplyr::select(-c(Sapwood_width_thermal_cm,Heartwood_radius_thermal,SH_ratio_thermal,Sapwood_area_thermal,Heartwood_area_thermal)) %>% 
  rename(Sapwood_width_cm = Sapwood_width_visual_cm,  Heartwood_width_cm = Heartwood_radius_visual,SH_ratio = SH_ratio_visual,Sapwood_area=Sapwood_area_visual, Heartwood_area = Heartwood_area_visual)%>%
  mutate(Category = "visual")


thermal_data_dry_2 <- rbind(thermal_dry,visual_dry) 


species_mean_2 <- thermal_data_dry_2 %>%  group_by(Species, Tree_ID, Category)%>%
  mutate_at(c("DBH_cm", "Core_radius_cm","Sapwood_width_cm","Heartwood_width_cm","Sapwood_area", "Heartwood_area"), as.numeric)%>%
  summarise_at(vars(DBH_cm, Core_radius_cm,Sapwood_width_cm,Heartwood_width_cm,Sapwood_area, Heartwood_area), mean, na.rm = TRUE)


standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

allmean = thermal_data_dry_2 %>% select(c(2,6,9,10,11,12,14))%>%
  mutate(across(-c(1,7), as.numeric))%>%group_by(Species, Category) %>% 
  summarise(across(everything(), list(mean = ~ mean(.), se = ~ standard_error(.))))%>%
  mutate(across(-c(1,2), ~ round(as.numeric(.), 3)))



thermal_data_monsoons <- readxl::read_xlsx('//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/Figures_Metadata_sapwood_sampling_080222.xlsx', sheet = 4)


thermal_data_wet_1 <- readxl::read_xlsx('//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/Figures_Metadata_sapwood_sampling_080222.xlsx', sheet = 1) %>% filter(weather == 'wet')%>%
  transform(Heartwood_radius_thermal = as.numeric(Core_radius_cm) - as.numeric(Sapwood_width_thermal_cm))%>%
  transform(Heartwood_radius_visual = as.numeric(Core_radius_cm) - as.numeric(Sapwood_width_visual_cm))%>%
  na.omit()%>%
  group_by(Species)%>%
  transform(cooksd = cooks.distance(lm(as.numeric(Sapwood_width_thermal_cm) ~ as.numeric(Sapwood_width_visual_cm),na.rm=TRUE)))%>%
  filter(as.numeric(cooksd) <= quantile(as.numeric(cooksd), 0.95,na.rm=TRUE), 
         as.numeric(cooksd) >= quantile(as.numeric(cooksd), 0.05,na.rm=TRUE))



thermal_wet <- thermal_data_wet_1 %>% dplyr::select(-c(Sapwood_width_visual_cm,Heartwood_radius_visual)) %>% 
  rename(Sapwood_width_cm = Sapwood_width_thermal_cm, Heartwood_width_cm = Heartwood_radius_thermal)%>%
  mutate(Category = "thermal")
visual_wet <-thermal_data_wet_1 %>% dplyr::select(-c(Sapwood_width_thermal_cm,Heartwood_radius_thermal)) %>% 
  rename(Sapwood_width_cm = Sapwood_width_visual_cm,  Heartwood_width_cm = Heartwood_radius_visual)%>%
  mutate(Category = "visual")

thermal_data_wet_2 <- rbind(thermal_wet,visual_wet)





















options("scipen"=100, "digits"=2)

































































##############thermal VS VISUAL COMAPRISONS###########
#Figure 1 scatterplot thermal vs visual
#ppm(file = 'P:/RaBET/Results/CELL_NDWI_2019_points_cfelipse.ppm')
windows()
par(mar=c(5,6,4,6))

####Sapwood width A##########
scatterplot(as.numeric(Sapwood_width_thermal_cm) ~ as.numeric(Sapwood_width_visual_cm) | Species, 
            regLine=FALSE, smooth=FALSE, data=thermal_data_dry, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal CSW (cm)", xlab = "Visual CSW (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)

pp <- thermal_data_dry %>% filter(Species == 'PIPO')
ps <- thermal_data_dry %>% filter(Species == 'PISF')
pm <- thermal_data_dry %>% filter(Species == 'PSME')

pp1 = lm(as.numeric(pp$Sapwood_width_thermal_cm) ~ as.numeric(pp$Sapwood_width_visual_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Sapwood_width_thermal_cm) ~ as.numeric(pm$Sapwood_width_visual_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Sapwood_width_thermal_cm) ~ as.numeric(ps$Sapwood_width_visual_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)

my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('topleft', legend = rp, bty = 'n',cex=1,text.width = 2)



#####Heartwood Radius B########
windows()
par(mar=c(5,6,4,6))
scatterplot(as.numeric(Heartwood_radius_thermal) ~ as.numeric(Heartwood_radius_visual) | Species, 
            regLine=FALSE, smooth=FALSE, data=thermal_data_dry, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal NCW (cm)", xlab = "Visual NCW (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)

pp <- thermal_data_dry %>% filter(Species == 'PIPO')
ps <- thermal_data_dry %>% filter(Species == 'PISF')
pm <- thermal_data_dry %>% filter(Species == 'PSME')

pp1 = lm(as.numeric(pp$Heartwood_radius_thermal) ~ as.numeric(pp$Heartwood_radius_visual))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Heartwood_radius_thermal) ~ as.numeric(pm$Heartwood_radius_visual))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Heartwood_radius_thermal) ~ as.numeric(ps$Heartwood_radius_visual))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)

my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('topleft', legend = rp, bty = 'n',cex=1,text.width = 2)





#####SH_ratio C########
windows()
par(mar=c(5,6,4,6))
scatterplot(as.numeric(SH_ratio_thermal) ~ as.numeric(SH_ratio_visual) | Species, 
            regLine=FALSE, smooth=FALSE, data=thermal_data_dry, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal CS-NCR", xlab = "Visual CS-NCR",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)

pp <- thermal_data_dry %>% filter(Species == 'PIPO')
ps <- thermal_data_dry %>% filter(Species == 'PISF')
pm <- thermal_data_dry %>% filter(Species == 'PSME')

pp1 = lm(as.numeric(pp$SH_ratio_thermal) ~ as.numeric(pp$SH_ratio_visual))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$SH_ratio_thermal) ~ as.numeric(pm$SH_ratio_visual))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$SH_ratio_thermal) ~ as.numeric(ps$SH_ratio_visual))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)

my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('topleft', legend = rp, bty = 'n',cex=1,text.width = 2)






###Part of figure 1 or could be figure 2 MEAN CORE - tree datapoint
#####Sapwood Area D########
windows()
par(mar=c(5,6,4,6))
scatterplot(as.numeric(Sapwood_area_thermal) ~ as.numeric(Sapwood_area_visual) | Species, 
            regLine=FALSE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal CSA (cm)", xlab = "Visual CSA (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)

pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Sapwood_area_thermal) ~ as.numeric(pp$Sapwood_area_visual))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Sapwood_area_thermal) ~ as.numeric(pm$Sapwood_area_visual))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Sapwood_area_thermal) ~ as.numeric(ps$Sapwood_area_visual))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)

my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('topleft', legend = rp, bty = 'n',cex=1,text.width = 2)





























thermal_only_data_2 <-thermal_data_dry_2 %>% filter(Category == "thermal")


#Grouped density plot thermal vs visual per species
###each core
#sapwood width
Sapwood_width_cm_ggplot<- ggplot(thermal_only_data_2, aes(x=as.numeric(Sapwood_width_cm), colour=Species, fill=Species)) + 
  geom_density(alpha=0.5, size =1)+
  theme(legend.position="left")+
  scale_fill_manual(values=c("red", "blue",'green')) +  # Add more colors if needed
  scale_colour_manual(values=c("red", "blue",'green')) +  # Add more colors if needed
  
  labs(x="Conducting Sapwood Width (cm)")+ 
  theme(
    legend.position = "left",
    axis.text=element_text(size=20, face="bold"),
    axis.title=element_text(size=20, face="bold"),
    axis.text.x=element_text(size=20, angle=90),
    axis.text.y=element_text(size=20, face="bold"),
    strip.text=element_text(size=20, face="bold"))

#treemean heartwood radius
Heartwood_width_cm_ggplot<-ggplot(thermal_only_data_2, aes(x=as.numeric(Heartwood_width_cm), colour=Species, fill=Species)) + 
  geom_density(alpha=0.5, size =1)+
  theme(legend.position="none")+
  scale_fill_manual(values=c("red", "blue",'green')) +  # Add more colors if needed
  scale_colour_manual(values=c("red", "blue",'green')) +  # Add more colors if needed
  
  labs(x="Non Conducting Wood (cm)")+ 
  theme(
    legend.position = "none",
    axis.text=element_text(size=20, face="bold"),
    axis.title=element_text(size=20, face="bold"),
    axis.text.x=element_text(size=20, angle=90),
    axis.text.y=element_text(size=20, face="bold"),
    strip.text=element_text(size=20, face="bold"))

#SHRatio
SH_ratio_ggplot<-ggplot(thermal_only_data_2, aes(x=as.numeric(SH_ratio), colour=Species, fill=Species)) + 
  geom_density(alpha=0.5, size =1)+
  theme(legend.position="none")+
  scale_fill_manual(values=c("red", "blue",'green')) +  # Add more colors if needed
  scale_colour_manual(values=c("red", "blue",'green')) +  # Add more colors if needed
  
  labs(x="Conducting - Non Conducting Ratio")+    
  theme(
    legend.position = "none",
    axis.text=element_text(size=20, face="bold"),
    axis.title=element_text(size=20, face="bold"),
    axis.text.x=element_text(size=20, angle=90),
    axis.text.y=element_text(size=20, face="bold"),
    strip.text=element_text(size=20, face="bold"))
species_mean_3<-species_mean_2 %>% filter(Category == 'thermal')

Sapwood_area_ggplot <- ggplot(species_mean_3, aes(x=as.numeric(Sapwood_area), colour=Species, fill=Species)) + 
  geom_density(alpha=0.5, size=1) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("red", "blue",'green')) +  # Add more colors if needed
  scale_colour_manual(values=c("red", "blue",'green')) +  # Add more colors if needed
  
  labs(x="Conducting Sapwood Area (cm)") +
  theme(
    legend.position = "none",
    axis.text=element_text(size=20, face="bold"),
    axis.title=element_text(size=20, face="bold"),
    axis.text.x=element_text(size=20, angle=90),
    axis.text.y=element_text(size=20, face="bold"),
    strip.text=element_text(size=20, face="bold"))

windows()

grid.arrange(Sapwood_width_cm_ggplot, Heartwood_width_cm_ggplot,SH_ratio_ggplot, Sapwood_area_ggplot,nrow = 2)



##TREE metrics
#Sapwood area
windows()



Sapwood_area_ggplot



##DBH density plot
windows()
DBH_density_ggplot<-ggplot(data=thermal_dry_species_mean, aes(x=DBH_cm, fill=Species)) +
  geom_density(alpha=0.5, size =1) + 
  xlab("DBH (cm)")+
  theme_minimal()+theme(legend.text = element_text(size=12),strip.text = element_text(size=15),axis.text=element_text(size=13),
                        axis.title=element_text(size=12,face="bold"))
DBH_density_ggplot









#detach(package:plyr)

###BARPLOTS COMPARING cores for each tree
thermal_data_dry_2_sd <-thermal_data_dry_2%>%
  group_by(Category, Tree_ID, Species) %>% summarise(mean = mean(as.numeric(Sapwood_width_cm)),
                                                     std = sd(as.numeric(Sapwood_width_cm)))

pp_thermal_data_dry_2 <- thermal_data_dry_2 %>% filter(Species == 'PIPO') %>%
  left_join(thermal_data_dry_2_sd, by = c("Species","Tree_ID","Category"))

ps_thermal_data_dry_2 <- thermal_data_dry_2 %>% filter(Species == 'PISF')%>%
  left_join(thermal_data_dry_2_sd, by = c("Species","Tree_ID","Category"))

pm_thermal_data_dry_2 <- thermal_data_dry_2 %>% filter(Species == 'PSME')%>%
  left_join(thermal_data_dry_2_sd, by = c("Species","Tree_ID","Category"))

pp_thermal_data_dry_2_barplot <- pp_thermal_data_dry_2 %>% 
  group_by(Category) %>% 
  mutate(se = sd(as.numeric(Sapwood_width_cm))/sqrt(length(as.numeric(Sapwood_width_cm)))) %>% 
  ggplot(aes(x = reorder(as.numeric(DBH_cm),-as.numeric(DBH_cm)), y = as.numeric(Sapwood_width_cm), fill = Category)) + 
  geom_bar(stat="identity", alpha=0.5, 
           position=position_dodge()) +theme_bw()+
  scale_fill_manual(values=c("red",
                             "darkblue"))+
  geom_errorbar(aes(ymin= mean-std, ymax=mean+std), width=.2,position=position_dodge(.9))+
  xlab("DBH (cm)")+
  ylab("CSW")+ggtitle("PIPO") +
  theme(legend.position = "none" ,axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"), axis.text.x = element_text(size=15,angle = 90),
        axis.text.y = element_text(size=15),strip.text = element_text(size = 20,face="bold"))+
  ylim(0,15)


ps_thermal_data_dry_2_barplot <- ps_thermal_data_dry_2 %>% 
  group_by(Category) %>% 
  mutate(se = sd(as.numeric(Sapwood_width_cm))/sqrt(length(as.numeric(Sapwood_width_cm)))) %>% 
  ggplot(aes(x = reorder(as.numeric(DBH_cm),-as.numeric(DBH_cm)), y = as.numeric(Sapwood_width_cm), fill = Category)) + 
  geom_bar(stat="identity", alpha=0.5, 
           position=position_dodge()) +theme_bw()+
  scale_fill_manual(values=c("red",
                             "darkblue"))+
  geom_errorbar(aes(ymin= mean-std, ymax=mean+std), width=.2,position=position_dodge(.9))+
  xlab("DBH (cm)")+
  ylab("CSW")+ggtitle("PISF") +
  theme(legend.position = "none" ,axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"), axis.text.x = element_text(size=15,angle = 90),
        axis.text.y = element_text(size=15),strip.text = element_text(size = 20,face="bold"))+
  ylim(0,15)




pm_thermal_data_dry_2_barplot <- pm_thermal_data_dry_2 %>% 
  group_by(Category) %>% 
  mutate(se = sd(as.numeric(Sapwood_width_cm))/sqrt(length(as.numeric(Sapwood_width_cm)))) %>% 
  ggplot(aes(x = reorder(as.numeric(DBH_cm),-as.numeric(DBH_cm)), y = as.numeric(Sapwood_width_cm), fill = Category)) + 
  geom_bar(stat="identity", alpha=0.5, 
           position=position_dodge()) +theme_bw()+  scale_fill_manual(values=c("red",
                                                                               "darkblue"))+
  geom_errorbar(aes(ymin= mean-std, ymax=mean+std), width=.2,position=position_dodge(.9))+
  xlab("DBH (cm)")+
  ylab("CSW")+ggtitle('PSME') +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"), axis.text.x = element_text(size=15,angle = 90),
        axis.text.y = element_text(size=15),strip.text = element_text(size = 20,face="bold"))+
  ylim(0,15)


windows()
grid.arrange(pp_thermal_data_dry_2_barplot, ps_thermal_data_dry_2_barplot,pm_thermal_data_dry_2_barplot, nrow = 1)



















#######core differences

##lineplot
windows()
ggplot(data=core_difference, aes(x = reorder(as.numeric(DBH_cm),-as.numeric(DBH_cm)), y=thermal_core_diff, group = Species, colour = Species)) +
  geom_line(size=1) +
  geom_point( size=6, shape=19, fill="white")+
  xlab("DBH (cm)")+
  ylab("Difference")+
  theme(axis.text.x = element_text(angle = 90))



#barplot
windows()
ggplot(data=core_difference, aes(x = reorder(as.numeric(DBH_cm),-as.numeric(DBH_cm)), y=thermal_core_diff, fill = Species, colour = Species)) +
  geom_bar(stat="identity", alpha=0.5, 
           position=position_dodge()) +
  xlab("DBH (cm)")+
  ylab("Difference")+
  theme(axis.text.x = element_text(angle = 90))


#  annotate("text", x= 50, y=-3, label= "PIPO = - 0.103 mm") + 
#  annotate("text", x = 50, y=-2.5, label = "PISF = 0.04 mm")+
#  annotate("text", x = 50, y=-2, label = "PM = - 0.316 mm")





### histogram thermal difference per species

# change fill and outline color manually 
ggplot(core_difference, aes(x = thermal_core_diff)) +
  geom_histogram(aes(color = Species,), 
                 position = "identity", bins = 30, alpha = 0.4) +
  facet_wrap(~Species)+
  scale_color_manual(values = c("red", "green","blue")) +
  scale_fill_manual(values = c("red", "green","blue"))+
  xlab("Thermal Core Difference")+
  ylab("Count")+
  
  
  
  
  
  
  
  
  
  
  



















  
  







###########DBH and thermal comaprisons# speices##########
########################Figure 1 scatterplot hearthwood thermal vs DBH####################
#Figure4
#ppm(file = 'P:/RaBET/Results/CELL_NDWI_2019_points_cfelipse.ppm')
#ppm(file = 'P:/RaBET/Results/CELL_NDWI_2019_points_cfelipse.ppm')
windows()
scatterplot(as.numeric(Sapwood_width_thermal_cm) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal CSW (cm)", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)


pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Sapwood_width_thermal_cm) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Sapwood_width_thermal_cm) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Sapwood_width_thermal_cm) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('topright', legend = rp, bty = 'n',cex=1.2)



windows()
scatterplot(as.numeric(Heartwood_radius_thermal) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal NCW (cm)", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)

pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Heartwood_radius_thermal) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Heartwood_radius_thermal) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Heartwood_radius_thermal) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',2)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('bottomright', legend = rp, bty = 'n',cex=1.2)



















#####SH_ratio C########
windows()
scatterplot(as.numeric(SH_ratio_thermal) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal CS-NCR", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)

pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$SH_ratio_thermal) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$SH_ratio_thermal) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$SH_ratio_thermal) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('topright', legend = rp, bty = 'n',cex=1.2)








#ppm(file = 'P:/RaBET/Results/CELL_NDWI_2019_points_cfelipse.ppm')
windows()
scatterplot(as.numeric(Heartwood_area_thermal) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal NCA", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)


pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Heartwood_area_thermal) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Heartwood_area_thermal) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Heartwood_area_thermal) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('bottomright', legend = rp, bty = 'n',cex=1.2)






#ppm(file = 'P:/RaBET/Results/CELL_NDWI_2019_points_cfelipse.ppm')
windows()
scatterplot(as.numeric(Sapwood_area_thermal) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Thermal CSA", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)


pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Sapwood_area_thermal) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Sapwood_area_thermal) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Sapwood_area_thermal) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('bottomright', legend = rp, bty = 'n',cex=1.2)




















windows()
grid.arrange(CSW, NCW,ratio,NCA,CSA, nrow = 3)















































































































#####sapwood width#########

#####################SPECIES tree Sapwood Area vs DBH##########################

####################DBH and visual############


###########DBH and visual comaprisons###########
########################Figure 1 scatterplot hearthwood visual vs DBH####################
#Figure4
#ppm(file = 'P:/RaBET/Results/CELL_NDWI_2019_points_cfelipse.ppm')
windows()
scatterplot(as.numeric(Heartwood_radius_visual) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Heartwood Visual Radius", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)

pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Heartwood_radius_visual) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Heartwood_radius_visual) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Heartwood_radius_visual) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('bottomright', legend = rp, bty = 'n',cex=1.2)






#####################SPECIES tree Sapwood Area vs DBH##########################

#ppm(file = 'P:/RaBET/Results/CELL_NDWI_2019_points_cfelipse.ppm')
windows()
scatterplot(as.numeric(Sapwood_area_visual) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Visual Sapwood Area", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)


pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Sapwood_area_visual) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Sapwood_area_visual) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Sapwood_area_visual) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('bottomright', legend = rp, bty = 'n',cex=1.2)




#sapwood width
windows()
scatterplot(as.numeric(Sapwood_width_visual_cm) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Visual Sapwood Width (cm)", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)


pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Sapwood_width_visual_cm) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Sapwood_width_visual_cm) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Sapwood_width_visual_cm) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('bottomright', legend = rp, bty = 'n',cex=1.2)





#####SH_ratio C########
windows()
scatterplot(as.numeric(SH_ratio_visual) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Visual Sapwood-Heartwood Ratio", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)

pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$SH_ratio_visual) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$SH_ratio_visual) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$SH_ratio_visual) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('bottomright', legend = rp, bty = 'n',cex=1.2)























#####################SPECIES tree Heartwood Area vs DBH##########################

#ppm(file = 'P:/RaBET/Results/CELL_NDWI_2019_points_cfelipse.ppm')
windows()
scatterplot(as.numeric(Heartwood_area_visual) ~ as.numeric(DBH_cm) | Species, 
            regLine=TRUE, smooth=FALSE, data=species_mean_1, pch = c(1,2,3), 
            col = c('red','blue', 'green'), ylab = "Visual NCA", xlab = "DBH (cm)",
            cex.lab=1.5,cex.main=2, cex.axis=1.5, cex = 2)


pp <- species_mean_1 %>% subset(Species == 'PIPO')
ps <- species_mean_1 %>% subset(Species == 'PISF')
pm <- species_mean_1 %>% subset(Species == 'PSME')

pp1 = lm(as.numeric(pp$Heartwood_area_visual) ~ as.numeric(pp$DBH_cm))
ppsum = summary(pp1)
r2_pp = ppsum$adj.r.squared

pm1 = lm(as.numeric(pm$Heartwood_area_visual) ~ as.numeric(pm$DBH_cm))
pmsum = summary(pm1)
r2_pm = pmsum$adj.r.squared

ps1 = lm(as.numeric(ps$Heartwood_area_visual) ~ as.numeric(ps$DBH_cm))
pssum = summary(ps1)
r2_ps = pssum$adj.r.squared

my.p_pm = format(pmsum$coefficients[2,4], scientific = TRUE)
my.p_pp = format(ppsum$coefficients[2,4], scientific = TRUE)
my.p_ps = format(ppsum$coefficients[2,4], scientific = TRUE)


my.slope_pm = pmsum$coefficients[2,1] 
my.slope_pp = ppsum$coefficients[2,1] 
my.slope_ps = pssum$coefficients[2,1] 


rp = vector('expression',12)
rp[1] = substitute(expression('PIPO'),)[2]
rp[2] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pp,dig=3)))[2]
rp[3] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pp, digits = 3)))[2]
rp[4] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pp, digits = 3)))[2]
rp[5] = substitute(expression('PSME'),)[2]

rp[6] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2_pm,dig=3)))[2]
rp[7] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p_pm, digits = 3)))[2]
rp[8] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope_pm, digits = 3)))[2]

rp[9] = substitute(expression('PISF'),)[2]
rp[10] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2_ps,dig=3)))[2]
rp[11] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p_ps, digits = 3)))[2]
rp[12] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.slope_ps, digits = 3)))[2]


legend('topleft', legend = rp, bty = 'n',cex=1,text.width = 2)














#T.test



diff_pp <- core_difference %>% filter(Species == 'PIPO')
diff_ps <- core_difference %>% filter(Species == 'PISF')
diff_pm <- core_difference %>% filter(Species == 'PSME')



t.test(as.numeric(diff_pp$Sapwood_width_thermal_1), as.numeric(diff_pp$Sapwood_width_thermal_2), paired = TRUE, alternative = "two.sided")
t.test(as.numeric(diff_ps$Sapwood_width_thermal_1), as.numeric(diff_ps$Sapwood_width_thermal_2), paired = TRUE, alternative = "two.sided")



wilcox.test(as.numeric(diff_pm$Sapwood_width_thermal_1), as.numeric(diff_pm$Sapwood_width_thermal_2), alternative = "two.sided")



























































####OTHER PLOTS##############
windows()
Sapwood_area_boxplot <- ggplot(species_mean_2, aes(x=Tree_ID, y=Sapwood_area, fill=Category)) + 
  geom_boxplot() +
  facet_wrap(~Species)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Sapwood Area") 
Sapwood_area_boxplot




##Figure 2 Boxplot Species broken down by species and sapwood width
windows()
# One box per treatment
sapwood_width_boxplot <- ggplot(species_mean_1, aes(x=Tree_ID, y=Sapwood_width_cm, fill=Category)) + 
  geom_boxplot() +
  facet_wrap(~Species)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Sapwood Width") 
sapwood_width_boxplot

heartwood_radius_boxplot <- ggplot(species_mean_1, aes(x=Tree_ID, y=Heartwood_width_cm, fill=Category)) + 
  geom_boxplot() +
  facet_wrap(~Species)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Heartwood Radius") 
heartwood_radius_boxplot






attach(mtcars)
windows()
par(mfrow=c(1,2))

boxplot(Sapwood_width_visual_cm ~ Species, 
        data = thermal_data_dry, 
        col=c("pink","light blue", "grey"), 
        ylab = "Visual Width", 
        xlab = "Species", 
        main = "Visual",
        ylim = c(0,15),
        border = NA, 
        xaxt='n', yaxt = "n", frame = FALSE)
grid(nx=17, ny=17)
boxplot(Sapwood_width_visual_cm ~ Species, 
        data = thermal_data_dry, 
        col=c("pink","light blue", "grey"), 
        ylab = "Visual Width", 
        xlab = "Species", 
        main = "Visual",
        ylim = c(0,15),add = TRUE, ann = FALSE)

boxplot(as.numeric(Sapwood_width_thermal_cm) ~ Species, 
        data = thermal_data_dry, 
        col=c("pink","light blue", "grey"), 
        ylab = "thermal Width", 
        xlab = "Species", 
        main = "thermal",
        ylim = c(0,15),
        border = NA, 
        xaxt='n', yaxt = "n", frame = FALSE)
grid(nx=17, ny=17)
boxplot(as.numeric(Sapwood_width_thermal_cm) ~ Species, 
        data = thermal_data_dry, 
        col=c("pink","light blue", "grey"), 
        ylab = "thermal Width", 
        xlab = "Species", 
        main = "thermal",
        ylim = c(0,15),add = TRUE, ann = FALSE)

