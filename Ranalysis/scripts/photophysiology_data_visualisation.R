#### Install useful packages
## Matt Doherty


# some packages may not be in your library if using a different computer than in the research office
# if you highlight library, press ctrl + f you can then replace this with the function "install.packages"
# after doing this run the whole script, press ctrl+ f again, change install.packages to library and run the script
# all functions should then be running smoothly!!

library("devtools")
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("lubridate")
library("chron")
library("plyr")
library("dplyr")
library("tidyr")
library("tidyverse")
library("broom")
library("ggpubr")
library("minpack.lm")
library("ggpmisc")

library("purrr")
library("magrittr")
library("hms")
library("stringr")
library("forcats")
library("data.table")
library("car")

library("RColorBrewer")
library("cowplot")
library("grid")


## Eilat photophysiology
options(scipen = 999)
fire.data<- read_csv("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Ocullina_Translocation/Ranalysis/data/Annotated FIRe data.csv") ## you need to path this to your data file

#FIRe Data

fire.data$depth <- as_factor(fire.data$depth)
fire_data <- fire.data %>% na.omit()

# fvfm

quartiles <- quantile(fire_data$fv_fm, probs = c(.25,.75))
IQR_cor <- IQR(fire_data$fv_fm)

lower <- quartiles[1] - 1.5*IQR_cor
upper <- quartiles[2] + 1.5*IQR_cor

fvfm_no_outlier <- subset(fire_data, fire_data$fv_fm >lower & fire_data$fv_fm< upper)


fvfm <- ggplot(data = fvfm_no_outlier,aes(x= depth, y=fv_fm, fill= Species))+
  geom_boxplot() + theme_classic() + stat_boxplot(geom = 'errorbar', width = 0.75, position = "dodge") +
  labs(title="Quantum yield of photochemistry in PSII", y="Fv’/Fm’") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(limits = rev) +
  coord_flip() 

fvfm

### sigma

quartiles <- quantile(fire_data$Sigma, probs = c(.25,.75))
IQR_cor <- IQR(fire_data$Sigma)

lower <- quartiles[1] - 1.5*IQR_cor
upper <- quartiles[2] + 1.5*IQR_cor

sigma_no_outlier <- subset(fire_data, fire_data$Sigma >lower & fire_data$Sigma< upper)


sigma <- ggplot(data = sigma_no_outlier,aes(x= depth, y=Sigma, fill= Species))+
  geom_boxplot() + theme_classic() + stat_boxplot(geom = 'errorbar', width = 0.75, position = "dodge") +
  labs(title="Functional absorption cross-section of PSII", y= "σPSII’(A2)") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(limits = rev) +
  coord_flip() 
sigma

#### pmax

quartiles <- quantile(fire_data$Pmax.e.s, probs = c(.25,.75))
IQR_cor <- IQR(fire_data$Pmax.e.s)

lower <- quartiles[1] - 1.5*IQR_cor
upper <- quartiles[2] + 1.5*IQR_cor

Pmax.e.s_no_outlier <- subset(fire_data, fire_data$Pmax.e.s >lower & fire_data$Pmax.e.s< upper)

pmax <- ggplot(data = Pmax.e.s_no_outlier,aes(x= depth, y=Pmax.e.s, fill= Species))+
  geom_boxplot() + theme_classic() + stat_boxplot(geom = 'errorbar', width = 0.75, position = "dodge") +
  labs(title="Maximum photosynthetic rate ", y="Pmax (electron s-1 PSII-1)") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(limits = rev) +
  coord_flip()  
pmax

#### p

quartiles <- quantile(fire_data$p, probs = c(.25,.75))
IQR_cor <- IQR(fire_data$p)

lower <- quartiles[1] - 1.5*IQR_cor
upper <- quartiles[2] + 1.5*IQR_cor

p_no_outlier <- subset(fire_data, fire_data$p >lower & fire_data$p< upper)

p <- ggplot(data = p_no_outlier,aes(x= depth, y=p, fill= Species)) +
  geom_boxplot() +  theme_classic() + stat_boxplot(geom = 'errorbar', width = 0.75, position = "dodge") +
  labs(title="connectivity parameter ", y="p") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(limits = rev) +
  coord_flip()
p

ggarrange(fvfm, sigma, pmax, p)

