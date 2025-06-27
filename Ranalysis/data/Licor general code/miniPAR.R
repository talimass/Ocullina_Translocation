
##### from new data for graphs

library(anytime)
#the whole folder
setwd("/Users/talimass/Desktop/FIRE general code/Licor general code")
list.files("/Users/talimass/Desktop/FIRE general code/Licor general code")


# prepare each dataset from different depths

m10=read.csv("10m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m10$Depth="10m"
m10$Depth=as.factor(m10$Depth)
m10$Local_timespamp=as.POSIXct(m10$Unix_Timestamp)
str(m10)

m25=read.csv("25m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m25$Depth="25m"
m25$Depth=as.factor(m25$Depth)
m25$Local_timespamp=as.POSIXct(m25$Unix_Timestamp)
str(m25)

m45=read.csv("45m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m45$Depth="45m"
m45$Depth=as.factor(m45$Depth)
m45$Local_timespamp=as.POSIXct(m45$Unix_Timestamp)
str(m45)

# merge all depths

m_all=rbind(m10, m25, m45)
str(m_all)

# add variable of the date only (up to day)
m_all <- m_all %>%
  mutate(date = as.Date(Local_timespamp))
str(m_all)

# add column for time up to hour
m_all$Local_timespamp_hour=as.character(m_all$Local_timespamp)
substr(m_all$Local_timespamp_hour, 15,16) <- "00"
m_all$Local_timespamp_hour <- as.POSIXct(m_all$Local_timespamp_hour)

str(m_all)

# examine the first and last date for each of the depths

m_all %>%
  group_by(Depth) %>%
  summarise(min=min(Local_timespamp),
    max=max(Local_timespamp))

# pase data to contain only dates measured for all depths

library(dplyr)

m_all_filt=m_all %>% 
  #  filter(date >= as.POSIXct("2021-07-19 11:30:30"))
  filter(Local_timespamp >= "2024-06-20 11:19:00",
         Local_timespamp <="2024-09-29 13:39:00")
str(m_all_filt)

# calculate the mean temperature and mean PAR per hour
m_all_filt_means_hour = as.data.frame(m_all_filt %>% 
                              group_by(Depth,Local_timespamp_hour) %>% 
                              summarise(mean_hour_Temperature = mean(Temperature), mean_hour_PAR = mean(PAR)) %>%
                              ungroup())

str(m_all_filt_means_hour)


# calculate the mean temperature and mean PAR per day
m_all_filt_means_day = as.data.frame(m_all_filt %>% 
  group_by(Depth,date) %>% 
  summarise(mean_day_Temperature = mean(Temperature), mean_day_PAR = mean(PAR)) %>%
  ungroup())

str(m_all_filt_means_day)

# subset the merged table (m_all) to the rows with the maximum PAR value
PAR_max <- as.data.frame(m_all_filt %>%
  group_by(Depth,date) %>%
  slice(which.max(PAR)) %>%
  ungroup())

str(PAR_max)

#### plots


library(ggplot2)
str(m_all_filt_means_hour)

scaleFactor <- max(m_all_filt_means_hour$mean_hour_Temperature) / max(m_all_filt_means_hour$mean_hour_PAR)

ggplot(m_all_filt_means_hour, aes(x=Local_timespamp_hour)) +
  geom_line(aes(y=mean_hour_Temperature),  col="blue") +
  geom_line(aes(y=mean_hour_PAR * scaleFactor),  col="red") +
  scale_y_continuous(name="mean_hour_Temperature", sec.axis=sec_axis(~./scaleFactor, name="mean_hour_PAR")) +
  scale_x_datetime(labels = date_format("%Y %b"))+
  theme(
    axis.title.y.left=element_text(color="blue"),
    axis.text.y.left=element_text(color="blue"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red")
  )+
  facet_wrap(Depth~., ncol=1, scales="free")


str(m_all_filt_means_day)

scaleFactor <- max(m_all_filt_means_day$mean_day_Temperature) / max(m_all_filt_means_day$mean_day_PAR)

ggplot(m_all_filt_means_day, aes(x=date)) +
  geom_line(aes(y=mean_day_Temperature),  col="blue") +
  geom_line(aes(y=mean_day_PAR * scaleFactor),  col="red") +
  scale_y_continuous(name="mean_day_Temperature", sec.axis=sec_axis(~./scaleFactor, name="mean_day_PAR")) +
  scale_x_date(labels = date_format("%Y %b"))+
  theme(
    axis.title.y.left=element_text(color="blue"),
    axis.text.y.left=element_text(color="blue"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red")
  )+
  facet_wrap(Depth~., ncol=1, scales="free")



str(PAR_max)

scaleFactor <- max(PAR_max$PAR) / max(PAR_max$Temperature)

ggplot(PAR_max, aes(x=date)) +
  geom_line(aes(y=PAR),  col="blue") +
  geom_line(aes(y=Temperature  * scaleFactor),  col="red") +
  scale_y_continuous(name="PAR", sec.axis=sec_axis(~./scaleFactor, name="Temperature")) +
  scale_x_date(labels = date_format("%Y %b"))+
  theme(
    axis.title.y.left=element_text(color="blue"),
    axis.text.y.left=element_text(color="blue"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red")
  )+
  facet_wrap(Depth~., ncol=1, scales="free")

