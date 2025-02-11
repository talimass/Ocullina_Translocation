# Adapted from Carpenter et al 2022
# Light and photoacclimatization drive distinct differences between shallow and mesophotic coral communities


#LI-CORE
```{r}

library(ggplot2)

licore_data <- read.csv('LICOR_06242021.csv')
licore_data$Depth <- as.factor(licore_data$Depth)
licore <-ggplot(licore_data, aes(Depth, INPUT1, fill=Reading)) +
  geom_boxplot()
licore

licore_percent_all <- licore_data %>%
  mutate(percent_of_max = INPUT1/311.12200) #Solar noon surface = 311.12200

licore_percent_plot <- ggplot(licore_percent_all, aes(x=Depth,y=percent_of_max)) +
  geom_boxplot()
licore_percent_plot

licore_means <- licore_data %>%            #  THIS IS THE DATA SET USED TO MAKE PAR DATA IN MS
  group_by(Depth, Reading) %>%
  summarise(PAR_means = mean(INPUT1),
            std_dev = sd(INPUT1)) %>%
  ungroup()

licore_percent <- read.csv('licor_percents.csv')

summary(licore_means)

licore_aov <- aov (INPUT1 ~ Reading * Depth, data=licore_data)
anova(licore_aov)
```

