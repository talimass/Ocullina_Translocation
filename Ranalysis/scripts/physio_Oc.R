library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(car)
library(cowplot)
library(gridExtra)
library(dplyr)
library(multcompView)
library(tidyverse)
library(patchwork)
library(effsize)

setwd("/home/gospozha/haifa/O.patagonica/physio/")

#### physiology ####
# opening file, checking, replacing characters to vectors when needed
physio <- read.csv('physiological_data.csv', stringsAsFactors = F)
str(physio)
#View(physio)
physio$Depth = as.factor(physio$Depth)
physio <- na.omit(physio)
# adding the same theme to each plot
mytheme = theme_bw()+
  theme(plot.title = element_text(size = 10), legend.position = "none",
        axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 9))

# function to compute Hedges' g for pairwise comparisons
pairwise_hedges <- function(df, value_var, group_var) {
  # unique groups
  groups <- unique(df[[group_var]])
  
  # pairwise combinations of groups
  combs <- combn(groups, 2, simplify = FALSE)
  
  # loop over comparisons and calculate Hedges' g
  results <- lapply(combs, function(x) {
    g <- cohen.d(df[[value_var]][df[[group_var]] == x[1]],
                 df[[value_var]][df[[group_var]] == x[2]],
                 hedges.correction = TRUE)
    data.frame(
      Comparison = paste(x[1], "vs", x[2]),
      g = g$estimate,
      lower = g$conf.int[1],
      upper = g$conf.int[2]
    )
  })
  
  do.call(rbind, results)
}


#### growth ####
physio.growth <- physio[physio$growrth_cm.2>0,]
shapiro.test(physio.growth$growrth_cm.2)
leveneTest((growrth_cm.2)~Depth,d=physio.growth)
shapiro.test(sqrt(physio.growth$growrth_cm.2)) # data is normal (p.value > 0.05)
leveneTest(sqrt(growrth_cm.2)~Depth,d=physio.growth)  # no heteroscedasticity of variance (p.value > 0.05)

# anova + Tukey post-hoc 
model <- lm(data = physio.growth, sqrt(growrth_cm.2) ~ Depth)
summary(model)  
plot(model) # Q-Q plot is ok
anova(model) # Depth is not significant factor 

# boxplot

growth <- ggplot(physio.growth, aes(y = (growrth_cm.2), x = Depth)) +
  geom_boxplot(outlier.size = 0.5, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  labs(x = "Depth (m)", y= ~cm^2, title = "Coral growth rate")+
  mytheme

growth


# Pairwise Hedges' g
physio.growth$sqrt_growth <- sqrt(physio.growth$growrth_cm.2)
hedges_df <- pairwise_hedges(physio.growth, "sqrt_growth", "Depth")
hedges_df

# a column with formatted text labels
hedges_df$label <- sprintf("g = %.2f [%.2f, %.2f]",
                           hedges_df$g,
                           hedges_df$lower,
                           hedges_df$upper)

# forest plot with text labels
growth_hedges <- ggplot(hedges_df, aes(x = Comparison, y = g)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = label), 
            hjust = 0.7, vjust = -2, size = 3.2) +  
  labs(y = "Hedges' g (Effect Size, sqrt growth)", x = "Depth comparison", 
       title="Depth effect size on coral growth") +
  coord_flip() +
  mytheme 

growth_hedges

#### protein ####

shapiro.test(physio$protein_ug_ml) #  normal  
leveneTest(protein_ug_ml~Depth,d=physio) # variance is ok

ggplot(physio, aes(x = protein_ug_ml, fill=Depth, alpha=0.1)) +
  geom_density() +
  mytheme

# anova + Tukey post-hoc 
model <- lm(data = physio, protein_ug_ml ~ Depth)
summary(model)  
plot(model) # Q-Q plot is ok
anova(model) # Depth is significant factor 

# Tukey posthoc test
posthoc <- TukeyHSD(aov(model)) 
posthoc$Depth # shows pairwise comparisons
letters <- multcompLetters4(aov(model), posthoc)
letters.df <- data.frame(letters$Depth$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Depth <- rownames(letters.df) 
placement <- physio %>%
  group_by(Depth) %>%
  summarise(quantile(protein_ug_ml, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement)

# double check with non-parametric test - the same results
kruskal.test(protein_ug_ml ~ Depth, data = physio) # Depth is a significant factor 
dunn.res <- dunnTest(protein_ug_ml ~ Depth,
                     data=physio,
                     method="bonferroni")
dunn.res

# boxplot with letters

protein <- ggplot(physio, aes(y = (protein_ug_ml), x = Depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  labs(x = "Depth (m)", y= "Protein concentration (ug/ml)", title = "Coral host protein concentration")+
  mytheme +
  geom_text(data = letters.df, aes(x = Depth, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

protein

hedges_df <- pairwise_hedges(physio, "protein_ug_ml", "Depth")
hedges_df

# a column with formatted text labels
hedges_df$label <- sprintf("g = %.2f [%.2f, %.2f]",
                           hedges_df$g,
                           hedges_df$lower,
                           hedges_df$upper)

# forest plot with text labels
protein_hedges <- ggplot(hedges_df, aes(x = Comparison, y = g)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = label), 
            hjust = 0.7, vjust = -2, size = 3.2) +  
  labs(y = "Hedges' g (Effect Size, protein)", x = "Depth comparison", 
       title="Depth effect size on protein concentration") +
  coord_flip() +
  mytheme 

protein_hedges

#### cell/cm2 ####
# normality, homoscedasticity
shapiro.test(physio$zoox_cm.2) #  normal  
leveneTest(zoox_cm.2~Depth,d=physio) # variance is ok

ggplot(physio, aes(x = (zoox_cm.2), fill=Depth, alpha=0.1)) +
  geom_density() +
  mytheme

# anova  
model <- lm(data = physio, (zoox_cm.2) ~ Depth)
summary(model)  
plot(model)
anova(model) # Depth is not significant factor 

# double check with kruskal-wallis
kruskal.test(zoox_cm.2 ~ Depth, data = physio)  # Depth is not significant factor 

# posthoc: Dunn's Test with Bonferroni correction for p-values
dunn.res <- dunnTest(zoox_cm.2 ~ Depth,
                     data=physio,
                     method="bonferroni")

# boxplot 

cellcm <- ggplot(physio, aes(y = (zoox_cm.2), x = Depth)) +
  geom_boxplot(outlier.size = 0.5, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  labs(x = "Depth (m)", y= ~cells ~x~cm^-2, title = "Symbiont cell count per surface area")+
  mytheme 

cellcm



# Pairwise Hedges' g
hedges_df <- pairwise_hedges(physio, "zoox_cm.2", "Depth")
hedges_df

# a column with formatted text labels
hedges_df$label <- sprintf("g = %.2f [%.2f, %.2f]",
                           hedges_df$g,
                           hedges_df$lower,
                           hedges_df$upper)

# forest plot with text labels
cellcm_hedges <- ggplot(hedges_df, aes(x = Comparison, y = g)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = label), 
            hjust = 0.7, vjust = -2, size = 3.2) +  
  labs(y = "Hedges g (Effect Size, cell/cm^2)", x = "Depth comparison", 
       title="Depth effect size on symbiont cell count") +
  coord_flip()+
  mytheme 

cellcm_hedges

#### chl/ug ####
# normality, homoscedasticity
shapiro.test((physio$chlorophyl_ug_algea)) #  not normal  
leveneTest(chlorophyl_ug_algea~Depth,d=physio) # variance is ok

shapiro.test(sqrt(physio$chlorophyl_ug_algea)) #  not normal  

ggplot(physio, aes(x = (chlorophyl_ug_algea), fill=Depth, alpha=0.1)) +
  geom_density() +
  mytheme

kruskal.test(chlorophyl_ug_algea ~ Depth, data = physio)

# boxplot

chl <- ggplot(physio, aes(y = chlorophyl_ug_algea, x = Depth)) +
  geom_boxplot(outlier.size = 0.5, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  labs(x = "Depth (m)", title="Chlorophyll per ug algae", y="Chlorophyll (~chlorophyll[a] ~x ~ug-1)") +
  mytheme 
 
chl


# Pairwise Hedges' g
hedges_df <- pairwise_hedges(physio, "chlorophyl_ug_algea", "Depth")
hedges_df

# a column with formatted text labels
hedges_df$label <- sprintf("g = %.2f [%.2f, %.2f]",
                           hedges_df$g,
                           hedges_df$lower,
                           hedges_df$upper)

# forest plot with text labels
chl_hedges <- ggplot(hedges_df, aes(x = Comparison, y = g)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = label), 
            hjust = 0.7, vjust = -2, size = 3.2) +  
  labs(y = "Hedges' g (Effect Size, Chl/ug)", x = "Depth comparison", 
       title="Depth effect size on chlorophyll per ug algae") +
  coord_flip() +
  mytheme 

chl_hedges


#### combining plots ####

combined <- (growth + protein) / (cellcm + chl)
combined
ggsave("boxplots_combined.jpg", combined,  width = 8, height = 8)


combined_hedges <- (growth_hedges + protein_hedges) / (cellcm_hedges + chl_hedges)
combined_hedges
ggsave("hedges_combined.jpg", combined_hedges,  width = 8, height = 8)


#### photophysiology ####

fire.data <- read.csv("Annotated_FIRe_data.csv")
fire.data$depth = factor(fire.data$depth, levels = c("4", "10", "25", "4_October"))
fire.data$Colony <- as.factor(fire.data$Colony)
fire.data <- na.omit(fire.data)

mycolors <- c('#C77CFF',  '#F8766D', '#00BA38',"#00BFC4"  )

# reduce dataset 
fire <- fire.data %>%
  drop_na() %>%
  group_by(depth, Colony) %>%
  summarise(
    `Fv.Fm` = mean(`fv_fm`, na.rm = TRUE),
    Sigma = mean(Sigma, na.rm = TRUE),
    `Pmax.e.s` = mean(`Pmax.e.s`, na.rm = TRUE),
    p = mean(p, na.rm = TRUE),
    .groups = "drop"  # optional: ungroups the result
  )

#### Functional absorption cross-section of PSII ####

# checking for normality and homoscedasticity
shapiro.test((fire$Sigma))  # data is normal
leveneTest(Sigma~depth, d=fire) # no heteroscedasticity of variance

ggplot(fire, aes(x = (Sigma), fill=depth, alpha=0.1)) +
  geom_density() +
  mytheme

# plot doesn't look normal - using Kruscal test

kruskal.test(Sigma ~ depth, data = fire) #  not significant

# anova also not significant
model <- lm(data = fire, Sigma ~ depth)
summary(model)  
plot(model) # Q-Q plot is ok
anova(model) # Depth is not significant factor 

# boxplot

sigma <- ggplot(fire, aes(y = (Sigma), x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  scale_fill_manual(values = mycolors )+
  labs(title="Functional absorption cross-section of PSII", y= "σPSII’(A2)", x= "Depth (m)") +
  mytheme

sigma


#### Quantum yield of photochemistry in PSII ####

# checking for normality and homoscedasticity

shapiro.test((fire$Fv.Fm))  # data is not normal (p.value < 0.05)
leveneTest(Fv.Fm~depth, d=fire) # heteroscedasticity is ok

# Kruskal-Wallis test for nonparametric data
kruskal.test(Fv.Fm ~ depth, data = fire) #  significant


# posthoc: Dunn's Test with Bonferroni correction for p-values
dunn.res <- dunnTest(Fv.Fm ~ depth,
                     data=fire,
                     method="bonferroni")

diff <- dunn.res$res$P.adj < 0.05
Names <- gsub(" ", "", dunn.res$res$Comparison)
names(diff) <- Names 
# assigning letters to groups ("compact letter display")
letters <- multcompLetters(diff)
# creating df of letters and their positions to add them to the plot
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$depth <- rownames(letters.df) 
placement <- fire %>% 
  group_by(depth) %>%
  summarise(quantile(Fv.Fm, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters

FvFm = ggplot(fire, aes(y = Fv.Fm, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  scale_fill_manual(values = mycolors )+
  labs(title="Algal quantum yield of photochemistry in PSII", y="Fv’/Fm’") +
  mytheme +
  geom_text(data = letters.df, aes(x = depth, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

FvFm

#### Maximum photosynthetic rate ####
# checking for normality and homoscedasticity
shapiro.test((fire$Pmax.e.s))  # data is not normal (p.value < 0.05)
leveneTest(Pmax.e.s~depth, d=fire) # heteroscedasticity is ok

kruskal.test(Pmax.e.s ~ depth, data = fire) #  significant

# posthoc: Dunn's Test with Bonferroni correction for p-values
dunn.res <- dunnTest(Pmax.e.s ~ depth,
                     data=fire,
                     method="bonferroni")

diff <- dunn.res$res$P.adj < 0.05
Names <- gsub(" ", "", dunn.res$res$Comparison)
names(diff) <- Names 
# assigning letters to groups ("compact letter display")
letters <- multcompLetters(diff)
# creating df of letters and their positions to add them to the plot
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$depth <- rownames(letters.df) 
placement <- fire %>% 
  group_by(depth) %>%
  summarise(quantile(Pmax.e.s, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters

Pmax = ggplot(fire, aes(y = Pmax.e.s, x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  scale_fill_manual(values = mycolors )+
  labs(title="Maximum photosynthetic rate", y="Pmax (electron s-1 PSII-1)") +
  mytheme +
  geom_text(data = letters.df, aes(x = depth, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

Pmax


#### connectivity parameter ####
# checking for normality and homoscedasticity
shapiro.test((fire$p))  # data is not normal (p.value < 0.05)
leveneTest(p~depth, d=fire) # heteroscedasticity

kruskal.test(p ~ depth, data = fire) # not significant

p <- ggplot(fire, aes(y = (p), x = depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  scale_fill_manual(values = mycolors )+
  labs(title="Functional absorption cross-section of PSII", y= "σPSII’(A2)", x= "Depth (m)") +
  mytheme

p

#### combined plots ####
combined <- (FvFm + sigma) / (Pmax + p)
combined
ggsave("boxplots_fire_combined.jpg", combined,  width = 8, height = 8)

combined <- (growth + protein) / (cellcm + chl) / (FvFm + sigma) / (Pmax + p) 
combined
ggsave("boxplots_all_combined.jpg", combined,  width = 8, height = 10)

