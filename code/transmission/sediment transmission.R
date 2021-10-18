library(ggplot2)
library(ggpubr)
library(rcompanion)
library(MASS)

####data import####
sed <- read.csv("sctld_sediment.csv", head=T)
str(sed)
head(sed)
# keeps species/treatment order as imported
sed$Species.Treatment=factor(sed$Species.Treatment, levels=unique(sed$Species.Treatment)) 

####data normality/transformation####
# plots histogram and q-q plot for test of normality assumptions
hist(sed$total.days)
qqnorm(sed$total.days)
qqline(sed$total.days)

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(sed$total.days)
# not normal

#BoxCox transformation finds best exponent value for data transformation
bc <- boxcox(sed$total.days ~ sed$Species.Treatment)
lambda <- bc$x[which.max(bc$y)]

####statistical tests####
# ANOVA on untransformed data
anova <- aov (total.days ~ Species*Treatment, data =sed)
summary(anova)

# ANOVA on Box-Cox transformed data
anova.bc <- aov(((total.days^lambda-1)/lambda) ~ Species*Treatment, data=sed)
summary(anova.bc)

# Q-Q plot comparison between untransformed and transformed data
qq <- par(pty = "s", mfrow = c(1, 2))
qqnorm(anova$residuals); qqline(anova$residuals)
qqnorm(anova.bc$residuals); qqline(anova.bc$residuals)
# despite deviations from normality in raw data, Q-Q plot of ANOVA residuals appears more normal than with Box-Cox transformation
# proceeding with ANOVA of raw data
# both factors (species and treatment) significant

# Tukey post hoc tests
tukey <- TukeyHSD(anova)
tukey

letters <- data.frame(tukey$`Species:Treatment`)
# need to export table for reformatting
write.csv(letters, file = "post hoc results.csv")
# add column label 'comparison' to first column, and find/replace all ':' with '.'
letters <- read.csv("post hoc results.csv", head = T)
letters

# creates compact letter display of significant pairwise differences for figure
cld <- cldList(p.adj ~ comparison, data = letters, threshold = 0.05)
cld

####time to transmission figure####
sed$Species.Treatment=factor(sed$Species.Treatment, levels=unique(sed$Species.Treatment)) 
# boxplots comparing time to transmission among species/treatments
transmission <-
  ggboxplot(
    sed,
    x = "Species.Treatment",
    y = "total.days",
    color = "grey30",
    # palette = c("cornflowerblue", "orange2"),
    fill = "Species",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Treatment",
           y = "Days",
           fill = 'Species') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right")+
  geom_text(data=cld, aes(x = Group, y=0, label=Letter)) 
transmission

ggsave("time_to_transmission.pdf", plot= transmission, width=6, height=4, units="in", dpi=300)

####transmission rate figure####
rate <- read.csv("transmission_rate.csv", head = T)
rate
rate$treatment=factor(rate$treatment, levels=unique(rate$treatment)) 
transrate <- ggplot(rate, aes(fill=forcats::fct_rev(condition), y=rate, x=treatment)) + 
   scale_fill_manual(values=c("#F8766d","#F8766d")) +
  geom_col(width = 0.5) +
  theme_bw() 
transrate

ggsave("transmission_rate.pdf", plot= transrate, width=7.25, height=1.5, units="in", dpi=300)