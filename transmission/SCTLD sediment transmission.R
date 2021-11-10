#### packages ####

# install.packages("survival")
# install.packages("survminer")
library(ggplot2)
library(ggpubr)
library(rcompanion)
library(MASS)
library(stringr)
library(survival)
library(survminer)


#### data import ####

sed <- read.csv("sctld_sediment.csv", head=T)
sed$Genotype <- as.factor(sed$Genotype)
str(sed)
head(sed)
# keeps species/treatment order as imported
sed$Species.Treatment=factor(sed$Species.Treatment, levels=unique(sed$Species.Treatment)) 


#### data normality/transformation ####

# plots histogram and q-q plot for test of normality assumptions
hist(sed$total.days)
qqnorm(sed$total.days)
qqline(sed$total.days)

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(sed$total.days)
# not normal

# BoxCox transformation finds best exponent value for data transformation
bc <- boxcox(sed$total.days ~ sed$Species.Treatment)
lambda <- bc$x[which.max(bc$y)]


#### statistical tests ####

# ANOVA on untransformed data
anova <- aov (total.days ~ Species*Genotype*Treatment, data =sed)
summary(anova)

# ANOVA on Box-Cox transformed data
anova.bc <- aov(((total.days^lambda-1)/lambda) ~ Species*Genotype*Treatment, data=sed)
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

# creating a dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
letters <- data.frame(tukey$`Species:Treatment`)
letters$Var <- rownames(letters)
names(letters)[5] <- "comparison"
letters$comparison = str_replace_all(letters$comparison,":",".")
letters

# creates compact letter display of significant pairwise differences for figure
cld <- cldList(p.adj ~ comparison, data = letters, threshold = 0.05)
cld


#### time to transmission figure ####

sed$Species.Treatment=factor(sed$Species.Treatment, levels=unique(sed$Species.Treatment)) 
# boxplots comparing time to transmission among species/treatments
transmission <-
  ggboxplot(
    sed,
    x = "Species.Treatment",
    y = "total.days",
    color = "grey30",
    palette = c("#EBED9D", "#868659","#CECE88","#EBED9D", "#868659","#CECE88"),
    fill = "Species.Treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Treatment",
           y = "Days",
           fill = 'Treatment') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right")+
  geom_text(data=cld, aes(x = Group, y=0, label=Letter)) 
transmission

ggsave("time_to_transmission.pdf", plot= transmission, width=6, height=4, units="in", dpi=300)


#### transmission rate figure ####

rate <- read.csv("transmission_rate.csv", head = T)
rate
rate$treatment=factor(rate$treatment, levels=unique(rate$treatment)) 
transrate <- ggplot(rate, aes(fill=forcats::fct_rev(condition), y=rate, x=treatment)) + 
   scale_fill_manual(values=c("#EBED9D", "#EBED9D")) +
  geom_col(width = 0.5) +
  theme_bw() 
transrate

ggsave("transmission_rate.pdf", plot= transrate, width=7.25, height=1.5, units="in", dpi=300)


#### survivorship ####

survivorship <- read.csv("sediment_survivorship.csv", head = T)

# subsetting data to remove healthy sediment treatment since no observed lesions
survivorship <- subset(survivorship, Treatment!="hs")
str(survivorship)
survivorship$Treatment <- factor(survivorship$Treatment, levels=c("dc","bds","ids")) 

# subset by species
ofav <- subset(survivorship, Species=="Of")
mcav <- subset(survivorship, Species=="Mc")

# create survival objects for each species (using the Kaplan-Meier method)
survOf <- Surv(time = ofav$days, event = ofav$status)
survOf

survMc <- Surv(time = mcav$days, event = mcav$status)
survMc

# run survival model for each species
fitOfav <- survfit(survOf ~ Treatment, data = ofav)
summary(fitOfav)

fitMcav <- survfit(survMc ~ Treatment, data = mcav)
summary(fitMcav)

# Kaplan-Meier plots for each species
Fill.colour<-c("#EBED9D", "#868659","#CECE88")

survival_Ofav<-ggsurvplot(fitOfav, data = ofav, pval = TRUE, 
                              conf.int = T, risk.table=T, palette=Fill.colour,
                              break.time.by=5, xlim=c(0,25), risk.table.y.text = FALSE,
                              risk.table.title="Number of fragments at risk") + ggtitle("O. faveolata") 
survival_Ofav

survival_Mcav<-ggsurvplot(fitMcav, data = mcav, pval = TRUE, 
                          conf.int = T, risk.table=T, palette=Fill.colour,
                          break.time.by=5, xlim=c(0,25), risk.table.y.text = FALSE,
                          risk.table.title="Number of fragments at risk") + ggtitle("M. cavernosa") 
survival_Mcav

# hazard ratio by disease treatments
hazOf <- coxph(survOf ~ Treatment, data = ofav)
summary(hazOf)

hazMc <- coxph(survMc ~ Treatment, data = mcav)
summary(hazMc)

# plots
hazard_Ofav <- ggforest(hazOf, data = ofav)

hazard_Mcav <- ggforest(hazMc, data = mcav)

#### survivorship figure ####
survival_multiplot<-ggarrange(survival_Ofav$plot,
                              survival_Mcav$plot, 
                              survival_Ofav$table, 
                              survival_Mcav$table,
                              hazard_Ofav,
                              hazard_Mcav,
                              heights = c(2, 0.5, 0.75),
                              ncol = 2, nrow = 3)
survival_multiplot

ggsave("survivorship.pdf", survival_multiplot, width=10, height=10,dpi = 300)

