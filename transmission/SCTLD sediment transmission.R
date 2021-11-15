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

transmission <- read.csv("sctld sediment transmission.csv", head=T)
str(transmission)
head(transmission)

# subsets the data frame to remove healthy sediment (since no disease observed)
transmission <- subset(transmission, treatment!="hs")
transmission

#### data normality/transformation ####

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(transmission$lesion)
# not normal

# BoxCox transformation finds best exponent value for data transformation
bc <- boxcox(transmission$lesion ~ transmission$species_treatment)
lambda <- bc$x[which.max(bc$y)]

# repeating Shaprio test
shapiro.test((transmission$lesion^lambda-1)/lambda)
# still not normal

# plots histogram and q-q plots for raw and Box-Cox transformed data
pdf("sctld sediment transmission normality.pdf")
par(mfrow=c(2,2))
hist(transmission$lesion)
qqnorm(transmission$lesion)
qqline(transmission$lesion)

hist((transmission$lesion^lambda-1)/lambda)
qqnorm((transmission$lesion^lambda-1)/lambda)
qqline((transmission$lesion^lambda-1)/lambda)
dev.off()
# Box-Cox transformation does not appear to have helped


#### statistical tests ####

# ANOVA on untransformed data
anova <- aov (lesion ~ species*colony*treatment, data =transmission)
summary(anova)

# ANOVA on Box-Cox transformed data
anova.bc <- aov(((lesion^lambda-1)/lambda) ~ species*colony*treatment, data=transmission)
summary(anova.bc)

# Q-Q plot comparison between untransformed and transformed data
pdf("sctld sediment transmission ANOVA.pdf")
par(mfrow=c(1,2))
qq <- par(pty = "s", mfrow = c(1, 2))
qqnorm(anova$residuals); qqline(anova$residuals)
qqnorm(anova.bc$residuals); qqline(anova.bc$residuals)
dev.off()

# despite deviations from normality in raw data, Q-Q plot of ANOVA residuals appears more normal than with Box-Cox transformation
# proceeding with ANOVA of raw data

# Tukey post hoc tests
tukey <- TukeyHSD(anova)
tukey

# creating a dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
letters <- data.frame(tukey$`species:treatment`)
letters$Var <- rownames(letters)
names(letters)[5] <- "comparison"
letters$comparison = str_replace_all(letters$comparison,":","_")
letters

# creates compact letter display of significant pairwise differences for figure
cld <- cldList(p.adj ~ comparison, data = letters, threshold = 0.05)
cld


#### time to transmission figure ####

transmission$species_treatment=factor(transmission$species_treatment, levels=c("Of_dc","Of_bds","Of_ids","Mc_dc","Mc_bds","Mc_ids")) 
# boxplots comparing time to transmission among species/treatments
transmission_plot_lesion <-
  ggboxplot(
    transmission,
    x = "species_treatment",
    y = "lesion",
    color = "grey30",
    palette = c("#EBED9D", "#868659","#CECE88","#EBED9D", "#868659","#CECE88"),
    fill = "species_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Time to transmission (d)",
           fill = 'Treatment') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right")+
  geom_text(data=cld, aes(x = Group, y=0, label=Letter)) 
transmission_plot_lesion

ggsave("sctld sediment transmission.pdf", plot= transmission_plot_lesion, width=6, height=4, units="in", dpi=300)


#### transmission rate figure ####

rate <- read.csv("sctld sediment rate.csv", head = T)
rate
rate$treatment=factor(rate$treatment, levels=unique(rate$treatment)) 
transmission_plot_rate <- ggplot(rate, aes(fill=forcats::fct_rev(condition), y=rate, x=treatment)) + 
   scale_fill_manual(values=c("#EBED9D", "#EBED9D")) +
  geom_col(width = 0.5) +
  theme_bw() 
transmission_plot_rate

ggsave("sctld sediment rate.pdf", plot= transmission_plot_rate, width=7.25, height=1.5, units="in", dpi=300)


#### survivorship ####

transmission$treatment=factor(transmission$treatment, levels=c("dc","bds","ids")) 

# subset dataframe by species
ofav <- subset(transmission, species=="Of")
mcav <- subset(transmission, species=="Mc")

# create survival objects for each species (using the Kaplan-Meier method)
survOf <- Surv(time = ofav$days, event = ofav$status)
survOf

survMc <- Surv(time = mcav$days, event = mcav$status)
survMc

# run survival model for each species
fitOf <- survfit(survOf ~ treatment, data = ofav)
summary(fitOf)

fitMc <- survfit(survMc ~ treatment, data = mcav)
summary(fitMc)

# Kaplan-Meier plots for each species
fill.color<-c("#EBED9D", "#868659","#CECE88")

survival_Ofav<-ggsurvplot(fitOf, data = ofav, pval = TRUE, xlab="Days", ylab="Health probability",
                              conf.int = T, risk.table=T, palette=fill.color,
                              break.time.by=5, xlim=c(0,28), risk.table.y.text = FALSE) + ggtitle("O. faveolata") 
survival_Ofav

survival_Mcav<-ggsurvplot(fitMc, data = mcav, pval = TRUE, xlab="Days", ylab="Health probability",
                          conf.int = T, risk.table=T, palette=fill.color,
                          break.time.by=5, xlim=c(0,28), risk.table.y.text = FALSE) + ggtitle("M. cavernosa") 
survival_Mcav

# hazard ratio by disease treatments
hazOf <- coxph(survOf ~ treatment, data = ofav)
summary(hazOf)

hazMc <- coxph(survMc ~ treatment, data = mcav)
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

ggsave("sctld sediment survivorship.pdf", survival_multiplot, width=10, height=10,dpi = 300)
