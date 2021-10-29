#### packages ####

library(readr)
library(tidyverse)
library(multcomp)
library(emmeans)
library(readr)
library(betareg)
library(emmeans)
library(gridExtra)


#### data import and summary stats ####

histo <- read_csv("SCTLD sediment histology.csv")
histo$Result <- as.numeric(histo$Result)
str(histo)

dis <- vector()

for (val in 1:length(histo$Treatment)) {
  if (is.na(histo$Treatment[val])) {
    dis[val] <- NA
  }
  else if (histo$Treatment[val] == 'IDS') {
    dis[val] <- 1
  } else if (histo$Treatment[val] == 'DC') {
    dis[val] <- 1
  } else if (histo$Treatment[val] == 'BDS') {
    dis[val] <- 1
  } else if (histo$Treatment[val] == 'HS') {
    dis[val] <- 0
  } else {
  }
}
histo$Diseased <- dis

histo %>%
  group_by(Species, Genotype_Replicate, Treatment, Measurement) %>%
  summarize(mean = mean(Result, na.rm = TRUE), sd = sd(Result), diseased = mean(Diseased)) -> summarized_data


#### symbiont:vacuole ratio ######

# filtering data
summarized_data %>%
  filter(Measurement == "Zooxarea") -> vac1

vac1$Species <- factor( as.character(vac1$Species), levels=c("Ofav","Mcav") )
vac1$Treatment <- factor( as.character(vac1$Treatment), levels=c("HS","DC","BDS","IDS") )

# general linear models for both species
# Mcav
fitzoox <- glm(as.logical(diseased) ~ mean, family = "binomial", data = vac1[which(vac1$Species == 'Mcav'),])
summary(fitzoox)

# Ofav
fitzoox <- glm(as.logical(diseased) ~ mean, family = "binomial", data = vac1[which(vac1$Species == 'Ofav'),])
summary(fitzoox)

# making predictions of symbiont:vacuole ratio to health condition
xfitzoox1 <- seq(0,1,0.001)
yfitzoox1 <- predict(fitzoox, list(mean = xfitzoox1),type = "response")

# plot of symbiont:vacuole ratio prediction of disease condition
pdf("histology_symbiontvacuole_prediction.pdf", width=6, height=6)
plot(as.logical(diseased) ~ mean, data = vac1, pch  = 16, 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), 
     xlab = "Symbiont:vacuole ratio", ylab = "Probability diseased", xlim = c(0.25,0.75))
legend(0.25, 0.2, legend=c("HS", "DC", "BDS", "IDS"),
       col=c("#EBED9D", "grey50", "#868659", "#CECE88"), lty=1, cex=0.8)
lines(xfitzoox1,yfitzoox1)

thresh <- max(yfitzoox1[yfitzoox1 <= 0.5])
loc = yfitzoox1[yfitzoox1 == thresh]
loc <- as.numeric(attributes(loc)$names)
maxprop <- xfitzoox1[loc]
segments(0,.5,maxprop,.5,col = 'black', lty = 'dashed')
segments(maxprop,0,maxprop,.5,col = 'black', lty = 'dashed')

xH <- vac1$mean[which(vac1$Treatment == 'HS')]
yH <- predict(fitzoox, list(mean = xH), type = "response")
points(xH,yH,col = '#EBED9D')
sdH <- vac1$sd[which(vac1$Treatment == 'HS')]
sdlowH <- xH-sdH
sdhighH <- xH+sdH
segments(sdlowH,yH,sdhighH,yH,col = '#EBED9D')

xD <- vac1$mean[which(vac1$Treatment == 'DC')]
yD <- predict(fitzoox, list(mean = xD), type = "response")
points(xD,yD,col = 'grey50')
sdC <- vac1$sd[which(vac1$Treatment == 'DC')]
sdlowC <- xD-sdC
sdhighC <- xD+sdC
segments(sdlowC,yD,sdhighC,yD, col = 'grey50')

xB <- vac1$mean[which(vac1$Treatment == 'BDS')]
yB <- predict(fitzoox, list(mean = xB), type = "response")
points(xB,yB,col = '#868659')
sdS <- vac1$sd[which(vac1$Treatment == 'BDS')]
sdlowS <- xB-sdS
sdhighS <- xB+sdS
segments(sdlowS,yB,sdhighS,yB, col = '#868659')

xI <- vac1$mean[which(vac1$Treatment == 'IDS')]
yI <- predict(fitzoox, list(mean = xI), type = "response")
points(xI,yI,col = '#CECE88')
sdD <- vac1$sd[which(vac1$Treatment == 'IDS')]
sdlowD <- xI-sdD
sdhighD <- xI+sdD
segments(sdlowD,yI,sdhighD,yI, col = '#CECE88')
dev.off()

# beta regression for treatment and species comparisons
y <- vac1$mean
n.obs <- sum(!is.na(y))
y2 <- (y * (n.obs -1) + 0.5)/n.obs
vac1$y2 <- y2

beta_vac <- betareg(y2 ~ Treatment + Species, vac1)
# overall factor p values
joint_tests(beta_vac)

# individual species tests
vac_ofav <- betareg(mean ~ Treatment, vac1[which(vac1$Species == 'Ofav'),])
summary(vac_ofav)
vac_mcav <- betareg(mean ~ Treatment, vac1[which(vac1$Species == 'Mcav'),])
summary(vac_mcav)

# boxplot of symbiont:vacuole ratio
specnames <- c("Orbicella faveolata", "Montastraea cavernosa")
names(specnames) <- c("Ofav", "Mcav")

pdf("histology_symbiontvacuole_boxplot.pdf", width=6, height=4)
ggplot(aes(y = mean, x = Treatment, fill = Treatment),  data = vac1) +
  geom_boxplot() + scale_fill_manual(values=c('#EBED9D','grey50','#868659','#CECE88')) +
  facet_wrap(~Species, labeller = labeller(Species = specnames)) + theme_classic() +
  xlab(element_blank()) + ylab("Symbiont:vacuole ratio") + theme(
    axis.text.x=element_blank(),
    strip.text.x = element_text(
      size = 11, face = "bold.italic",
      
    ))
dev.off()


#### exocytosis ####

# filtering data
summarized_data %>%
  filter(Measurement == "Exopercent") -> exo1

exo1$Species <- factor( as.character(exo1$Species), levels=c("Ofav","Mcav") )
exo1$Treatment <- factor( as.character(exo1$Treatment), levels=c("HS","DC","BDS","IDS") )

# beta regression for treatment and species comparisons
y <- exo1$mean
n.obs <- sum(!is.na(y))
y2 <- (y * (n.obs -1) + 0.5)/n.obs
exo1$y2 <- y2

beta_exo <- betareg(y2 ~ Treatment + Species, exo1)
joint_tests(beta_exo)

# individual species tests
exo_ofav <- betareg(y2 ~ Treatment, exo1[which(exo1$Species == 'Ofav'),])
summary(exo_ofav)
exo_mcav <- betareg(y2 ~ Treatment, exo1[which(exo1$Species == 'Mcav'),])
summary(exo_mcav)

# boxplot of exocytosis
specnames <- c("Orbicella faveolata", "Montastraea cavernosa")
names(specnames) <- c("Ofav", "Mcav")

pdf("histology_exocytosis_boxplot.pdf", width=6, height=4)
ggplot(aes(y = mean, x = Treatment, fill = Treatment),  data = exo1) +
  geom_boxplot() + scale_fill_manual(values=c('#EBED9D','grey50','#868659','#CECE88')) +
  facet_wrap(~Species, labeller = labeller(Species = specnames)) + theme_classic() +
  xlab(" ") + ylab("Proportion of Exocytosis") + theme(
    strip.text.x = element_text(
      size = 11, face = "bold.italic"
    ))
dev.off()


#### gastrodermal separation ####

# data filtering
summarized_data %>%
  filter(Measurement == 'Gastrosep') -> gastro1

gastro1$Species <- factor( as.character(gastro1$Species), levels=c("Ofav","Mcav") )
gastro1$Treatment <- factor( as.character(gastro1$Treatment), levels=c("HS","DC","BDS","IDS") )

# ANOVA for treatment and species comparisons
anova_gastro <- aov(mean ~ Treatment*Species, gastro1)
summary(anova_gastro)

# boxplot of gastrodermal separation
specnames <- c("Orbicella faveolata", "Montastraea cavernosa")
names(specnames) <- c("Ofav", "Mcav")

pdf("histology_gastrosep_boxplot.pdf", width=6, height=4)
ggplot(aes(y = mean, x = Treatment, fill = Treatment),  data = gastro1) +
  geom_boxplot() + scale_fill_manual(values=c('#EBED9D','grey50','#868659','#CECE88')) +
  facet_wrap(~Species, labeller = labeller(Species = specnames)) + theme_classic() +
  xlab(" ") + ylab("Mean gastrodermal separation (microns)") + theme(
    strip.text.x = element_text(
      size = 11, face = "bold.italic"
    ))
dev.off()


#### multiplot ####

vac1$Species <- factor( as.character(vac1$Species), levels=c("Ofav","Mcav") )
vac1$Treatment <- factor( as.character(vac1$Treatment), levels=c("HS","DC","BDS","IDS") )

exo1$Species <- factor( as.character(exo1$Species), levels=c("Ofav","Mcav") )
exo1$Treatment <- factor( as.character(exo1$Treatment), levels=c("HS","DC","BDS","IDS") )

gastro1$Species <- factor( as.character(gastro1$Species), levels=c("Ofav","Mcav") )
gastro1$Treatment <- factor( as.character(gastro1$Treatment), levels=c("HS","DC","BDS","IDS") )

specnames <- c("Orbicella faveolata", "Montastraea cavernosa")
names(specnames) <- c("Ofav", "Mcav")

vacplot <- ggplot(aes(y = mean, x = Treatment, fill = Treatment),  data = vac1) +
  geom_boxplot() + scale_fill_manual(values=c('#EBED9D','grey50','#868659','#CECE88')) +
  facet_wrap(~Species, labeller = labeller(Species = specnames)) + theme_classic() +
  xlab(element_blank()) + ylab("Symbiont:vacuole ratio") + theme(
    axis.text.x=element_blank(),
    strip.text.x = element_text(
      size = 11, face = "bold.italic"
    ))

exoplot <- ggplot(aes(y = mean, x = Treatment, fill = Treatment),  data = exo1) +
  geom_boxplot() + scale_fill_manual(values=c('#EBED9D','grey50','#868659','#CECE88')) +
  facet_wrap(~Species, labeller = labeller(Species = specnames)) + theme_classic() +
  xlab(" ") + ylab("Proportion of exocytosis") + theme(
    axis.text.x=element_blank(),
    strip.text.x = element_blank()
    )

gastroplot <- ggplot(aes(y = mean, x = Treatment, fill = Treatment),  data = gastro1) +
  geom_boxplot() + scale_fill_manual(values=c('#EBED9D','grey50','#868659','#CECE88')) +
  facet_wrap(~Species, labeller = labeller(Species = specnames)) + theme_classic() +
  xlab(" ") + ylab("Gastrodermal separation (μm)") + theme(
    axis.text.x=element_blank(),
    strip.text.x = element_blank()
    )


plot=grid.arrange(vacplot, exoplot, gastroplot, ncol=1, nrow=3, widths=c(6), heights=c(4,4,4))

ggsave("histology_boxplots.pdf", plot= plot, width=6, height=12, units="in", dpi=300)