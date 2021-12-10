#### packages ####

# BiocManager::install("phyloseq")
# install.packages("corncob")

library(corncob)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(reshape2)
library(cowplot)
library(grid)
library(scales)
library(RColorBrewer)
library(vegan)
library(ape)


#### data import ####

asv <- read.table("sctld sediment microbial ASV.txt", sep = "\t", header = TRUE, row.names = 1)
taxa <- as.matrix(read.table("sctld sediment microbial taxa.txt", sep = "\t", header = TRUE, row.names = 1))
metadata <- read.table("sctld sediment microbial metadata.txt", header = TRUE,row.names = 'sample_id')
asv <- t(asv)

# order the metadata so it matches the ASV table
idx <- match(rownames(asv),rownames(metadata))
metadata <- metadata[idx,]

metadata$condition <- factor( as.character(metadata$condition), levels=c("C","NAI","TL") )
metadata$treatment <- factor( as.character(metadata$treatment), levels=c("HS","BDS","IDS") )

ASV = otu_table(asv, taxa_are_rows = FALSE)
TAX = tax_table(taxa)
META = sample_data(metadata)

ps <- phyloseq(ASV, TAX, META)

length(get_taxa_unique(ps, "Family")) 
# number of unique family samples = 313
 
length(get_taxa_unique(ps, "Order")) 
# number of unique order samples = 213

get_taxa_unique(ps, "Phylum") 


#### differential abundance ####

# corncob differential abundance test 
set.seed(1)

# by sample condition
da_condition <- differentialTest(formula = ~ condition,
	phi.formula = ~ condition,
	formula_null = ~ 1,
	phi.formula_null = ~ condition,
	test = "Wald", boot = FALSE,
	data = ps,
	fdr_cutoff = 0.05)

# test statistics (number of significant taxa and test outputs for each DA taxa)
summary(da_condition$significant_taxa)
da_condition$significant_models

# lists the significant ASVs
da_condition$significant_taxa 

# plot of differentially-abundant ASVs in NAI and TL samples compared to control
pdf("sctld sediment microbial DA condition.pdf", width=12, height=4)
plot(da_condition)
dev.off()

# create genus-level-table
ps_genus <- ps %>% tax_glom("Genus")

da_condition_genus <- differentialTest(formula = ~ condition,
	phi.formula = ~ condition,
	formula_null = ~ 1,
	phi.formula_null = ~ condition,
	test = "Wald", boot = FALSE,
	data = ps_genus,
	fdr_cutoff = 0.05)

# test statistics (number of significant taxa and test outputs for each DA taxa)
summary(da_condition_genus$significant_taxa)
da_condition_genus$significant_models

# lists the significant ASVs
da_condition_genus$significant_taxa 

# plot of differentially-abundant ASVs in NAI and TL samples compared to control
pdf("sctld sediment microbial DA condition genus.pdf", width=10, height=6)
plot(da_condition_genus)
dev.off()

# by treatment
da_treatment <- differentialTest(
  formula = ~ treatment,
  phi.formula = ~ treatment,
  formula_null = ~ 1,
  phi.formula_null = ~ treatment,
  test = "Wald",
  boot = FALSE,
  data = ps,
  fdr_cutoff = 0.05
)

# test statistics (number of significant taxa and test outputs for each DA taxa)
summary(da_treatment$significant_taxa)
da_treatment$significant_models

# lists the significant ASVs
da_treatment$significant_taxa 

# plot of differentially-abundant ASVs in treatments BDS and IDS compared to HS
pdf("sctld sediment microbial DA treatment.pdf", width=16, height=40)
plot(da_treatment)
dev.off()

da_treatment_genus <- differentialTest(
  formula = ~ treatment,
  phi.formula = ~ treatment,
  formula_null = ~ 1,
  phi.formula_null = ~ treatment,
  test = "Wald",
  boot = FALSE,
  data = ps_genus,
  fdr_cutoff = 0.05
)

# test statistics (number of significant taxa and test outputs for each DA taxa)
summary(da_treatment_genus$significant_taxa)
da_treatment_genus$significant_models

# lists the significant ASVs
da_treatment_genus$significant_taxa 

# plot of differentially-abundant ASVs in treatments BDS and IDS compared to HS
pdf("sctld sediment microbial DA treatment genus.pdf", width=12, height=18)
plot(da_treatment_genus)
dev.off()

# filtering dataset to remove BDS (i.e., look at ASVs related to SCTLD exposure while eliminating sediment origin effects; HS vs IDS)
ps_sub <- subset_samples(ps, treatment != 'BDS')

# by treatment
da_treatment_sub <- differentialTest(
  formula = ~ treatment,
  phi.formula = ~ treatment,
  formula_null = ~ 1,
  phi.formula_null = ~ treatment,
  test = "Wald",
  boot = FALSE,
  data = ps_sub,
  fdr_cutoff = 0.05
)

# test statistics (number of significant taxa and test outputs for each DA taxa)
summary(da_treatment_sub$significant_taxa)
da_treatment_sub$significant_models

# lists the significant ASVs
da_treatment_sub$significant_taxa 

# plot of differentially-abundant ASVs in treatment IDS compared to HS
pdf("sctld sediment microbial DA treatment subset.pdf", width=16, height=2)
plot(da_treatment_sub)
dev.off()

# create genus-level-table
ps_genus_sub <- ps_sub %>% tax_glom("Genus")

da_treatment_genus_sub <- differentialTest(
  formula = ~ treatment,
  phi.formula = ~ treatment,
  formula_null = ~ 1,
  phi.formula_null = ~ treatment,
  test = "Wald",
  boot = FALSE,
  data = ps_genus_sub,
  fdr_cutoff = 0.05
)

# test statistics (number of significant taxa and test outputs for each DA taxa)
summary(da_treatment_genus_sub$significant_taxa)
da_treatment_genus_sub$significant_models

# lists the significant ASVs
da_treatment_genus_sub$significant_taxa 

# plot of differentially-abundant ASVs in treatment IDS compared to HS
pdf("sctld sediment microbial DA treatment genus subset.pdf", width=16, height=2)
plot(da_treatment_genus_sub)
dev.off()

save(ps,ps_genus,da_condition,da_condition_genus,da_treatment,da_treatment_genus,ps_sub,ps_genus_sub,da_treatment_sub,da_treatment_genus_sub,file="sctld sediment microbial corncob.RData")
load("sctld sediment microbial corncob.RData")


#### PCoA ####

diss <- vegdist(asv, "bray")
pcoa <- pcoa(diss)
scores <- pcoa$values
vectors <- pcoa$vectors

pdf(file="SCTLD sediment microbial PCoA.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(vectors[,1], vectors[,2],col=c('#EBED9D','#868659','#CECE88')[as.numeric(as.factor(metadata$treatment))],pch=c(15,1,19)[as.numeric((as.factor(metadata$condition)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
# ordispider(vectors, metadata$treatment, label=F, col=c('#EBED9D','#868659','#CECE88'))
legend("top", legend=c("HS", "BDS","IDS"), fill = c('#EBED9D','#868659','#CECE88'), bty="n")
legend("bottom", legend=c("control","NAI","TL"), pch=c(15,1,19), bty="n")
plot(vectors[,1], vectors[,2],col=c("green","orange","red")[as.numeric(as.factor(metadata$condition))],pch=c(15,17,25)[as.numeric(as.factor(metadata$treatment))], xlab="Coordinate 1", ylab="Coordinate 2", main="Condition")
# ordispider(vectors, metadata$condition, label=F, col=c("green","orange","red"))
legend("top", legend=c("control", "NAI", "TL"), fill = c("green","orange","red"), bty="n")
legend("bottom", legend=c("HS","BDS","IDS"), pch=c(15,17,25), bty="n")
dev.off()

# filtering dataset to remove BDS
metadata_sub <- subset(metadata, treatment !='BDS')
asv_sub <- subset(asv, row.names(asv) %in% row.names(metadata_sub))

metadata_sub$condition <- factor( as.character(metadata_sub$condition), levels=c("C","NAI","TL") )
metadata_sub$treatment <- factor( as.character(metadata_sub$treatment), levels=c("HS","IDS") )

diss_sub <- vegdist(asv_sub, "bray")
pcoa_sub <- pcoa(diss_sub)
scores_sub <- pcoa_sub$values
vectors_sub <- pcoa_sub$vectors

pdf(file="SCTLD sediment microbial PCoA subset.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(vectors_sub[,1], vectors_sub[,2],col=c('#EBED9D','#CECE88')[as.numeric(as.factor(metadata_sub$treatment))],pch=c(15,1,19)[as.numeric((as.factor(metadata_sub$condition)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
legend("top", legend=c("HS","IDS"), fill = c('#EBED9D','#CECE88'), bty="n")
legend("bottom", legend=c("control","NAI","TL"), pch=c(15,1,19), bty="n")
plot(vectors_sub[,1], vectors_sub[,2],col=c("green","orange","red")[as.numeric(as.factor(metadata_sub$condition))],pch=c(15,25)[as.numeric(as.factor(metadata_sub$treatment))], xlab="Coordinate 1", ylab="Coordinate 2", main="Condition")
legend("top", legend=c("control", "NAI", "TL"), fill = c("green","orange","red"), bty="n")
legend("bottom", legend=c("HS","IDS"), pch=c(15,25), bty="n")
dev.off()


#### PERMANOVA ####
perm_condition <- adonis(diss ~ species*condition, metadata)
perm_condition

pair_condition <- pairwise.adonis(diss, metadata$condition, p.adjust.m='fdr')
pair_condition

perm_treatment <- adonis(diss ~ species*treatment, metadata)
perm_treatment

pair_treatment <- pairwise.adonis(diss, metadata$treatment, p.adjust.m='fdr')
pair_treatment

# filtered dataset to remove BDS
perm_condition_sub <- adonis(diss_sub ~ species*condition, metadata_sub)
perm_condition_sub

perm_treatment_sub <- adonis(diss_sub ~ species*treatment, metadata_sub)
perm_treatment_sub


#### bubble plot ####

sign.genera <- read.table("sctld sediment microbial DA genera.txt", sep = "\t", header = TRUE)

# convert data frame from a "wide" format to a "long" format
sign.genera.long = melt(sign.genera, id = c("Sample"))
sign.genera.long <- sign.genera.long %>% separate(Sample, 
                c("Sample", "Treatment","Condition"))

# pick colors (you need 16 distinct colors) 
colors = c("#833FE2", "#DD579C", "#CEA984", "#8DAAE4", 
			"#8ADF97", "#DF7573", "#836678", "#E4CED1", 
			"#7E74D5", "#73E2CF", "#E8924A", "#CBE545", 
			"#70B8C6", "#DBA9D3", "#C7E5E6", "#72E356")
						
show_col(colors)

# keep the order of samples from your data
sign.genera.long$Condition <- factor(sign.genera.long$Condition,levels=unique(sign.genera.long$Condition))
 
bubble <- ggplot(sign.genera.long, aes(x = Sample, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0, 7), range = c(0,17), breaks = c(0.5,2.5,5)) + 
  facet_grid(. ~ Condition + Treatment, scales="free", space="free") +
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
  axis.text.x = element_text(color = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
  axis.text.y = element_text(color = "black", face = "bold", size = 12), 
  legend.text = element_text(size = 10, face ="bold", color ="black"), 
  legend.title = element_text(size = 12, face = "bold"), 
  panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 1.2), 
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  legend.position = "top",legend.justification="left", legend.direction = "horizontal") +  
  scale_fill_manual(values = colors, guide = "none") + 
  scale_y_discrete(limits = rev(levels(sign.genera.long$variable))) 
bubble

ggsave("SCTLD sediment microbial abundance.pdf", plot= bubble, width=10, height=8, units="in", dpi=300)