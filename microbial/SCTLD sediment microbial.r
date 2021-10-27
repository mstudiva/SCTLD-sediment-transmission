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
library(reshape2)
library(cowplot)
library(grid)
library(scales)
library(RColorBrewer)
library(vegan)
library(ape)


#### data import ####

asv <- read.table("ASV_table.txt", sep = "\t", header = TRUE, row.names = 1)
taxa <- as.matrix(read.table("taxa.txt", sep = "\t", header = TRUE, row.names = 1))
metadata <- read.table("metadata.txt", header = TRUE,row.names = 'sample_id')
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


#### statistical tests ####

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

# lists the significant ASVs
da_condition$significant_taxa 

# plot of differentially-abundant ASVs in NAI and TL samples compared to control
pdf("microbial_DA_condition.pdf", width=12, height=4)
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

# lists the significant ASVs
da_condition_genus$significant_taxa 

# plot of differentially-abundant ASVs in NAI and TL samples compared to control
pdf("microbial_DA_condition_genus.pdf", width=10, height=6)
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

# lists the significant ASVs
da_treatment$significant_taxa 

# plot of differentially-abundant ASVs in treatments BDS and IDS compared to HS
pdf("microbial_DA_treatment.pdf", width=16, height=40)
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

# lists the significant ASVs
da_treatment_genus$significant_taxa 

# plot of differentially-abundant ASVs in treatments BDS and IDS compared to HS
pdf("microbial_DA_treatment_genus.pdf", width=12, height=18)
plot(da_treatment_genus)
dev.off()

save(da_condition,da_condition_genus,da_treatment,da_treatment_genus,file="microbial_corncob.RData")

#### PCoA ####

diss <- vegdist(asv, "bray")
pcoa <- pcoa(diss)
scores <- pcoa$values
vectors <- pcoa$vectors

pdf(file="SCTLD_sediment_microbial_PCoA.pdf", width=12, height=6)
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


#### bubble plot ####

sign.genera <- read.table("sign.genera.txt", sep = "\t", header = TRUE)

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

ggsave("SCTLD_sediment_microbial_abundance.pdf", plot= bubble, width=10, height=8, units="in", dpi=300)