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


#### data import ####

asv <- read.table("ASV_table.txt", sep = "\t", header = TRUE, row.names = 1)
taxa <- as.matrix(read.table("taxa.txt", sep = "\t", header = TRUE, row.names = 1))
metadata <- read.table("metadata.txt", header = TRUE,row.names = 'sample_id')
asv <- t(asv)

# order the metadata so it matches the ASV table
idx <- match(rownames(asv),rownames(metadata))
metadata <- metadata[idx,]

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

# corncob differencial abundance test 
set.seed(1)

da_analysis <- differentialTest(formula = ~ condition,
	phi.formula = ~ condition,
	formula_null = ~ 1,
	phi.formula_null = ~ condition,
	test = "Wald", boot = FALSE,
	data = ps,
	fdr_cutoff = 0.05)

da_analysis
da_analysis$significant_taxa 
plot(da_analysis)

# create genus-level-table
ps_genus <- ps %>% tax_glom("Genus")

da_genus <- differentialTest(formula = ~ condition,
	phi.formula = ~ condition,
	formula_null = ~ 1,
	phi.formula_null = ~ condition,
	test = "Wald", boot = FALSE,
	data = ps_genus,
	fdr_cutoff = 0.05)

plot(da_genus)


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
 
xx = ggplot(sign.genera.long, aes(x = Sample, y = variable)) + 
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
xx

ggsave("microbial_abundance.pdf", plot= xx, width=10, height=8, units="in", dpi=300)