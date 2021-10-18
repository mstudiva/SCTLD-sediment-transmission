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
library(ggplot)


asv <- read.table("ASV_table.txt", sep = "\t", header = TRUE, row.names = 1)
taxa <- as.matrix(read.table("taxa.txt", sep = "\t", header = TRUE, row.names = 1))
metadata <- read.table("metadata.txt", header = TRUE,row.names = 'sample_id')
asv <- t(asv)

#order the metadata so it matches the ASV table

idx <- match(rownames(asv),rownames(metadata))

metadata <- metadata[idx,]

ASV = otu_table(asv, taxa_are_rows = FALSE)
TAX = tax_table(taxa)
META = sample_data(metadata)

ps <- phyloseq(ASV, TAX, META)

length(get_taxa_unique(ps, "Family")) 

#number of unique family samples = 313
 
length(get_taxa_unique(ps, "Order")) 

#number of unique order samples = 213

get_taxa_unique(ps, "Phylum") 

#corncob-different test 

set.seed(1)

da_analysis <- differentialTest(formula = ~ fate,
	phi.formula = ~ fate,
	formula_null = ~ 1,
	phi.formula_null = ~ fate,
	test = "Wald", boot = FALSE,
	data = ps,
	fdr_cutoff = 0.05)

da_analysis
da_analysis$significant_taxa 
plot(da_analysis)

#to create genus-level-table

ps_genus <- ps %>% tax_glom("Genus")

da_genus <- differentialTest(formula = ~ fate,
	phi.formula = ~ fate,
	formula_null = ~ 1,
	phi.formula_null = ~ fate,
	test = "Wald", boot = FALSE,
	data = ps_genus,
	fdr_cutoff = 0.05)

plot(da_genus)


#to generate Figure 5 

sign.genera <- read.table("sign.genera.txt", sep = "\t", header = TRUE)

#convert data frame from a "wide" format to a "long" format

sign.genera.long = melt(sign.genera, id = c("Sample"))

#pick colors (you need 16 distinct colors) 

colours = c("#833FE2", "#DD579C", "#CEA984", "#8DAAE4", 
			"#8ADF97", "#DF7573", "#836678", "#E4CED1", 
			"#7E74D5", "#73E2CF", "#E8924A", "#CBE545", 
			"#70B8C6", "#DBA9D3", "#C7E5E6", "#72E356")
						
show_col(colours)

#keep the order of samples from your data

sign.genera.long$Sample <- factor(sign.genera.long$Sample,levels=unique(sign.genera.long$Sample))
 
xx = ggplot(sign.genera.long, aes(x = Sample, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0, 7), range = c(0,17), breaks = c(0.5,5)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
  axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
  axis.text.y = element_text(colour = "black", face = "bold", size = 12), 
  legend.text = element_text(size = 10, face ="bold", colour ="black"), 
  legend.title = element_text(size = 12, face = "bold"), 
  panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
  legend.position = "top",legend.justification="left", legend.direction = "horizontal") +  
  scale_fill_manual(values = colours, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(sign.genera.long$variable))) 
  
xx

