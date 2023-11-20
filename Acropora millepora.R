library(plotrix)
library(microbiome)
library(gridExtra)
library(nlme)
library(emmeans)
library(dplyr)
library(ggplot2)
library(Rmisc)

##############################################################################
#####                       ALPHA DIVERSITY                                ###
##############################################################################

sort(sample_sums(AM)) 

# rarefy
AM_rare<- rarefy_even_depth(AM, sample.size = 2184, rngseed = 1)

# extract metadata from subsetted file
divMeta <- meta(AM_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(AM_rare))
#write.csv(divMeta, "AM_alphadiv.csv")

#Observed ASVs_ANOVA
Obs <- lme(Observed ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) #no significant differences
summarySE(divMeta, measurevar = "Observed")

#Simpsons_ANOVA
Sim <- lme(InvSimpson ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Sim)
anova(Sim)
summarySE(divMeta, groupvars = "PMA", measurevar = "InvSimpson")


#Shannon_ANOVA
Shan <- lme(Shannon ~ PMA * Genotype, random = ~1|Tank,
            data = divMeta, na.action = na.omit)
summary(Shan)
anova(Shan)

#              numDF denDF  F-value p-value
#Genotype         3    34  6.438680  0.0014

pairs(emmeans(Shan, "Genotype"), adjust="bonferroni")

#contrast      estimate    SE df t.ratio p.value
#Black - Blue    0.6243 0.275 34   2.268  0.1788
#Black - Cream  -0.4352 0.268 34  -1.622  0.6846
#Black - Green  -0.3424 0.268 34  -1.276  1.0000
#Blue - Cream   -1.0595 0.268 34  -3.958  0.0022
#Blue - Green   -0.9667 0.268 34  -3.611  0.0058
#Cream - Green   0.0928 0.261 34   0.356  1.0000

#reorder levels of PMA
divMeta$PMA <- factor(divMeta$PMA, levels = c("Untreated","PMA-treated"))

##BoxPlot###

AM.obs <- ggplot(divMeta, aes(x=Genotype, y=Observed, fill=Genotype), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#FFFFFF", "#4e79a7", "#F5F5DC", "#59a14f")) + lims(y=c(0,150))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("A. Observed ASVs")
print(AM.obs)

AM.sim <- ggplot(divMeta, aes(x=Genotype, y=InvSimpson, fill=Genotype)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#FFFFFF", "#4e79a7", "#F5F5DC", "#59a14f"))+ lims(y=c(0,4.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("B. inverse Simpson")
print(AM.sim)

AM.shan <- ggplot(divMeta, aes(x=Genotype, y=Shannon, fill=Genotype)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#FFFFFF", "#4e79a7", "#F5F5DC", "#59a14f"))+ lims(y=c(0,3.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("C. Shannon")
print(AM.shan)

grid.arrange(AM.obs,AM.sim,AM.shan, nrow=1, ncol=3)


rm(AM_rare)
rm(AM.obs)
rm(AM.shan)
rm(AM.sim)
rm(Obs)
rm(Shan)
rm(Sim)
rm(divMeta)


##############################################################################
#####                        BETA DIVERSITY                                ###
##############################################################################

##############################
##          PERMANOVA       ##
##############################

library(vegan)
library(pairwiseAdonis)
library(microbiome)

comp <- microbiome::transform(AM, transform = "compositional")

#Generate unweighted Unifrac distance matrix
wunifrac_dist_matrix <- phyloseq::distance(comp, method = "wunifrac")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#weighted Unifrac
dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$PMA)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.5208
#fail to reject the assumption of homogeneity of dispersion by PMA treatment

dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.2677
#fail to reject the assumption of homogeneity of dispersion by colony/genotype

#Because the weighted unifrac matrix met the assumption of homogeneity for PMA, 
#we continue with that distance matrix

#ADONIS test; strata = sets random effect
vegan::adonis2(wunifrac_dist_matrix ~ phyloseq::sample_data(comp)$PMA*phyloseq::sample_data(comp)$Genotype, strata = phyloseq::sample_data(comp)$Tank, permutations = 9999) 

#                                                                      Df SumOfSqs      R2      F Pr(>F)    
#phyloseq::sample_data(comp)$PMA                                       1  0.01917 0.02576 1.3830 0.1986
#phyloseq::sample_data(comp)$Genotype                                  3  0.13175 0.17702 3.1675 0.1458
#phyloseq::sample_data(comp)$PMA:phyloseq::sample_data(comp)$Genotype  3  0.02490 0.03345 0.5986 0.6326
#Residual                                                             41  0.56847 0.76377              
#Total                                                                48  0.74429 1.00000                                                                  47  0.176569  1.00000                                                              37  1.03150 1.00000                                                              38  12.6609 1.00000 

rm(dispr.wunifrac)
rm(wunifrac_dist_matrix)

##############################################################################
#####                               PCoA                                   ###
##############################################################################


wunifrac <- ordinate(comp, method = "PCoA", distance = "wunifrac")
W <- plot_ordination(physeq = comp,
                     ordination = wunifrac,
                     type = "samples",
                     color = "PMA",
                     shape = "PMA") + 
  theme_bw() + scale_color_manual(values=c("#6A6AAD","#ADAD6A")) +
  geom_point(size=3) 

W + stat_ellipse(type = "t", linetype = 2, level = 0.70) + ggtitle("A. Weighted Unifrac PCoA")

rm(W)
rm(wunifrac)
rm(comp)


##############################################################################
#####                             BARPLOT                                  ###
############################################################################## 

#Merge by groups for figure

all <- merge_samples(AM, "Genotype")
all <- prune_taxa((taxa_sums(all) > 0), all)


##Genus level##
# To represent at the Genus level (instead of ASV level)
Genus <- tax_glom(all,taxrank = "Genus")

# Instead of making copies of this code - I just run the same code and change the phyloseq object that is being transformed
#I also manually reorganized the colours so each group (all, gender, stage) would match.
# Transform counts in relative abundance and select most abundant families
Genus <- transform_sample_counts(Genus, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$Genus <- as.character(genus$Genus)

#rename Genera with <2% abundance
genus$Genus[genus$Abundance < 2]<- "< 2% Abundance"
write.csv(genus, "genus_AM.csv", row.names=FALSE)
genus <- read.csv("genus_AM.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))


BP_AM <- ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(Genus, Abundance))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = c("#9A879D","#E56B70","#1B998B","#f7fff7","#0B5351","#92140C","#7D84B2","#FB9A99","#E3D8F1","#48304D",
                               "#A6CEE3","#6A3D9A","#1F78B4","#FF7F00","#CAB2D6","#33A02C","#1F3172", "#E7298A","#FDBF6F")) + ggtitle("Acropora millepora") + coord_flip()

print(BP_AM)
rm(Genus)
rm(genus)
rm(HowMany)
rm(all)
rm(BP_AM)

#Endozoicomonas ASVs by Genotype 

merge <- merge_samples(AM, "Genotype")
comp <- transform_sample_counts(merge, function(x) 100 * x/sum(x))
Endo <- subset_taxa(comp, Genus=="Endozoicomonas")
Endo <- prune_taxa((taxa_sums(Endo) > 0), Endo)
endo <- psmelt(Endo)
endo$OTU <- as.character(endo$OTU)
write.csv(endo, "Endo_AM.csv", row.names=FALSE)
endo <- read.csv("Endo_AM.csv")

#How many levels in endo
HowMany <- length(levels(as.factor(endo$OTU)))


#reorder levels of Genotype
endo$Sample <- factor(endo$Sample, levels = c("Blue","Green","Cream","Black"))



plot_AM <- ggplot(endo, aes(x = Sample, y = Abundance, fill = reorder(OTU, Abundance))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance (%) of Endozoicomonas \n") + scale_fill_manual(values = c("#9A879D","#E56B70","#343434","#1B998B","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#FB9A99","#E3D8F1","#48304D",
                                                                                     "#A6CEE3","#6A3D9A","#1F78B4","#FF7F00","#CAB2D6","#33A02C","#1F3172", "#E7298A","#FDBF6F")) +
  ggtitle("Acropora millepora") + coord_flip()

print(plot_AM)

rm(merge)
rm(Endo)
rm(endo)
rm(comp)
rm(HowMany)
rm(plot_AM)


########################################
### Shared Core Taxa - Venn Diagram  ###
########################################

#load packages
library(eulerr)
library(microbiome)

## Shared taxa between PMAs
#simple way to count the number of samples in each group
table(meta(AM)$PMA, useNA = "always")

#convert to relative abundance
AM.rel <- microbiome::transform(AM, "compositional")

#make a list of PMA types
PMA.all <- unique(as.character(meta(AM.rel)$PMA))

#Write a for loop to go through each of the PMAs one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in PMA.all){ # for each variable n in PMA.all
  #print(paste0("Identifying Core Taxa for ", n))
  
  AM.sub <- subset_samples(AM.rel, PMA.all == n) # Choose sample from PMA.all by n
  
  core_m <- core_members(AM.sub, # AM.sub is phyloseq object selected with only samples from n PMA type
                         detection = 1/100, # 1% in any sample 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each substrate.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

#print(list_core)

mycols <- c("#6A6AAD","#ADAD6A") 
plot(venn(list_core), fills = mycols)

