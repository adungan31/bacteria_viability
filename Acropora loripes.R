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

sort(sample_sums(AL)) 

# rarefy
AL_rare<- rarefy_even_depth(AL, sample.size = 737, rngseed = 1)

# extract metadata from subsetted file
divMeta <- meta(AL_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(AL_rare))
write.csv(divMeta, "AL_alphadiv.csv")

#Observed ASVs_ANOVA
Obs <- lme(Observed ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) #no significant differences
#summary data for observed ASVs
summarySE(divMeta, measurevar = "Observed")

#Simpsons_ANOVA
Sim <- lme(InvSimpson ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Sim)
anova(Sim)

#              numDF denDF   F-value p-value
#PMA             1    35 12.519971  0.0012


#Shannon_ANOVA
Shan <- lme(Shannon ~ PMA * Genotype, random = ~1|Tank,
            data = divMeta, na.action = na.omit)
summary(Shan)
anova(Shan)

#              numDF denDF  F-value p-value
#PMA              1    35  15.75356  0.0003


#reorder levels of PMA
divMeta$PMA <- factor(divMeta$PMA, levels = c("Untreated","PMA-treated"))

##BoxPlot###

AL.obs <- ggplot(divMeta, aes(x=PMA, y=Observed, fill=PMA), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD")) + lims(y=c(0,60))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("A. Observed ASVs")
print(AL.obs)

AL.sim <- ggplot(divMeta, aes(x=PMA, y=InvSimpson, fill=PMA)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD"))+ lims(y=c(0,18))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("B. Inverse Simpson")
print(AL.sim)

AL.shan <- ggplot(divMeta, aes(x=PMA, y=Shannon, fill=PMA)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD"))+ lims(y=c(0,4))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("C. Shannon")
print(AL.shan)


grid.arrange(AL.obs,AL.sim,AL.shan, nrow=1, ncol=3)


rm(AL_rare)
rm(AL.obs)
rm(AL.shan)
rm(AL.sim)
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

comp <- microbiome::transform(AL, transform = "compositional")

#Generate weighted Unifrac distance matrix
unifrac_dist_matrix <- phyloseq::distance(comp, method = "wunifrac")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#Weighted Unifrac
dispr.unifrac <- vegan::betadisper(unifrac_dist_matrix, phyloseq::sample_data(comp)$PMA)
dispr.unifrac
plot(dispr.unifrac)
anova(dispr.unifrac)#p=0.5852
#fail to reject the assumption of homogeneity of dispersion by PMA treatment

dispr.unifrac <- vegan::betadisper(unifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype)
dispr.unifrac
plot(dispr.unifrac)
anova(dispr.unifrac)#p=0.2777
#fail to reject the assumption of homogeneity of dispersion by colony/genotype

#Because the unifrac matrix failed to meet the assumption of homogeneity for PMA and genotype we 
#continue with that distance matrix in the adonis test but add in a permutation test for multivariate 
#dispersion to account for unequal dispersion between groups

#ADONIS test; strata = sets random effect
vegan::adonis2(unifrac_dist_matrix ~ phyloseq::sample_data(comp)$PMA*phyloseq::sample_data(comp)$Genotype, strata = phyloseq::sample_data(comp)$Tank, permutations = 9999) 

#                                                                      Df SumOfSqs      R2      F Pr(>F)    
#phyloseq::sample_data(comp)$PMA                                       1  0.17107 0.34604 24.1511 0.0001 ***
#phyloseq::sample_data(comp)$Genotype                                  4  0.04389 0.08877  1.5489 0.1352    
#phyloseq::sample_data(comp)$PMA:phyloseq::sample_data(comp)$Genotype  4  0.02441 0.04937  0.8614 0.5610    
#Residual                                                             36  0.25500 0.51582                   
#Total                                                                45  0.49436 1.00000                                                                   47  0.176569  1.00000                                                              37  1.03150 1.00000                                                              38  12.6609 1.00000 


rm(dispr.unifrac)
rm(unifrac_dist_matrix)

##############################################################################
#####                               PCoA                                   ###
##############################################################################

unifrac <- ordinate(comp, method = "PCoA", distance = "wunifrac")
W <- plot_ordination(physeq = comp,
                     ordination = unifrac,
                     type = "samples",
                     color = "PMA",
                     shape = "PMA") + 
  theme_bw() + scale_color_manual(values=c("#6A6AAD","#ADAD6A")) +
  geom_point(size=3) 

W + ggtitle("A. Weighted Unifrac PCoA - A. loripes") + stat_ellipse(type = "t", linetype = 2, level = 0.70) 

rm(W)
rm(unifrac)
rm(comp)


##############################################################################
#####                             BARPLOT                                  ###
############################################################################## 

#Merge by groups for figure

all <- merge_samples(AL, "PMA")
all <- prune_taxa((taxa_sums(all) > 0), all)


##Genus level##
# To represent at the Genus level (instead of ASV level)
Genus <- tax_glom(all,taxrank = "Genus")

# Instead of making copies of this code - I just run the same code and change the phyloseq object that is being transformed
#I also manually reorganized the colours so each group (all, PMA, stage) would match.
# Transform counts in relative abundance and select most abundant families
Genus <- transform_sample_counts(Genus, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$Genus <- as.character(genus$Genus)

#rename Genera with <1% abundance
genus$Genus[genus$Abundance < 1]<- "< 1% Abundance"
write.csv(genus, "genus_AL.csv", row.names=FALSE)
genus <- read.csv("genus_AL.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))


BP <- ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(Genus, Abundance))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = c("#9A879D","#E56B70","#A6CEE3","#48304D","#f7fff7","#0B5351",
                               "#92140C","#7D84B2","#4B7F52","#D2F9F9","#FB9A99","#E3D8F1")) + ggtitle("B. Genus level Barplot - A. loripes")

BP + coord_flip()
rm(Genus)
rm(genus)
rm(HowMany)
rm(all)
rm(BP)

##############################################################################
#####                        INDICATOR SPECIES                             ###
############################################################################## 

library(indicspecies)
library(plyr)
library(phylosmith)
library(microbiome)


###For this to work correctly, the order of your samples MUST MATCH the order of the classifications you assign.

##ASV LEVEL###

# This package does not work with data in Phyloseq format, so we will have to adjust that.
rel <- transform_sample_counts(AL, function(x) 100 * x/sum(x))

#order metadata by variable of interest and check metadata table to confirm
rel <- set_sample_order(rel, 'PMA')
meta <- microbiome::meta(rel)
count(meta$PMA)#use these numbers and order for the 'group' below (152 female, 26 male)

# Community data matrix: convert from phyloseq object
asv.table <- as(otu_table(rel), "matrix")

# Convert in order to have taxa in columns and species in rows:
asv.table <- t(asv.table)

# save row names to file
sampNames <- row.names(asv.table)

# remove rownames from otu table
row.names(asv.table) <- NULL

# convert otu table to data frame
asv.table <- as.data.frame(asv.table)

# add row names to data frame as column of data
asv.table <- cbind(sampNames, asv.table)


# Defining classification of samples (groups have to be in discrete categories)

group <- c(rep("PMA-treated", 23),rep("Untreated", 23))


# Performing the indicator value analysis (omit sampNames column in asv table, otherwise indval fails)
indval = multipatt(asv.table[,-1], group, control = how(nperm=9999))


# "Indval.g" allows to correct for different sample size per group
# duleg = TRUE would not allow to consider combination of groups; so we would obtain species associated either
# with the Low, Medium or High treatment. If we want combinations of groups to be taken into account,
# we can use the parameter "restcomb = c(...)" instead of "duleg = TRUE", where we list which groups are to be considered.
# The list goes in that way: 1 = Low, 2 = Medium, 3 = High, 4 = Low + Medium, 5 = Low + High, 6 = Medium + High
# So if we want to consider Medium + High in addition to them individually, we would type:
#indval2 = multipatt(asv.table, groups, control = how(nperm=999), restcomb = c(1, 2, 3, 6), func = "IndVal.g")

# Displaying the data: we can set the thresholds of specificity, sensitivity and significance
summary(indval, At = 0.6, Bt = 0.6, alpha = 0.05, func = "IndVal.g")


# As explained in the tutorial:
# The indicator value index is the product of two components, called "A" and "B"
# Component "A" is the probability that a sample belongs to its target group
# given the fact that the ASV has been found. This conditional probability is called the specificity
# of the ASV as indicator of the sample group. Component "B" is the probability of finding the ASV
# in samples belonging to the sample group. This second conditional probability is called the fidelity
# of the ASV as indicator of the target sample group.

# For easier reporting, transfer the output in a file:
capture.output(summary(indval, At = 0.6, Bt = 0.6, invdalcomp = TRUE, alpha = 0.05, func = 'Indval.g'), file = "IndvalOutput_PMA_AL.txt")


##Genus LEVEL###

Genus <- tax_glom(AL,taxrank = "Genus")

# This package does not work with data in Phyloseq format, so we will have to adjust that.
rel <- transform_sample_counts(Genus, function(x) 100 * x/sum(x))

#order metadata by variable of interest and check metadata table to confirm
rel <- set_sample_order(rel, 'PMA')
meta <- microbiome::meta(rel)
count(meta$PMA)#use these numbers and order for the 'group' below (152 female, 26 male)

# Community data matrix: convert from phyloseq object
asv.table <- as(otu_table(rel), "matrix")

# Convert in order to have taxa in columns and species in rows:
asv.table <- t(asv.table)

# save row names to file
sampNames <- row.names(asv.table)

# remove rownames from otu table
row.names(asv.table) <- NULL

# convert otu table to data frame
asv.table <- as.data.frame(asv.table)

# add row names to data frame as column of data
asv.table <- cbind(sampNames, asv.table)


# Defining classification of samples (groups have to be in discrete categories)

group <- c(rep("PMA-treated", 23),rep("Untreated", 23))


# Performing the indicator value analysis (omit sampNames column in asv table, otherwise indval fails)
indval = multipatt(asv.table[,-1], group, control = how(nperm=9999))

summary(indval, At = 0.6, Bt = 0.6, alpha = 0.05, func = "IndVal.g")

capture.output(summary(indval, At = 0.6, Bt = 0.6, invdalcomp = TRUE, alpha = 0.05, func = 'Indval.g'), file = "IndvalOutput_PMA_AL_Genus.txt")

rm(rel)
rm(Genus)
rm(asv.table)
rm(indval)
rm(meta)
rm(group)
rm(sampNames)


#Endozoicomonas ASVs by PMA

merge <- merge_samples(AL, "PMA")
comp <- transform_sample_counts(merge, function(x) 100 * x/sum(x))
Endo <- subset_taxa(comp, Genus=="Endozoicomonas")
Endo <- prune_taxa((taxa_sums(Endo) > 0), Endo)
endo <- psmelt(Endo)
endo$OTU <- as.character(endo$OTU)
write.csv(endo, "Endo.csv", row.names=FALSE)
endo <- read.csv("Endo.csv")

#How many levels in endo
HowMany <- length(levels(as.factor(endo$OTU)))

#reorder levels of PMA
endo$Sample <- factor(endo$Sample, levels = c("PMA-treated","Untreated"))

plot_AL <- ggplot(endo, aes(x = Sample, y = Abundance, fill = reorder(OTU, Abundance))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance (%) of Endozoicomonas \n") + scale_fill_manual(values = c("#9A879D","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#FB9A99","#E3D8F1","#E56B70","#343434","#1B998B","#48304D",
                                                                                     "#A6CEE3","#6A3D9A","#1F78B4","#FF7F00","#CAB2D6","#33A02C","#1F3172", "#E7298A","#FDBF6F")) +
  ggtitle("Acropora loripes") + coord_flip()

print(plot_AL) 

rm(merge)
rm(Endo)
rm(endo)
rm(comp)
rm(HowMany)
rm(plot_AL)

#create a list of values for indicator taxa (by gestation stage)
S <- c("99bca50b5dccc7d415db3473750009c7","fc9710608cae96aeaa60fc2d2882150c")

merge <- merge_samples(AL, "PMA")


# This package does not work with data in Phyloseq format, so we will have to adjust that.
rel <- transform_sample_counts(merge, function(x) 100 * x/sum(x))

# retain only the selected ASVs
phyS <- prune_taxa(S, rel)

# export otu table as data frame
phyS.df <- as.data.frame(phyS@otu_table)
phyS.df <- t(phyS.df)

# reshape data
library(reshape2)
data <- data.frame(Family = paste(phyS@tax_table[,3],phyS@tax_table[,4],phyS@tax_table[,5],phyS@tax_table[,6],sep = "_"), phyS.df)
write.csv(data,"IndicSpecies_IDs_Alor.csv", row.names = TRUE)



#create a list of values for indicator taxa (by gestation stage)
S <- "9f6605fed5b773446f8f5c2233a7a823"

merge <- merge_samples(PA, "PMA")


# This package does not work with data in Phyloseq format, so we will have to adjust that.
rel <- transform_sample_counts(merge, function(x) 100 * x/sum(x))

# retain only the selected ASVs
phyS <- prune_taxa(S, rel)

# export otu table as data frame
phyS.df <- as.data.frame(phyS@otu_table)
phyS.df <- t(phyS.df)

# reshape data
library(reshape2)
data <- data.frame(Family = paste(phyS@tax_table[,3],phyS@tax_table[,4],phyS@tax_table[,5],phyS@tax_table[,6],sep = "_"), phyS.df)
write.csv(data,"IndicSpecies_IDs_Pacuta.csv", row.names = TRUE)

########################################
### Shared Core Taxa - Venn Diagram  ###
########################################

#load packages
library(eulerr)
library(microbiome)

## Shared taxa between PMAs
#simple way to count the number of samples in each group
table(meta(AL)$PMA, useNA = "always")
#Female   Male   <NA> 
#23     23      0 

#convert to relative abundance
AL.rel <- microbiome::transform(AL, "compositional")

#make a list of PMA types
PMA.all <- unique(as.character(meta(AL.rel)$PMA))

#Write a for loop to go through each of the PMAs one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in PMA.all){ # for each variable n in PMA.all
  #print(paste0("Identifying Core Taxa for ", n))
  
  AL.sub <- subset_samples(AL.rel, PMA.all == n) # Choose sample from PMA.all by n
  
  core_m <- core_members(AL.sub, # AL.sub is phyloseq object selected with only samples from n PMA type
                         detection = 1/100, # 1% in any sample 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each substrate.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)

mycols <- c("#6A6AAD","#ADAD6A") 
plot(venn(list_core), fills = mycols)



