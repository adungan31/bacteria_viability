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

sort(sample_sums(PA)) 

# rarefy
PA_rare<- rarefy_even_depth(PA, sample.size = 883, rngseed = 1)

# extract metadata from subsetted file
divMeta <- meta(PA_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(PA_rare))
#write.csv(divMeta, "PA_alphadiv.csv")

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


#Shannon_ANOVA
Shan <- lme(Shannon ~ PMA * Genotype, random = ~1|Tank,
            data = divMeta, na.action = na.omit)
summary(Shan)
anova(Shan)


#reorder levels of PMA
divMeta$PMA <- factor(divMeta$PMA, levels = c("Untreated","PMA-treated"))

##BoxPlot###

PA.obs <- ggplot(divMeta, aes(x=PMA, y=Observed, fill=PMA), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD")) + lims(y=c(0,40))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("A. Observed ASVs")
print(PA.obs)

PA.sim <- ggplot(divMeta, aes(x=PMA, y=InvSimpson, fill=PMA)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD"))+ lims(y=c(0,18))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("B. inverse Simpson")
print(PA.sim)

PA.shan <- ggplot(divMeta, aes(x=PMA, y=Shannon, fill=PMA)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD"))+ lims(y=c(0,4))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("C. Shannon")
print(PA.shan)


grid.arrange(PA.obs,PA.sim,PA.shan, nrow=1, ncol=3)


rm(PA_rare)
rm(PA.obs)
rm(PA.shan)
rm(PA.sim)
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

comp <- microbiome::transform(PA, transform = "compositional")

#Generate unweighted Unifrac distance matrix
wunifrac_dist_matrix <- phyloseq::distance(comp, method = "wunifrac")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#weighted Unifrac
dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$PMA)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.5053
#fail to reject the assumption of homogeneity of dispersion by PMA treatment

dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.1963
#fail to reject the assumption of homogeneity of dispersion by colony/genotype


#ADONIS test; strata = sets random effect
vegan::adonis2(wunifrac_dist_matrix ~ phyloseq::sample_data(comp)$PMA*phyloseq::sample_data(comp)$Genotype, strata = phyloseq::sample_data(comp)$Tank, permutations = 9999) 

#                                                                     Df SumOfSqs     R2      F Pr(>F)   
#phyloseq::sample_data(comp)$PMA                                       1  0.07170 0.03657 1.8771 0.0106 *
#phyloseq::sample_data(comp)$Genotype                                  4  0.22073 0.11260 1.4448 0.4598  
#phyloseq::sample_data(comp)$PMA:phyloseq::sample_data(comp)$Genotype  4  0.14008 0.07146 0.9169 0.6832  
#Residual                                                             40  1.52778 0.77937                
#Total                                                                49  1.96028 1.00000                                                           49  2.27409 1.0000                                                                   47  0.176569  1.00000                                                              37  1.03150 1.00000                                                              38  12.6609 1.00000 

rm(dispr.wunifrac)
rm(wunifrac_dist_matrix)

##############################################################################
#####                               PCoA                                   ###
##############################################################################

#remove outlier sample for visualization
remove <- microbiome::remove_samples("Pa1.1A_FWD09", PA)
comp <- microbiome::transform(remove, transform = "compositional")

wunifrac <- ordinate(comp, method = "PCoA", distance = "wunifrac")
W <- plot_ordination(physeq = comp,
                     ordination = wunifrac,
                     type = "samples",
                     color = "PMA",
                     shape = "PMA") + 
  theme_bw() + scale_color_manual(values=c("#6A6AAD","#ADAD6A")) +
  geom_point(size=3) 

W  + ggtitle("A. Pocillopora acuta Weighted Unifrac PCoA") + stat_ellipse(type = "t", linetype = 2)

rm(W)
rm(wunifrac)
rm(comp)
rm(remove)


##############################################################################
#####                             BARPLOT                                  ###
############################################################################## 

#Merge by groups for figure

all <- merge_samples(PA, "PMA")
all <- prune_taxa((taxa_sums(all) > 0), all)


##Genus level##
# To represent at the Genus level (instead of ASV level)
Genus <- tax_glom(all,taxrank = "Genus")


# Transform counts in relative abundance and select most abundant families
Genus <- transform_sample_counts(Genus, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$Genus <- as.character(genus$Genus)

#rename Genera with <3% abundance
genus$Genus[genus$Abundance < 3]<- "< 3% Abundance"
write.csv(genus, "PA_genus.csv", row.names=FALSE)
genus <- read.csv("PA_genus.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

# Define 9 distinguishable colors
my_colors <- c("#9A879D","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#D2F9F9","#FB9A99","#E3D8F1")



BP <- ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(Genus, Abundance))) +  
  geom_bar(stat = "identity")  +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = my_colors) + ggtitle("B. P.acuta Genus level Barplot")

BP + coord_flip()


##ASV level##
# To represent at the Genus level (instead of ASV level)

# Transform counts in relative abundance and select most abundant families
Genus <- transform_sample_counts(all, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$OTU <- as.character(genus$OTU)
write.csv(genus, "PA_ASV.csv", row.names=FALSE)
genus <- read.csv("PA_ASV.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$OTU)))

# Define 14 distinguishable colors
my_colors <- c("#9A879D","#E56B70","#A6CEE3","#343434","#1B998B","#48304D","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#D2F9F9","#FB9A99","#E3D8F1")


BP <- ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(OTU, Abundance))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = my_colors) + ggtitle("B. P.acuta ASV level Barplot")

BP + coord_flip()

rm(Genus)
rm(genus)
rm(HowMany)
rm(all)
rm(my_colors)
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
rel <- transform_sample_counts(PA, function(x) 100 * x/sum(x))

#order metadata by variable of interest and check metadata table to confirm
rel <- set_sample_order(rel, 'PMA')
meta <- microbiome::meta(rel)
table(meta$PMA)#use these numbers and order for the 'group' below 

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

group <- c(rep("PMA-treated", 25),rep("Untreated", 25))


# Performing the indicator value analysis (omit sampNames column in asv table, otherwise indval fails)
indval = multipatt(asv.table[,-1], group, control = how(nperm=9999))

capture.output(summary(indval, At = 0.6, Bt = 0.6, invdalcomp = TRUE, alpha = 0.05, func = 'Indval.g'), file = "IndvalOutput_PMA_PA.txt")


##Genus LEVEL###

Genus <- tax_glom(PA,taxrank = "Genus")

# This package does not work with data in Phyloseq format, so we will have to adjust that.
rel <- transform_sample_counts(Genus, function(x) 100 * x/sum(x))

#order metadata by variable of interest and check metadata table to confirm
rel <- set_sample_order(rel, 'PMA')

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


# Performing the indicator value analysis (omit sampNames column in asv table, otherwise indval fails)
indval = multipatt(asv.table[,-1], group, control = how(nperm=9999))

capture.output(summary(indval, At = 0.6, Bt = 0.6, invdalcomp = TRUE, alpha = 0.05, func = 'Indval.g'), file = "IndvalOutput_PMA_PA_Genus.txt")

rm(rel)
rm(Genus)
rm(asv.table)
rm(indval)
rm(meta)
rm(group)
rm(sampNames)

#Endozoicomonas ASVs by PMA 

merge <- merge_samples(PA, "PMA")
comp <- transform_sample_counts(merge, function(x) 100 * x/sum(x))
Endo <- subset_taxa(comp, Genus=="Endozoicomonas")
Endo <- prune_taxa((taxa_sums(Endo) > 0), Endo)
endo <- psmelt(Endo)
endo$OTU <- as.character(endo$OTU)
write.csv(endo, "Endo_PA.csv", row.names=FALSE)
endo <- read.csv("Endo_PA.csv")

#How many levels in endo
HowMany <- length(levels(as.factor(endo$OTU)))

endo$Sample <- factor(endo$Sample, levels = c("D","C","B","A"))



plot <- ggplot(endo, aes(x = Sample, y = Abundance, fill = reorder(OTU, Abundance))) +  
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
  ggtitle("Endozoicomonas only Barplot - A.tenuis")

plot + coord_flip()

rm(merge)
rm(Endo)
rm(endo)
rm(comp)
rm(HowMany)
rm(plot)

########################################
### Shared Core Taxa - Venn Diagram  ###
########################################

#load packages
library(eulerr)
library(microbiome)

## Shared taxa between PMAs
#simple way to count the number of samples in each group
table(meta(PA)$PMA, useNA = "always")

#convert to relative abundance
PA.rel <- microbiome::transform(PA, "compositional")

#make a list of PMA types
PMA.all <- unique(as.character(meta(PA.rel)$PMA))

#Write a for loop to go through each of the PMAs one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in PMA.all){ # for each variable n in PMA.all
  #print(paste0("Identifying Core Taxa for ", n))
  
  PA.sub <- subset_samples(PA.rel, PMA.all == n) # Choose sample from PMA.all by n
  
  core_m <- core_members(PA.sub, # PA.sub is phyloseq object selected with only samples from n PMA type
                         detection = 1/100, # 1% in any sample 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each substrate.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

#print(list_core)

mycols <- c("#6A6AAD","#ADAD6A") 
plot(venn(list_core), fills = mycols)
