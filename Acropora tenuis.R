library(plotrix)
library(microbiome)
library(gridExtra)
library(nlme)
library(emmeans)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(multcomp)
library(multcompView)

##############################################################################
#####                       ALPHA DIVERSITY                                ###
##############################################################################

sort(sample_sums(AT)) 

# rarefy
AT_rare<- rarefy_even_depth(AT, sample.size = 9895, rngseed = 1)

# extract metadata from subsetted file
divMeta <- meta(AT_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(AT_rare))
#write.csv(divMeta, "AT_alphadiv.csv")

#Observed ASVs_ANOVA
Obs <- lme(Observed ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) 
#              numDF denDF  F-value p-value
#(Intercept)      1    29 8.281097  0.0074
#PMA              1    29 0.255426  0.6171
#Genotype         3    29 8.242300  0.0004
#PMA:Genotype     3    29 0.117597  0.9491

summarySE(divMeta, groupvars = "Genotype", measurevar = "Observed")

# Obtain least-squares means and confidence intervals
lsm <- emmeans(Obs, "Genotype", adjust = "sidak")

# Get Tukey's letters
lsm_tukey <- cld(lsm, alpha = 0.05)

pairs(emmeans(Obs, "Genotype", adjust="sidak"))
#contrast estimate    SE df t.ratio p.value
#A - B      362.05  95.1 29   3.807  0.0036
#A - C      363.27 100.5 29   3.613  0.0059
#A - D      439.35  95.1 29   4.620  0.0004
#B - C        1.23  97.9 29   0.013  1.0000
#B - D       77.30  92.3 29   0.838  0.8360
#C - D       76.08  97.9 29   0.777  0.8640

#Simpsons_ANOVA
Sim <- lme(InvSimpson ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Sim)
anova(Sim)

#              numDF denDF   F-value p-value
#Genotype         3    29 15.427061  <.0001


lsm <- emmeans(Sim, "Genotype", adjust = "sidak")
lsm_tukey <- cld(lsm, alpha = 0.05)
pairs(emmeans(Sim, "Genotype", adjust="sidak"))
#contrast estimate   SE df t.ratio p.value
#A - B       86.64 15.0 29   5.788  <.0001
#A - C       83.24 15.8 29   5.260  0.0001
#A - D       84.39 15.0 29   5.638  <.0001
#B - C       -3.41 15.4 29  -0.221  0.9961
#B - D       -2.25 14.5 29  -0.155  0.9986
#C - D        1.15 15.4 29   0.075  0.9998

#Shannon_ANOVA
Shan <- lme(Shannon ~ PMA * Genotype, random = ~1|Tank,
            data = divMeta, na.action = na.omit)
summary(Shan)
anova(Shan)

#              numDF denDF  F-value p-value
#Genotype         3    29 14.26984  <.0001

lsm <- emmeans(Shan, "Genotype", adjust = "sidak")
lsm_tukey <- cld(lsm, alpha = 0.05)
pairs(emmeans(Shan, "Genotype", adjust="sidak"))
#contrast estimate   SE df t.ratio p.value
#A - B       2.529 0.437 29   5.791  <.0001
#A - C       2.220 0.462 29   4.809  0.0002
#A - D       2.347 0.437 29   5.375  0.0001
#B - C      -0.309 0.449 29  -0.687  0.9012
#B - D      -0.181 0.424 29  -0.428  0.9731
#C - D       0.127 0.449 29   0.283  0.9919

#reorder levels of PMA
divMeta$PMA <- factor(divMeta$PMA, levels = c("Untreated","PMA-treated"))

##BoxPlot###

AT.obs <- ggplot(divMeta, aes(x=Genotype, y=Observed, fill=Genotype), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#0B5351","#9A879D","#7D84B2","#D2F9F9","#48304D")) + lims(y=c(0,900))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("A. Observed ASVs")
print(AT.obs)

AT.sim <- ggplot(divMeta, aes(x=Genotype, y=InvSimpson, fill=Genotype)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#0B5351","#9A879D","#7D84B2","#D2F9F9","#48304D"))+ lims(y=c(0,150))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("B. Inverse Simpson")
print(AT.sim)

AT.shan <- ggplot(divMeta, aes(x=Genotype, y=Shannon, fill=Genotype)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#0B5351","#9A879D","#7D84B2","#D2F9F9","#48304D"))+ lims(y=c(0,7))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("C. Shannon")
print(AT.shan)


grid.arrange(AT.obs,AT.sim,AT.shan, nrow=1, ncol=3)


rm(AT_rare)
rm(AT.obs)
rm(AT.shan)
rm(AT.sim)
rm(Obs)
rm(Shan)
rm(Sim)
rm(divMeta)
rm(lsm_tukey)
rm(lsm)

##############################################################################
#####                        BETA DIVERSITY                                ###
##############################################################################

##############################
##          PERMANOVA       ##
##############################

library(vegan)
library(pairwiseAdonis)
library(microbiome)


comp <- microbiome::transform(AT, transform = "compositional")

#Generate unweighted Unifrac distance matrix
wunifrac_dist_matrix <- phyloseq::distance(comp, method = "wunifrac")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#weighted Unifrac
dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$PMA)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.3853
#fail to reject the assumption of homogeneity of dispersion by PMA treatment

dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.4436
#fail to reject the assumption of homogeneity of dispersion by colony/genotype

#Because the weighted unifrac matrix met the assumption of homogeneity for PMA, 
#we continue with that distance matrix

#ADONIS test; strata = sets random effect
vegan::adonis2(wunifrac_dist_matrix ~ phyloseq::sample_data(comp)$PMA*phyloseq::sample_data(comp)$Genotype, strata = phyloseq::sample_data(comp)$Tank, permutations = 999) 

#                                                                      Df SumOfSqs      R2      F Pr(>F)    
#phyloseq::sample_data(comp)$PMA                                       1  0.00715 0.00733  0.5032  0.655    
#phyloseq::sample_data(comp)$Genotype                                  3  0.49091 0.50339 11.5174  0.001 ***
#phyloseq::sample_data(comp)$PMA:phyloseq::sample_data(comp)$Genotype  3  0.02250 0.02307  0.5278  0.824    
#Residual                                                             32  0.45465 0.46621                   
#Total                                                                39  0.97521 1.00000                                                                47  0.176569  1.00000                                                              37  1.03150 1.00000                                                              38  12.6609 1.00000 

#pairwise comparisons by genotype
results<-pairwise.adonis(wunifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype, p.adjust.m = "holm", perm = 9999 )
write.csv(results, "pairwisePERMANOVAresults_AT.csv")

rm(dispr.wunifrac)
rm(meta)
rm(wunifrac_dist_matrix)
rm(remove)

##############################################################################
#####                               PCoA                                   ###
##############################################################################


wunifrac <- ordinate(comp, method = "PCoA", distance = "wunifrac")
W <- plot_ordination(physeq = comp,
                     ordination = wunifrac,
                     type = "samples",
                     color = "Genotype") + 
  theme_bw() + scale_color_manual(values=c("#0B5351","#9A879D","#7D84B2","#D2F9F9","#48304D")) +
  geom_point(size=3) 

W + stat_ellipse(type = "t", linetype = 2) + ggtitle("A. Weighted Unifrac PCoA - A. tenuis")

rm(W)
rm(wunifrac)
rm(comp)
rm(remove)


##############################################################################
#####                             BARPLOT                                  ###
############################################################################## 

#Merge by groups for figure

all <- merge_samples(AT, "Genotype")
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
write.csv(genus, "genus_AT.csv", row.names=FALSE)
genus <- read.csv("genus_AT.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

#reorder levels of Genotype
genus$Sample <- factor(genus$Sample, levels = c("D","C","B","A"))



BP_AT <- ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(Genus, Abundance))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance (%) of Bacterial Genera \n") +  
  scale_fill_manual(values = c("#9A879D","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#E56B70","#343434","#1B998B","#48304D","#E3D8F1")) + ggtitle("Acropora tenuis") + coord_flip()

print(BP_AT)
rm(Genus)
rm(genus)
rm(HowMany)
rm(all)
rm(BP_AT)


#Endozoicomonas ASVs by Genotype 

merge <- merge_samples(AT, "Genotype")
comp <- transform_sample_counts(merge, function(x) 100 * x/sum(x))
Endo <- subset_taxa(comp, Genus=="Endozoicomonas")
Endo <- prune_taxa((taxa_sums(Endo) > 0), Endo)
endo <- psmelt(Endo)
endo$OTU <- as.character(endo$OTU)
write.csv(endo, "Endo_AT.csv", row.names=FALSE)
endo <- read.csv("Endo_AT.csv")

#How many levels in endo
HowMany <- length(levels(as.factor(endo$OTU)))

endo$Sample <- factor(endo$Sample, levels = c("D","C","B","A"))



plot_AT <- ggplot(endo, aes(x = Sample, y = Abundance, fill = reorder(OTU, Abundance))) +  
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
  ggtitle("Acropora tenuis") + coord_flip()

print(plot_AT)

rm(merge)
rm(Endo)
rm(endo)
rm(comp)
rm(HowMany)
rm(plot_AT)

########################################
### Shared Core Taxa - Venn Diagram  ###
########################################

#load packages
library(eulerr)
library(microbiome)

## Shared taxa between PMAs
#simple way to count the number of samples in each group
table(meta(AT)$PMA, useNA = "always")

#convert to relative abundance
AT.rel <- microbiome::transform(AT, "compositional")

#make a list of PMA types
PMA.all <- unique(as.character(meta(AT.rel)$PMA))

#Write a for loop to go through each of the PMAs one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in PMA.all){ # for each variable n in PMA.all
  #print(paste0("Identifying Core Taxa for ", n))
  
  AT.sub <- subset_samples(AT.rel, PMA.all == n) # Choose sample from PMA.all by n
  
  core_m <- core_members(AT.sub, # AT.sub is phyloseq object selected with only samples from n PMA type
                         detection = 1/100, # 1% in any sample 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each substrate.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

#print(list_core)

mycols <- c("#6A6AAD","#ADAD6A") 
plot(venn(list_core), fills = mycols)
                                        
