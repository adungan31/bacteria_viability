library(plotrix)
library(microbiome)
library(gridExtra)
library(nlme)
library(emmeans)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(multcomp)

##############################################################################
#####                       ALPHA DIVERSITY                                ###
##############################################################################

sort(sample_sums(PD)) 

# rarefy
PD_rare<- rarefy_even_depth(PD, sample.size = 10489, rngseed = 1)

# extract metadata from subsetted file
divMeta <- meta(PD_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(PD_rare))
#write.csv(divMeta, "PD_alphadiv.csv")

#Observed ASVs_ANOVA
Obs <- lme(Observed ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) 
#            numDF denDF  F-value p-value
#Genotype         4    31 3.168752  0.0271

summarySE(divMeta, groupvars = "Genotype", measurevar = "Observed")

# Obtain least-squares means and confidence intervals
lsm <- emmeans(Obs, "Genotype", adjust = "sidak")

# Get Tukey's letters
lsm_tukey <- cld(lsm, alpha = 0.05)

pairs(emmeans(Obs, "Genotype", adjust="sidak"))
#contrast estimate    SE df t.ratio p.value
# A - C      197.62 63.8 33   3.098  0.0302

#Simpsons_ANOVA
Sim <- lme(InvSimpson ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Sim)
anova(Sim)

#              numDF denDF   F-value p-value
#Genotype         4    31  8.46252  0.0001

lsm <- emmeans(Sim, "Genotype", adjust = "sidak")
lsm_tukey <- cld(lsm, alpha = 0.05)
pairs(emmeans(Sim, "Genotype", adjust="sidak"))
#contrast estimate   SE df t.ratio p.value
#A - B     -1.5080 0.288 33  -5.242  0.0001
#B - C      1.2655 0.297 33   4.268  0.0014
#B - D      1.1383 0.288 33   3.957  0.0033
#B - E      1.4149 0.423 33   3.342  0.0166

#Shannon_ANOVA
Shan <- lme(Shannon ~ PMA * Genotype, random = ~1|Tank,
            data = divMeta, na.action = na.omit)
summary(Shan)
anova(Shan)

#              numDF denDF  F-value p-value
#Genotype         4    31  3.402266  0.0203

lsm <- emmeans(Shan, "Genotype", adjust = "sidak")
lsm_tukey <- cld(lsm, alpha = 0.05)
pairs(emmeans(Shan, "Genotype", adjust="sidak"))
#contrast estimate   SE df t.ratio p.value
#A - B     -1.0173 0.331 33  -3.073  0.0321
#B - C      1.0467 0.341 33   3.068  0.0325

#reorder levels of Genotype
divMeta$PMA <- factor(divMeta$Genotype, levels = c("E","D","C","B","A"))

##BoxPlot###

PD.obs <- ggplot(divMeta, aes(x=Genotype, y=Observed, fill=Genotype), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#0B5351","#9A879D","#7D84B2","#D2F9F9","#48304D")) + lims(y=c(0,500))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("D. Observed ASVs")
print(PD.obs)

PD.sim <- ggplot(divMeta, aes(x=Genotype, y=InvSimpson, fill=Genotype)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#0B5351","#9A879D","#7D84B2","#D2F9F9","#48304D"))+ lims(y=c(0,4.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("E. inverse Simpson")
print(PD.sim)

PD.shan <- ggplot(divMeta, aes(x=Genotype, y=Shannon, fill=Genotype)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#0B5351","#9A879D","#7D84B2","#D2F9F9","#48304D"))+ lims(y=c(0,3.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("F. Shannon")
print(PD.shan)


grid.arrange(AT.obs,AT.sim,AT.shan, PD.obs,PD.sim,PD.shan, nrow=2, ncol=3)


rm(PD_rare)
rm(PD.obs)
rm(PD.shan)
rm(PD.sim)
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

comp <- microbiome::transform(PD, transform = "compositional")

#Generate unweighted Unifrac distance matrix
wunifrac_dist_matrix <- phyloseq::distance(comp, method = "wunifrac")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#weighted Unifrac
dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$PMA)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.9836
#fail to reject the assumption of homogeneity of dispersion by PMA treatment

dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=3.266e-10
#reject the assumption of homogeneity of dispersion by colony/genotype

#ADONIS test; strata = sets random effect
vegan::adonis2(wunifrac_dist_matrix ~ phyloseq::sample_data(comp)$PMA*phyloseq::sample_data(comp)$Genotype, strata = phyloseq::sample_data(comp)$Tank, permutations = 999) 

#                                                                      Df SumOfSqs      R2      F Pr(>F)    
#phyloseq::sample_data(comp)$PMA                                       1  0.00124 0.00174  0.4348  0.616    
#phyloseq::sample_data(comp)$Genotype                                  4  0.60213 0.84855 52.8800  0.001 ***
#phyloseq::sample_data(comp)$PMA:phyloseq::sample_data(comp)$Genotype  4  0.00375 0.00528  0.3293  0.909    
#Residual                                                             36  0.10248 0.14442                   
#Total                                                                45  0.70959 1.00000 

#pairwise comparisons by genotype
results<-pairwise.adonis(wunifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype, p.adjust.m = "holm", perm = 9999 )
write.csv(results, "pairwisePERMANOVAresults_PD.csv")

rm(dispr.wunifrac)
rm(wunifrac_dist_matrix)


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

W + ggtitle("A. Weighted Unifrac PCoA") stat_ellipse(type = "t", linetype = 2)

rm(W)
rm(wunifrac)
rm(comp)


##############################################################################
#####                             BARPLOT                                  ###
############################################################################## 

#Merge by groups for figure

all <- merge_samples(PD, "Genotype")
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

#rename Genera with <1% abundance
genus$Genus[genus$Abundance < 1]<- "< 1% Abundance"
write.csv(genus, "PD_genus.csv", row.names=FALSE)
genus <- read.csv("PD_genus.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

#reorder levels of Genotype
genus$Sample <- factor(genus$Sample, levels = c("E","D","C","B","A"))


BP_PD <- ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(Genus, Abundance))) +  
  geom_bar(stat = "identity") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = c("#9A879D","#A6CEE3","#f7fff7","#0B5351","#92140C","#7D84B2","#E56B70","#E3D8F1","#1F78B4","#FF7F00","#B2DF8A","#1F3172", "#E7298A","#FFFF99","#FB9A99",
                               "#E31A1C","#33A02C",  "#CAB2D6","#FDBF6F")) + ggtitle("Platygyra daedalea") + coord_flip()

print(BP_PD)

#ASV level
Genus <- transform_sample_counts(all, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$OTU <- as.character(genus$OTU)
write.csv(genus, "PD_ASV.csv", row.names=FALSE)
genus <- read.csv("PD_ASV.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$OTU)))

#reorder levels of Genotype
genus$Sample <- factor(genus$Sample, levels = c("E","D","C","B","A"))


BP <- ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(OTU, Abundance))) +  
  geom_bar(stat = "identity") +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance (%) of Bacterial ASVs \n") +  
  scale_fill_manual(values = c("#9A879D","#A6CEE3","#f7fff7","#0B5351","#92140C","#7D84B2","#E56B70","#E3D8F1")) + ggtitle("B. ASV level Barplot")

BP + coord_flip()

rm(Genus)
rm(genus)
rm(HowMany)
rm(all)
rm(BP)

#Endozoicomonas ASVs by Genotype 

merge <- merge_samples(PD, "Genotype")
comp <- transform_sample_counts(merge, function(x) 100 * x/sum(x))
Endo <- subset_taxa(comp, Genus=="Endozoicomonas" | Genus=="Unknown Rhodanobacteraceae")
Endo <- prune_taxa((taxa_sums(Endo) > 0), Endo)
endo <- psmelt(Endo)
endo$OTU <- as.character(endo$OTU)
write.csv(endo, "Endo_PD.csv", row.names=FALSE)
endo <- read.csv("Endo_PD.csv")

#How many levels in endo
HowMany <- length(levels(as.factor(endo$OTU)))

endo$Sample <- factor(endo$Sample, levels = c("E","D","C","B","A"))



plot_PD <- ggplot(endo, aes(x = Sample, y = Abundance, fill = reorder(OTU, Abundance))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance (%) of Endozoicomonas and Unknown Rhodanobacteraceae ASVs \n") + scale_fill_manual(values = c("#9A879D","#92140C","#FB9A99","#E3D8F1","#E56B70","#343434","#1B998B","#48304D",
                                                                                     "#A6CEE3","#6A3D9A","#1F78B4","#FF7F00","#CAB2D6","#33A02C","#1F3172", "#E7298A","#FDBF6F")) +
  ggtitle("Platygyra daedalea") + coord_flip()

print(plot_PD) 

grid.arrange(plot_AL, plot_AM, plot_AT, plot_PD, nrow=4, ncol=1)


rm(merge)
rm(Endo)
rm(endo)
rm(comp)
rm(HowMany)
rm(plot_PD)

########################################
### Shared Core Taxa - Venn Diagram  ###
########################################

#load packages
library(eulerr)
library(microbiome)

## Shared taxa between PMAs
#simple way to count the number of samples in each group
table(meta(PD)$PMA, useNA = "always")

#convert to relative abundance
PD.rel <- microbiome::transform(PD, "compositional")

#make a list of PMA types
PMA.all <- unique(as.character(meta(PD.rel)$PMA))

#Write a for loop to go through each of the PMAs one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in PMA.all){ # for each variable n in PMA.all
  #print(paste0("Identifying Core Taxa for ", n))
  
  PD.sub <- subset_samples(PD.rel, PMA.all == n) # Choose sample from PMA.all by n
  
  core_m <- core_members(PD.sub, # PD.sub is phyloseq object selected with only samples from n PMA type
                         detection = 1/100, # 1% in any sample 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each substrate.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

#print(list_core)

mycols <- c("#6A6AAD","#ADAD6A") 
plot(venn(list_core), fills = mycols)
