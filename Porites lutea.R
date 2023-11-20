library(plotrix)
library(microbiome)
library(gridExtra)
library(nlme)
library(emmeans)
library(plyr)
library(dplyr)
library(ggplot2)
library(Rmisc)

##############################################################################
#####                       ALPHA DIVERSITY                                ###
##############################################################################

sort(sample_sums(PL)) 

# rarefy
PL_rare<- rarefy_even_depth(PL, sample.size = 4390, rngseed = 1)

# extract metadata from subsetted file
divMeta <- meta(PL_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(PL_rare))
#write.csv(divMeta, "PL_alphadiv.csv")

#Observed ASVs_ANOVA
Obs <- lme(Observed ~ PMA * Genotype, random = ~1|Tank,
           data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) #no significant differences
#get summary statistics
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

PL.obs <- ggplot(divMeta, aes(x=PMA, y=Observed, fill=PMA), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD")) + lims(y=c(0,500))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("A. Observed ASVs")
print(PL.obs)

PL.sim <- ggplot(divMeta, aes(x=PMA, y=InvSimpson, fill=PMA)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD"))+ lims(y=c(0,200))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("B. inverse Simpson")
print(PL.sim)

PL.shan <- ggplot(divMeta, aes(x=PMA, y=Shannon, fill=PMA)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#ADAD6A","#6A6AAD"))+ lims(y=c(0,7))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("C. Shannon")
print(PL.shan)


grid.arrange(PL.obs,PL.sim,PL.shan, nrow=1, ncol=3)


rm(PL_rare)
rm(PL.obs)
rm(PL.shan)
rm(PL.sim)
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

comp <- microbiome::transform(PL, transform = "compositional")

#Generate unweighted Unifrac distance matrix
wunifrac_dist_matrix <- phyloseq::distance(comp, method = "wunifrac")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#weighted Unifrac
dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$PMA)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.2137
#fail to reject the assumption of homogeneity of dispersion by PMA treatment

dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.035
#reject the assumption of homogeneity of dispersion by colony/genotype

#ADONIS test; strata = sets random effect
vegan::adonis2(wunifrac_dist_matrix ~ phyloseq::sample_data(comp)$PMA*phyloseq::sample_data(comp)$Genotype, strata = phyloseq::sample_data(comp)$Tank, permutations = 999) 

#phyloseq::sample_data(comp)$Genotype                                  4 0.034483 0.20116 2.8044  0.001 ***


#pairwise comparisons by genotype
#pairwise comparisons by genotype
results<-pairwise.adonis(wunifrac_dist_matrix, phyloseq::sample_data(comp)$Genotype, p.adjust.m = "holm", perm = 9999 )
write.csv(results, "pairwisePERMANOVAresults_PL.csv")

#   pairs Df   SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#1  A vs B  1 0.007655742 1.923347 0.09653736   0.069       0.69    
#2  A vs C  1 0.013625158 4.454697 0.19838599   0.005       0.05   .
#3  A vs D  1 0.005287729 1.258136 0.06533010   0.246       1.00    
#4  A vs E  1 0.013287256 4.058093 0.18397297   0.006       0.06    
#5  B vs C  1 0.017169017 6.988544 0.27966991   0.001       0.01   *
#6  B vs D  1 0.009526216 2.645465 0.12813783   0.013       0.13    
#7  B vs E  1 0.016499560 6.174075 0.25540066   0.001       0.01   *
#8  C vs D  1 0.007349035 2.743058 0.13223980   0.008       0.08    
#9  C vs E  1 0.006313240 3.606386 0.16691295   0.001       0.01   *
#10 D vs E  1 0.004955820 1.711975 0.08684951   0.078       0.78    

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

W +  ggtitle("A. Porites lutea Weighted Unifrac PCoA") + stat_ellipse(type = "t", linetype = 2)

rm(W)
rm(wunifrac)
rm(comp)


##############################################################################
#####                             BARPLOT                                  ###
############################################################################## 

#Merge by groups for figure

all <- merge_samples(PL, "Genotype")
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

#rename Genera with <5% abundance
genus$Genus[genus$Abundance < 5]<- "< 5% Abundance"
write.csv(genus, "PL_genus.csv", row.names=FALSE)
genus <- read.csv("PL_genus.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

#reorder levels of Genotype
genus$Sample <- factor(genus$Sample, levels = c("E","D","C","B","A"))



BP_PL <- ggplot(genus, aes(x = Sample, y = Abundance, fill = reorder(Genus, Abundance))) +  
  geom_bar(stat = "identity") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  theme(legend.position="bottom", axis.title.y = element_blank()) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = c("#9A879D","#A6CEE3","#f7fff7","#0B5351","#92140C","#7D84B2","#4B7F52","#E56B70","#343434","#1B998B","#48304D","#6A3D9A","#E3D8F1","#1F78B4","#FF7F00","#B2DF8A","#1F3172", "#E7298A","#FFFF99","#FB9A99",
                               "#E31A1C","#33A02C",  "#CAB2D6","#FDBF6F")) + ggtitle("Porites lutea") + coord_flip()

print(BP_PL)

grid.arrange(BP_AM, BP_AT, BP_PD, BP_PL, nrow=4, ncol=1)


rm(Genus)
rm(genus)
rm(HowMany)
rm(all)
rm(BP_PL)

########################################
### Shared Core Taxa - Venn Diagram  ###
########################################

#load packages
library(eulerr)
library(microbiome)

## Shared taxa between PMAs
#simple way to count the number of samples in each group
table(meta(PL)$PMA, useNA = "always")

#convert to relative abundance
PL.rel <- microbiome::transform(PL, "compositional")

#make a list of PMA types
PMA.all <- unique(as.character(meta(PL.rel)$PMA))

#Write a for loop to go through each of the PMAs one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in PMA.all){ # for each variable n in PMA.all
  #print(paste0("Identifying Core Taxa for ", n))
  
  PL.sub <- subset_samples(PL.rel, PMA.all == n) # Choose sample from PMA.all by n
  
  core_m <- core_members(PL.sub, # PL.sub is phyloseq object selected with only samples from n PMA type
                         detection = 1/100, # 1% in any sample 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each substrate.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

#print(list_core)

mycols <- c("#6A6AAD","#ADAD6A") 
plot(venn(list_core), fills = mycols)

## Shared taxa between PMAs
#simple way to count the number of samples in each group
table(meta(coral)$PMA, useNA = "always")

#convert to relative abundance
coral.rel <- microbiome::transform(coral, "compositional")

#make a list of PMA types
PMA.all <- unique(as.character(meta(coral.rel)$PMA))

#Write a for loop to go through each of the PMAs one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in PMA.all){ # for each variable n in PMA.all
  #print(paste0("Identifying Core Taxa for ", n))
  
  coral.sub <- subset_samples(coral.rel, PMA.all == n) # Choose sample from PMA.all by n
  
  core_m <- core_members(coral.sub, # PL.sub is phyloseq object selected with only samples from n PMA type
                         detection = 0, # 1% in any sample 
                         prevalence = 0) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each substrate.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

plot(venn(list_core), fills = mycols)

