library(decontam)
library(microbiome)

# identify contaminants first from PCR negatives
consList <- isContaminant(seqtab = phy, neg = "PCR_Neg", method = "prevalence")

# pull out the names of contaminants
cons <- rownames(consList)[consList$contaminant=="TRUE"]
cons <- as.character(cons) #9 contaminants identified from PCR_Negs

#after viewing potential contaminant list, remove ASVs that are likely NOT contaminants 
cons <- cons[!cons %in% c("e66ea90d7c437d7eeea9da6db1dd8924")]

# - - - - - - - - - - - - - - - - - - - - - - - - - #

# to get info on the contaminants, uncomment the following code to
# run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPer.csv and taxonomy.csv file data

# subset the non PCRNeg coral samples
#vvv <- subset_samples(phy, PCR_Neg == "FALSE")
# merge the samples
#yyy <- merge_samples(vvv, "PCR_Neg", fun = sum)
# transform counts to percentages
#yyy <- transform_sample_counts(yyy, function(x) 100 * x/sum(x))
# extract the cons percentage data
#zzz <- prune_taxa(x = yyy, taxa = cons)
# write otu table to dataframe
#xxx <- data.frame(t(zzz@otu_table))
# write xxx to csv
#write.csv(x = xxx, row.names = TRUE, file = "consPCRNeg.csv")
# subset the contaminant ASVs
#phyCons <- prune_taxa(phy, taxa = cons)
# write the contaminants to a file for reference
#contaminants <-phyCons@tax_table
#contaxa <- contaminants@.Data
#write.csv(contaxa, "contaxaPCR_Neg.csv")

# 8 contaminant ASVs in PCR negatives
# total contamination in the samples = 0.02%


# - - - - - - - - - - - - - - - - - - - - - - - - - #

# remove the contaminants from the main phy phyloseq file
phy <- remove_taxa(phy, taxa = cons)


#Remove PCR_Negs
phy <- subset_samples(phy, PCR_Neg == "FALSE") 


# identify contaminants from extract blanks
consList2 <- isContaminant(seqtab = phy, neg = "Ext_Blank", method = "prevalence")

# pull out the names of contaminants
cons2 <- rownames(consList2)[consList2$contaminant=="TRUE"]
cons2 <- as.character(cons2) #9 contaminants identified from extraction blanks

#after viewing potential contaminant list, remove ASVs that are likely NOT contaminants 
cons2 <- cons2[!cons2 %in% c("ab8f716aab78790040e8f8be1d97f42d")]

# - - - - - - - - - - - - - - - - - - - - - - - - - #

# to get info on the contaminants, uncomment the following code to
# run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPer.csv and taxonomy.csv file data

# subset the non extraction blank coral samples
#vvv2 <- subset_samples(phy, Ext_Blank == "FALSE")
# merge the samples
#yyy2 <- merge_samples(vvv2, group = "Ext_Blank", fun = sum)
# transform counts to percentages
#yyy2 <- transform_sample_counts(yyy2, function(x) 100 * x/sum(x))
# extract the cons percentage data
#zzz2 <- prune_taxa(x = yyy2, taxa = cons2)
# write otu table to dataframe
#xxx2 <- data.frame(t(zzz2@otu_table))
# write xxx to csv
#write.csv(x = xxx2, row.names = TRUE, file = "consExtraction.csv")
# subset the contaminant ASVs
#phyCons2 <- prune_taxa(phy, taxa = cons2)
# write the contaminants to a file for reference
#contaminants2 <-phyCons2@tax_table
#contaxa2 <- contaminants2@.Data
#write.csv(contaxa2, "contaxaExtraction.csv")


# 8 contaminant ASVs in Extraction blanks
# total contamination in the samples from extraction = 0.07%


# remove the contaminants from the main phy phyloseq file
phy <- remove_taxa(phy, taxa = cons2)

# Remove extraction Blanks from phy 
phy <- subset_samples(phy, Ext_Blank=="FALSE")

# identify contaminants from tissue blasting blanks
consList3 <- isContaminant(seqtab = phy, neg = "Blank", method = "prevalence")

# pull out the names of contaminants
cons3 <- rownames(consList3)[consList3$contaminant=="TRUE"]
cons3 <- as.character(cons3) #14 contaminants identified from tissue blasting blanks

#after viewing potential contaminant list, remove ASVs that are likely NOT contaminants 
cons3 <- cons3[!cons3 %in% c("ab8f716aab78790040e8f8be1d97f42d","ef05ab36d891baa224eacda105f67431","3bc9d25cd89ef52b2af6112f44fcbb8d")]

# - - - - - - - - - - - - - - - - - - - - - - - - - #

# to get info on the contaminants, uncomment the following code to
# run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPer.csv and taxonomy.csv file data

# subset the non extraction blank coral samples
#vvv3 <- subset_samples(phy, Blank == "FALSE")
# merge the samples
#yyy3 <- merge_samples(vvv3, group = "Blank", fun = sum)
# transform counts to percentages
#yyy3 <- transform_sample_counts(yyy3, function(x) 100 * x/sum(x))
# extract the cons percentage data
#zzz3 <- prune_taxa(x = yyy3, taxa = cons3)
# write otu table to dataframe
#xxx3 <- data.frame(t(zzz3@otu_table))
# write xxx to csv
#write.csv(x = xxx3, row.names = TRUE, file = "consBlanks.csv")
# subset the contaminant ASVs
#phyCons3 <- prune_taxa(phy, taxa = cons3)
# write the contaminants to a file for reference
#contaminants3 <-phyCons3@tax_table
#contaxa3 <- contaminants3@.Data
#write.csv(contaxa3, "contaxaBlanks.csv")


# 11 contaminant ASVs in processing blanks
# total contamination in the samples from blanks = 0.50%
# total contamination in the samples = 0.59% (n=27)

# remove the contaminants from the main phy phyloseq file
phy <- remove_taxa(phy, taxa = cons3)

# Remove extraction Blanks from phy 
phy <- subset_samples(phy, Blank=="FALSE")

#final dataset
coral <- prune_taxa((taxa_sums(phy) > 0), phy)# 8502 ASVs in 290 samples

#view reads by sample
sort(sample_sums(coral)) 
#6 samples with < 100 reads
#remove these samples as they will not be suitable for future analyses

coral <- prune_samples((sample_sums(coral) > 500), coral)
coral <- prune_taxa((taxa_sums(coral) > 0), coral)# 8501 ASVs in 284 samples
phy_tree(coral) <- root(phy_tree(coral),sample(taxa_names(coral),1), resolve.root = TRUE)


#separate each coral species into their own phyloseq object

AL <- subset_samples(coral, Species=="Acropora loripes")
AL <- prune_taxa((taxa_sums(AL) > 0), AL) #728 ASVs in 46 samples
phy_tree(AL) <- root(phy_tree(AL),sample(taxa_names(AL),1), resolve.root = TRUE)

AM <- subset_samples(coral, Species=="Acropora millepora")
AM <- prune_taxa((taxa_sums(AM) > 0), AM) #1340 ASVs in 49 samples
phy_tree(AM) <- root(phy_tree(AM),sample(taxa_names(AM),1), resolve.root = TRUE)

AT <- subset_samples(coral, Species=="Acropora tenuis")
AT <- prune_taxa((taxa_sums(AT) > 0), AT) #3651 ASVs in 40 samples
phy_tree(AT) <- root(phy_tree(AT),sample(taxa_names(AT),1), resolve.root = TRUE)

PA <- subset_samples(coral, Species=="Pocillopora acuta")
PA <- prune_taxa((taxa_sums(PA) > 0), PA) #623 ASVs in 50 samples
phy_tree(PA) <- root(phy_tree(PA),sample(taxa_names(PA),1), resolve.root = TRUE)

PD <- subset_samples(coral, Species=="Platygyra daedalea")
PD <- prune_taxa((taxa_sums(PD) > 0), PD) #3121 ASVs in 46 samples
phy_tree(PD) <- root(phy_tree(PD),sample(taxa_names(PD),1), resolve.root = TRUE)

PL <- subset_samples(coral, Species=="Porites lutea")
PL <- prune_taxa((taxa_sums(PL) > 0), PL) #3685 ASVs in 50 samples
phy_tree(PL) <- root(phy_tree(PL),sample(taxa_names(PL),1), resolve.root = TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - #

rm(phyCons)
rm(phyCons2)
rm(phyCons3)
rm(vvv)
rm(vvv2)
rm(vvv3)
rm(xxx3)
rm(xxx)
rm(xxx2)
rm(yyy)
rm(yyy2)
rm(yyy3)
rm(zzz3)
rm(zzz)
rm(zzz2)
rm(cons)
rm(cons2)
rm(cons3)
rm(contaminants3)
rm(contaminants)
rm(contaminants2)
rm(contaxa)
rm(contaxa2)
rm(consList)
rm(consList2)
rm(contaxa3)
rm(consList3)
