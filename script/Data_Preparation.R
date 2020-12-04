
####Data Preparation
setwd("F:\\#R\\AMF\\Miseq\\Medicago\\R_code_for_github")

metadata <- read.csv(file = "16S_metadata.csv", header = T, row.names = 1)
otu_rdp_cope <- read.csv(file = "feature_16s_copes.csv", header=T) ####该文件使用前注意去重
feature_table <- read.csv(file = "feature_table.csv", header = T)
colSums(feature_table[,5:78]) ###使用的是未标准化的OTU列表

otu_rdp_tax <- read.csv(file = "feature_RDP_taxonomy.csv", header = T)

###(1) Relative 16S OTU abundance (%)
#feature_table

###(2) Quantified 16S OTU abundance (g-1)
names(feature_table)
names(metadata)
match(names(feature_table)[4:78], row.names(metadata))

feature_table_quant <- feature_table
for (i in 1:75) {
  feature_table_quant[,3+i] <- feature_table[,3+i]*metadata$X16S.spike.copes..g.1.[i]/metadata$X16S.spike.reads[i]
}

feature_table_quant$SUM <- rowSums(feature_table_quant[,4:78])
#write.csv(feature_table_quant, file = "feature_table_quant.csv")

###(3) Relative bacteria OTU abundance (%)
match(feature_table$Feature_ID, otu_rdp_cope$Feature.ID)

feature_table_cope_adjusted <- feature_table
for (i in 1:4423) {
  feature_table_cope_adjusted[i,3:78] <- feature_table[i,3:78]/otu_rdp_cope$mean[i]
}
#write.csv(feature_table_cope_adjusted, file = "feature_table_cope_adjusted.csv")


###(4) Quantified bacteria OTU abundance (g-1 )
match(feature_table_quant$Feature_ID, otu_rdp_cope$Feature_ID)

feature_table_quant_cope_adjusted <- feature_table_quant
for (i in 1:4423) {
  feature_table_quant_cope_adjusted[i,3:78] <- feature_table_quant[i,3:78]/otu_rdp_cope$mean[i]
}
#write.csv(feature_table_quant_cope_adjusted, file = "feature_table_quant_cope_adjusted.csv")


###################################################################################################################
###################################################################################################################
###reload the feature tables
feature_table <- read.csv(file = "feature_table.csv", header = T)
feature_table_cope_adjusted <- read.csv(file = "feature_table_cope_adjusted.csv", header = T, row.names = 1)
feature_table_quant <- read.csv(file = "feature_table_quant.csv", header = T, row.names = 1)
feature_table_quant_cope_adjusted <- read.csv(file = "feature_table_quant_cope_adjusted.csv", header = T, row.names = 1)


row.names(feature_table) <- feature_table$OTU_ID
row.names(feature_table_cope_adjusted) <- feature_table_cope_adjusted$OTU_ID
row.names(feature_table_quant) <- feature_table_quant$OTU_ID
row.names(feature_table_quant_cope_adjusted) <- feature_table_quant_cope_adjusted$OTU_ID


###phylum level
###high_abundance phylum	p__Acidobacteria	p__Actinobacteria	p__Bacteroidetes	
###p__Chloroflexi	p__Firmicutes	p__Gemmatimonadetes	p__Verrucomicrobia	
###c__Alphaproteobacteria	c__Betaproteobacteria	c__Deltaproteobacteria	c__Gammaproteobacteria
###low_abundance named others


otu_rdp_phylum1 <- otu_rdp_tax[otu_rdp_tax$V7 %in% c("Acidobacteria", "Actinobacteria", 
                                                     "Bacteroidetes", "Chloroflexi","Firmicutes",
                                                     "Gemmatimonadetes","Verrucomicrobia"),]
otu_rdp_phylum2 <- otu_rdp_tax[otu_rdp_tax$V10 %in% c("Alphaproteobacteria","Betaproteobacteria",
                                                      "Deltaproteobacteria","Gammaproteobacteria"),]
otu_rdp_phylum2$V7 <- otu_rdp_phylum2$V10

otu_rdp_phylum <- rbind(otu_rdp_phylum1,otu_rdp_phylum2)

###(1) Relative 16S phylum abundance (%)
names(otu_rdp_phylum)
names(feature_table)

phylum_table <- merge.data.frame(otu_rdp_phylum[,c(1,9)], feature_table, by="Feature_ID")

phylum <- aggregate(x=phylum_table[,4:79],by=list(phylum_table$V7),FUN=sum)

low_abundance <- colSums(feature_table[,3:78])-colSums(phylum_table[,4:79])
phylum[12,2:77] <- low_abundance

#write.csv(phylum, file = "phylum.csv")


### (2) Relative bacteria phylum abundance (%)
phylum_table_cope_adjusted <- merge.data.frame(otu_rdp_phylum[,c(1,9)], feature_table_cope_adjusted, by="Feature_ID")

phylum_cope_adjusted <- aggregate(x=phylum_table_cope_adjusted[,4:79],by=list(phylum_table_cope_adjusted$V7),FUN=sum)

low_abundance_cope_adjusted <- colSums(feature_table_cope_adjusted[,3:78])-colSums(phylum_table_cope_adjusted[,4:79])

phylum_cope_adjusted[12,2:77] <- low_abundance_cope_adjusted

#write.csv(phylum_cope_adjusted, file = "phylum_cope_adjusted.csv")

###(3) Quantified 16S phylum abundance (g-1)
phylum_table_quant <- merge.data.frame(otu_rdp_phylum[,c(1,9)], feature_table_quant, by="Feature_ID")

phylum_quant <- aggregate(x=phylum_table_quant[,4:79],by=list(phylum_table_quant$V7),FUN=sum)

low_abundance_quant <- colSums(feature_table_quant[,3:78])-colSums(phylum_table_quant[,4:79])

phylum_quant[12,2:77] <- low_abundance_quant

#write.csv(phylum_quant, file = "phylum_quant.csv")


###(4) Quantified bacteria phylum abundance (g-1 )
phylum_table_quant_cope_adjusted <- merge.data.frame(otu_rdp_phylum[,c(1,9)], feature_table_quant_cope_adjusted, by="Feature_ID")

phylum_quant_cope_adjusted <- aggregate(x=phylum_table_quant_cope_adjusted[,4:79],by=list(phylum_table_quant_cope_adjusted$V7),FUN=sum)

low_abundance_quant_cope_adjusted <- colSums(feature_table_quant_cope_adjusted[,3:78])-colSums(phylum_table_quant_cope_adjusted[,4:79])

phylum_quant_cope_adjusted[12,2:77] <- low_abundance_quant_cope_adjusted

#write.csv(phylum_quant_cope_adjusted, file = "phylum_quant_cope_adjusted.csv")
