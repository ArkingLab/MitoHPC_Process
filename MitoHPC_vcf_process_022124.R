## For generating het / homo count, per sample and per pos file from MitoHPC output\
## Stepwise code for easier troubleshoot
## Can work in MitoHPC output directory for convenience
## Fixed filter criteria, not filter out any homoplasmic alt variants

## Latest update: 022124 for UKB 500k WGS

## Files used: 
# mutect2.mutect2.##.vcf
# mutect2.##.suspicious.ids
# mutect2.haplogroup1.tab
# cvg.tab
# count.tab
# new_mito_score.tsv (in case there are variants without MLC score)
# dropouts (latest 020624 version, UKB only)

## Filters for het variants:
# filter column: base|slippage|weak|position|strand
# No INDEL, homopolymer
# read depth < 300
# no mito_lc_score
# no bi-allelic het (same sample and POS but different ALT)

## Files output:
# het/homo variants file
# per sample het/homo counts, CN, MSS and filter
# per pos counts

## Only sections with MODIFY should be modified with 03, 05 or 10


################################# Run


## Install package if necessary
#install.packages("tidyverse")
#install.packages("data.table")

## Load
library(tidyverse)
library(data.table)
library(magrittr)
#library(r)

## Load in files from MitoHPC output (use mutect2.mutect2.##.vcf, ## being AF cutoff for heteroplasmy) (MODIFY)

# 03
#counts = read.table('mutect2.mutect2.03.vcf')
# 05
counts = read.table('mutect2.mutect2.05.vcf')
#10
#counts = read.table('mutect2.mutect2.10.vcf')


## Change col names for MitoHPC outputs
colnames(counts) = c('CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE')
counts$SAMPLE = as.numeric(substr(counts$SAMPLE,1,7))

## Get info from INFO column (get AF, genotype, read depth information from vcf file)
test = unlist(strsplit(counts$INFO, split=';'))
AF = test[grepl('AF',test)]
GT = test[grepl('GT',test)]
DP = test[grepl('DP',test)]
HG = test[grepl('HG',test)]
counts$AF = as.numeric(substr(AF,4,15))
counts$Genotype = substr(GT,4,15)
counts$Read_depth = as.numeric(substr(DP,4,15))
counts$Haplogroup = NA

## Filter out het vs homo, and apply filters to het only, to not filter out homo alt variants
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]
counts = countshet

## Filter for variant level FILTER (current 4 filters, base|slippage|weak|position)
counts = counts %>% filter(!str_detect(FILTER,"strict_strand|strand_bias|base_qual|map_qual|weak_evidence|slippage|position|Homopolymer"))

## Filter for Read_counts < 300
counts = counts[counts$Read_depth >= 300,]

## merge
counts = rbind(counts,countshomo)
rm(countshet,countshomo)

## Filter for variant level INFO (remove INDEL, HP region, flag for NUMT)
counts = counts[!grepl('INDEL',counts$INFO),]
counts = counts[!grepl('Homopolymer',counts$INFO),]

## NUMT
counts1 = counts[!grepl('NUMT',counts$INFO),]
counts2 = counts[grepl('NUMT',counts$INFO),]
counts1$filter_numt = 'notnumt'
counts2$filter_numt = 'numt'
counts = rbind(counts1,counts2)
rm(counts1)
rm(counts2)

## Get genes and complex information (extract gene and complex information for each variant)
countstest = counts[!grepl('DLOOP',counts$INFO),]
countstest = countstest[!grepl('TRN',countstest$INFO),]
countstest = countstest[!grepl('RNR',countstest$INFO),]
countstest = countstest[!grepl('CDS',countstest$INFO),]
countstest$GENE = NA
countstest$COMPLEX = NA
counts$GENE = NA
counts$COMPLEX = NA

countstrn = counts[grepl('TRN',counts$INFO),]
trnlist = unlist(strsplit(countstrn$INFO, split=';'))
trngene = trnlist[grepl('TRN',trnlist)]
countstrn$GENE = substr(trngene,5,15)
countstrn$COMPLEX = 'TRNA'

countscds = counts[grepl('CDS',counts$INFO),]
cdslist = unlist(strsplit(countscds$INFO, split=';'))
cdsgene = cdslist[grepl('CDS',cdslist)]
countscds$GENE = substr(cdsgene,5,15)
cdscomplex = cdslist[grepl('COMPLEX',cdslist)]
countscds$COMPLEX = substr(cdscomplex,9,15)

countsrnr = counts[grepl('RNR=',counts$INFO),]
rnrlist = unlist(strsplit(countsrnr$INFO, split=';'))
rnrgene = rnrlist[grepl('RNR=',rnrlist)]
countsrnr$GENE = substr(rnrgene,5,15)
countsrnr$COMPLEX = 'RRNA'

countsdloop = counts[grepl('DLOOP',counts$INFO),]
countsdloop$GENE = 'DLOOP'
countsdloop$COMPLEX = 'DLOOP'

countsnew = rbind(countstrn, countscds, countsrnr, countsdloop, countstest)
counts = countsnew
rm(countsnew)
rm(countscds)
rm(countsdloop)
rm(countsrnr)
rm(countstest)
rm(countstrn)

## unique (the unique combination of POS_REF_ALT for each variant)
counts$unique = paste(paste(counts$POS, counts$REF, sep = '_'), counts$ALT, sep = '_')

## new mito score (merge unique with mito_lc_score)
scores1 = fread('new_mito_score.tsv')
scores1only = scores1
scores1only$unique = paste(paste(scores1only$POS, scores1only$REF, sep = '_'), scores1only$ALT, sep = '_')
scores1only = scores1only[!duplicated(scores1only$unique),]
scores1only = scores1only[,c(1:3,5:6)]
countsscore = merge(counts, scores1only, by = c('POS','REF','ALT'), all.x = TRUE)
colnames(countsscore) = c(colnames(countsscore)[1:18],'mito_lc_consequence','mito_lc_score')
counts = countsscore
rm(countsscore)
rm(scores1)
rm(scores1only)

## remove mito_lc_score with NA
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]
counts = countshet
counts = counts[!is.na(counts$mito_lc_score),]
counts = rbind(counts,countshomo)

## Other scores (other score from MitoHPC output vcf file) (optional)
test = unlist(strsplit(counts$INFO, split=';'))
AP = test[grepl('AP',test)]
APS = test[grepl('APS',test)]
MLC = test[grepl('MLC',test)]
MCC = test[grepl('MCC',test)]

# ## Apogee score
# AP_label = AP[seq(1,length(AP),by=2)]
# AP_score = AP[seq(2,length(AP),by=2)]
# counts$AP_label = NA
# counts$AP_score = NA
# countsap1 = counts[grepl('APS',counts$INFO),]
# countsap2 = counts[!grepl('APS',counts$INFO),]
# countsap1$AP_label = substr(AP_label,4,500)
# countsap1$AP_score = as.numeric(substr(AP_score,5,20))
# counts = rbind(countsap1,countsap2)
# rm(countsap1)
# rm(countsap2)

# ## Missense score
# test = unlist(strsplit(counts$INFO, split=';'))
# MMC = test[grepl('MMC',test)]
# MMC_score = MMC[seq(2,length(MMC),by=3)]
# MMC_class = MMC[seq(1,length(MMC),by=3)]
# MMC_consq = MMC[seq(3,length(MMC),by=3)]
# counts$MMC_score = NA
# counts$MMC_class = NA
# counts$MMC_consq = NA
# countsmmc1 = counts[grepl('MMC',counts$INFO),]
# countsmmc2 = counts[!grepl('MMC',counts$INFO),]
# countsmmc1$MMC_score = as.numeric(substr(MMC_score,11,50))
# countsmmc1$MMC_class = substr(MMC_class,11,50)
# countsmmc1$MMC_consq = substr(MMC_consq,11,50)
# counts = rbind(countsmmc1,countsmmc2)
# rm(countsmmc1)
# rm(countsmmc2)

# ## protein_gene_missense_constraint; missense_OEUF (MCC)
# test = unlist(strsplit(counts$INFO, split=';'))
# MCC = test[grepl('MCC',test)]
# counts$MCC_score = NA
# countsmcc1 = counts[grepl('MCC',counts$INFO),]
# countsmcc2 = counts[!grepl('MCC',counts$INFO),]
# countsmcc1$MCC_score = as.numeric(substr(MCC,5,20))
# counts = rbind(countsmcc1,countsmcc2)
# rm(countsmcc1)
# rm(countsmcc2)

## Sep homo and het variants
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]

## Filter out those with two ALT allele that add up to 1 (biallelic het)
countshetbi = countshet %>% group_by(SAMPLE,POS) %>% slice_max(mito_lc_score, with_ties = FALSE) %>% ungroup()
countshet = countshetbi
rm(countshetbi)

## Calculate Median and Max (for het variants)
countshet = countshet %>% group_by(unique) %>%
  mutate(Median_AF_Het = median(AF), Max_AF_Het = max(AF)) %>%
  as.data.frame()
countshomo$Median_AF_Het = NA
countshomo$Max_AF_Het = NA
counts = rbind(countshet, countshomo)

## Transition and Transversion (mutation type for all variants)
counts$mutation = paste(counts$REF, counts$ALT, sep = 'to')
counts$Substitution = ifelse((counts$mutation == 'AtoG' | counts$mutation == 'GtoA' | counts$mutation == 'CtoT' | counts$mutation == 'TtoC'),'Transition','Transversion')
counts1 = counts[counts$Substitution == 'Transition',]
counts2 = counts[counts$Substitution == 'Transversion',]
counts2$Substitution = ifelse((counts2$mutation %in% c('AtoC','CtoA','AtoT','TtoA','GtoC','CtoG','GtoT','TtoG')),'Transversion',NA)
counts = rbind(counts1, counts2)

## Non-synonymous variants
counts$mutation_nonsynonymous = NA
counts$mutation_stop = NA
countsnons1 = counts[grepl('NONSYN',counts$INFO),]
countsnons2 = counts[!grepl('NONSYN',counts$INFO),]
countsnons1$mutation_nonsynonymous = 'NONSYN'
countsnons = rbind(countsnons1,countsnons2)
counts = countsnons

## STOP
countsstop1 = counts[grepl('STOP',counts$INFO),]
countsstop2 = counts[!grepl('STOP',counts$INFO),]
countsstop1$mutation_stop = 'STOP'
countsstop = rbind(countsstop1,countsstop2)
counts = countsstop

## Sep into het and homo, also id for downloading if necessary (for UKB)
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]

## Save variant files (MODIFY)
write.table(countshomo, file='cleaned_homo.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(countshet, file='cleaned_het.txt', sep="\t", row.names=FALSE, quote = FALSE)



########## Per sample counts (for each sample, get the het/homo counts, add in filter, CN, MSS, haplogroup)
persamplecount = counts %>% group_by(SAMPLE) %>%
  summarize(count_het = sum(AF != 1),
            count_homo = sum(AF == 1),
            MSS = sum(mito_lc_score[AF != 1])) %>%
  as.data.frame()

## Add in missed samples (those with no het/homo variants in MitoHPC output)
everyone = fread('cvg.tab')
everyone$SAMPLE = as.numeric(substr(everyone$Run,1,7))
allsamplesmissed = everyone$SAMPLE[!everyone$SAMPLE %in% persamplecount$SAMPLE]
missedsampletable = data.frame(SAMPLE = allsamplesmissed, count_het = 0, count_homo = 0, MSS = 0)
persamplecount = rbind(persamplecount, missedsampletable)

## Add flag for suspicious samples (from MitoHPC output suspicious file)

# Run the following in command line used_files directory to split sus.tab file for easier read in
# cat mutect2.05.suspicious.tab | grep mismatch_HG > mutect2.05.suspicious.tab1
# cat mutect2.05.suspicious.tab | grep multiple_NUMTs > mutect2.05.suspicious.tab2
# cat mutect2.05.suspicious.tab | grep -v mismatch_HG | grep -v multiple_NUMTs > mutect2.05.suspicious.tab3
sus1 = fread("mutect2.05.suspicious.tab1",header=FALSE)
sus2 = fread("mutect2.05.suspicious.tab2",header=FALSE)
sus3 = fread("mutect2.05.suspicious.tab3",header=FALSE)
# sus 1 is mismatchHG, sus 2 is multiple NUMTs, sus 3 is haplocheck fail and others
sus = c(sus1$V1,sus2$V1,sus3$V1)
sus = as.numeric(substr(sus,1,7))
persamplecount$sus = ifelse(persamplecount$SAMPLE %in% sus,'Yes','No')

## add filter for dropouts (optional)
dropouts = fread('dropout.txt')
persamplecount$dropout = ifelse(persamplecount$SAMPLE %in% dropouts$V1, 'Yes', 'No')

## add CN for each sample
# Use cvg.tab file from MitoHPC
cn = fread('count.tab')
cn$SAMPLE = as.numeric(substr(cn$Run,1,7))
persamplecount = merge(persamplecount, cn[,c(6,5)], by = 'SAMPLE')
persamplecount$lowCN = ifelse((persamplecount$CN <= 40),'Yes','No')

## haplogroup (read in haplogroup from MitoHPC output and match with sample)
haplogroup = fread('mutect2.haplogroup1.tab')
haplogroup$SAMPLE = as.numeric(substr(haplogroup$Run,1,7))
haplogroup = haplogroup[,c(3,2)]
persamplecount = merge(persamplecount,haplogroup,by='SAMPLE')

## Get het count and MSS for each complex/gene for each sample (optional)
# Complex
# MSS
persamplecomplexmss = countshet %>%
  group_by(SAMPLE, COMPLEX) %>%
  summarize(MSS = sum(mito_lc_score)) %>%
  pivot_wider(names_from = COMPLEX, values_from = MSS, values_fill = 0) %>%
  as.data.frame()
# remove column NA
persamplecomplexmss = persamplecomplexmss %>% select(-'NA')
colnames(persamplecomplexmss)[2:8] = paste('MSS_complex',colnames(persamplecomplexmss)[2:8],sep='_')

# het_count
persamplecomplexhetcount = countshet %>%
  group_by(SAMPLE, COMPLEX) %>%
  summarize(hetcount = n()) %>%
  pivot_wider(names_from = COMPLEX, values_from = hetcount, values_fill = 0) %>%
  as.data.frame()
persamplecomplexhetcount = persamplecomplexhetcount %>% select(-'NA')
colnames(persamplecomplexhetcount)[2:8] = paste('hetcount_complex',colnames(persamplecomplexhetcount)[2:8],sep='_')

# merge
persamplecomplex = merge(persamplecomplexhetcount, persamplecomplexmss, by='SAMPLE')
persamp1 = merge(persamplecount, persamplecomplex, by='SAMPLE', all.x=TRUE)
persamp1[,c(10:23)][is.na(persamp1[,c(10:23)])] = 0

# GENE
# MSS
persamplegenemss = countshet %>%
  group_by(SAMPLE, GENE) %>%
  summarize(MSS = sum(mito_lc_score)) %>%
  pivot_wider(names_from = GENE, values_from = MSS, values_fill = 0) %>%
  as.data.frame()
persamplegenemss = persamplegenemss %>% select(-'NA')
colnames(persamplegenemss)[2:39] = paste('MSS_gene',colnames(persamplegenemss)[2:39],sep='_')

# het_count
persamplegenehetcount = countshet %>%
  group_by(SAMPLE, GENE) %>%
  summarize(hetcount = n()) %>%
  pivot_wider(names_from = GENE, values_from = hetcount, values_fill = 0) %>%
  as.data.frame()
persamplegenehetcount = persamplegenehetcount %>% select(-'NA')
colnames(persamplegenehetcount)[2:39] = paste('hetcount_gene',colnames(persamplegenehetcount)[2:39],sep='_')

# merge
persamplegene = merge(persamplegenehetcount, persamplegenemss, by='SAMPLE')
persamp2 = merge(persamp1, persamplegene, by='SAMPLE', all.x=TRUE)
persamp2[,c(24:99)][is.na(persamp2[,c(24:99)])] = 0

# rename
persamplecount = persamp2

## Save per_sample table (MODIFY)
write.table(persamplecount, file='per_sample.txt', sep="\t", row.names=FALSE, quote = FALSE)


############## Per position counts (get het/homo counts for each POS_REF_ALT combination, order by POS)
## use only variants from filtered samples (no lowCN, dropout or suspicious ids from MitoHPC)
filteredsample = persamplecount %>%
  dplyr::filter(lowCN=='No',dropout=='No',sus=='No',count_het <= 5)
countspos = counts[counts$SAMPLE %in% filteredsample$SAMPLE,]

## Per position stats
perposcount = countspos %>% group_by(unique) %>%
  select(unique,POS,REF,ALT,AF,filter_numt,GENE,COMPLEX,mito_lc_consequence,mito_lc_score,AP_score,MCC_score,Median_AF_Het,
         Max_AF_Het,mutation,Substitution,mutation_nonsynonymous,mutation_stop) %>%
  mutate(count_het = sum(AF != 1),
         count_homo = sum(AF == 1)) %>%
  slice(1) %>%
  as.data.frame()
perposcount = perposcount[!duplicated(perposcount$unique),c(1:4,19:20,6:18)]
perposcount = perposcount[order(perposcount$POS),]

## Save (MODIFY)
write.table(perposcount, file='per_position.txt', sep="\t", row.names=FALSE, quote = FALSE)