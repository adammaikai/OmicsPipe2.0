## ---- snpAnnotation_DNA ----
snpAnnotation_DNA <- function(){

  
### Import files
# import vcf
#snpVCF <- read.csv(paste(config$VARIANT_RESULTS, '/', patientID, '/', patientID, '.vcf', sep=''), stringsAsFactors=FALSE, header=F, sep='\t', comment.char="#")
snpVCF <- read.csv(paste(config$VARIANT_RESULTS, '/', patientID,  '/intogen_input.vcf', sep=''), stringsAsFactors=FALSE, header=F, sep='\t', comment.char="#")

# import Intogen consequences and variant_genes
intogenSNPs <- read.csv2(paste(config$INTOGEN_RESULTS, '/', patientID, '/consequences.tsv', sep=''), stringsAsFactors=FALSE, header=T, sep="\t", na.strings="", as.is=T)
intogenGenes <- read.csv2(paste(config$INTOGEN_RESULTS, '/', patientID, '/variant_genes.tsv', sep=''), stringsAsFactors=FALSE, header=T, sep="\t", na.strings="", as.is=T)

# import COSMIC complete export
CosmicCompleteExport <- read.csv2(config$COSMIC, stringsAsFactors=FALSE, sep='\t')

# import ClinVar variant summary
clinVarSummary <- read.csv2(config$CLINVAR, stringsAsFactors=FALSE, sep='\t')

# import pharmGkb clinical annotations per rsID and info for each SNP
pharmgkbRSID <- read.csv2(config$PHARMGKB_rsID, stringsAsFactors=FALSE, sep=',')
pharmgkbAllele <- read.csv2(config$PHARMGKB_Allele, stringsAsFactors=FALSE, sep='\t')

# import DrugBank stable data file
drugBankTable <- read.csv2(config$DRUGBANK, stringsAsFactors=FALSE, sep='\t')

#import SNPeff variants
snp_eff <- read.csv2(paste(config$VARIANT_RESULTS, '/', patientID, '/', patientID, '_final_variants_filt_snpeff.vcf', sep=''), stringsAsFactors=FALSE, header=F, sep='\t', comment.char="#")

#import Annovar variants
annovar <- read.csv2(paste(config$VARIANT_RESULTS, '/', patientID, '/', patientID, '.hg19_multianno.vcf', sep=''), stringsAsFactors=FALSE, header=F, sep='\t', comment.char="#")

### Format VCF
# Rename columns
names(snpVCF) <- c("CHROM", "POS", "rsID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")

# Remove chr from chromosome column
snpVCF$CHROM <- gsub("chr", "", snpVCF$CHROM)

# Create Alternative allele columns for matching
snpVCF$ALT1 <- sapply(strsplit(as.character(snpVCF$ALT), ","), "[", 1)
snpVCF$ALT2 <- sapply(strsplit(as.character(snpVCF$ALT), ","), "[", 2)
snpVCF$ALT3 <- sapply(strsplit(as.character(snpVCF$ALT), ","), "[", 3)

x <- colsplit(snpVCF$SAMPLE, ":", c("genotype", "geno2", "geno3", "geno4", "geno5"))
y <- colsplit(x$genotype, "/", c("allele1value", "allele2value"))
snpVCF <- cbind(snpVCF, "genotype" = x$genotype, y)

# Fill Genotype allele columns and Delete allele value columns
snpVCF$allele1 <- ifelse(snpVCF$allele1value == 0, snpVCF$REF,
                           ifelse(snpVCF$allele1value == 1, snpVCF$ALT1, 
                                  ifelse(snpVCF$allele1value == 2, snpVCF$ALT2,
                                         ifelse(snpVCF$allele1value == 3, snpVCF$ALT3, NA))))
snpVCF$allele2 <- ifelse(snpVCF$allele2value == 0, snpVCF$REF,
                           ifelse(snpVCF$allele2value == 1, snpVCF$ALT1,
                                  ifelse(snpVCF$allele2value == 2, snpVCF$ALT2,
                                         ifelse(snpVCF$allele2value == 3, snpVCF$ALT3, NA))))
snpVCF$allele1value <- NULL
snpVCF$allele2value <- NULL


### Label Mutation Type
# Create EFF column in vcf with single EFF for each snp
library(data.table)
snpVCF$EFF <- gsub("EFF=|;.*$", "", regmatches(snpVCF$INFO, regexpr("EFF=.*\\)", snpVCF$INFO)))
#dt <- data.table(snpVCF)
#snpVCFexpanded <- dt[, list(EFF = unlist(strsplit(EFF, ","))), by = POS]
#snpVCF <- merge(snpVCFexpanded, snpVCF, by="POS")
#snpVCF$EFF.y <- NULL
#setnames(snpVCF, "EFF.x", "EFF")
#snpVCF <- as.data.frame(snpVCF)
#snpVCF <- snpVCF[,c("CHROM", "POS", "rsID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE", "ALT1", "ALT2", "ALT3", "allele1", "allele2", "EFF")]

# Divide EFF info into separate columns
snpVCF$Effect_type <- gsub("\\(", "", regmatches(snpVCF$EFF, regexpr("^.*\\(", snpVCF$EFF)))
snpVCF$Effect_info <- gsub("\\(|\\)", "", regmatches(snpVCF$EFF, regexpr("\\(.*\\)", snpVCF$EFF)))
snpVCF$Effect_impact <- sapply(strsplit(as.character(snpVCF$Effect_info), "\\|"), "[", 1)
snpVCF$Functional_class <- sapply(strsplit(as.character(snpVCF$Effect_info), "\\|"), "[", 2)
snpVCF$Codon_change <- sapply(strsplit(as.character(snpVCF$Effect_info), "\\|"), "[", 3)
snpVCF$AA_change <- sapply(strsplit(as.character(snpVCF$Effect_info), "\\|"), "[", 4)
snpVCF$Gene <- sapply(strsplit(as.character(snpVCF$Effect_info), "\\|"), "[", 6)
snpVCF$Coding <- sapply(strsplit(as.character(snpVCF$Effect_info), "\\|"), "[", 8)

# Aggregate Codon_change and AA_change and merge back with snpVCF
snpVCFchangesAgged <- aggregate(snpVCF[,c("Codon_change","AA_change")], snpVCF[,c("CHROM", "POS", "Gene")], function(x) paste(unique(x[!is.na(x)]), " ", sep="", collapse=""))
snpVCFchangesAgged <- merge(snpVCF, snpVCFchangesAgged, by=c("CHROM", "POS", "Gene"))

# Format snpVCFchangesAgged and create final Effect column
snpVCFchangesAgged$Codon_change.x <- NULL
snpVCFchangesAgged$AA_change.x <- NULL
setnames(snpVCFchangesAgged, "Codon_change.y", "Codon_change")
setnames(snpVCFchangesAgged, "AA_change.y", "AA_change")
snpVCFchangesAgged$Effect <- paste(snpVCFchangesAgged$Effect_impact, snpVCFchangesAgged$Functional_class, snpVCFchangesAgged$Effect_type, snpVCFchangesAgged$Codon_change, snpVCFchangesAgged$AA_change, sep = " ")
snpVCFchangesAgged$Effect <- gsub("\\>\\s*$", "", snpVCFchangesAgged$Effect)
snpVCFchangesAgged$Effect <- paste(snpVCFchangesAgged$Effect, ";  ", sep="")
snpVCFchangesAgged$Effect <- gsub("\\s+\\<", ", ", snpVCFchangesAgged$Effect)
snpVCFchangesAgged$Effect <- paste(snpVCFchangesAgged$Gene, snpVCFchangesAgged$Effect, sep=": ")

# Aggregate by position to give final VCF
snpVCFchangesAgged <- snpVCFchangesAgged[order(snpVCFchangesAgged$Effect_impact, snpVCFchangesAgged$Coding),]
snpVCFbyPos <- aggregate(snpVCFchangesAgged[,c("rsID", "Gene", "Effect", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE", "ALT1", "ALT2", "ALT3", "allele1", "allele2")], snpVCFchangesAgged[,c("CHROM", "POS")], function(x) paste(unique(x[!is.na(x)]), " ", sep="", collapse="")) 

# Create separate gene columns
snpVCFbyPos$GENE1 <- sapply(strsplit(as.character(snpVCFbyPos$Gene), "\\s"), "[", 1)
snpVCFbyPos$GENE2 <- sapply(strsplit(as.character(snpVCFbyPos$Gene), "\\s"), "[", 2)
snpVCFbyPos$GENE3 <- sapply(strsplit(as.character(snpVCFbyPos$Gene), "\\s"), "[", 3)
snpVCFbyPos$GENE4 <- sapply(strsplit(as.character(snpVCFbyPos$Gene), "\\s"), "[", 4)
snpVCFbyPos$GENE5 <- sapply(strsplit(as.character(snpVCFbyPos$Gene), "\\s"), "[", 5)

# Remove hanging spaces created by aggregate function
for(i in (1:ncol(snpVCFbyPos))){
  snpVCFbyPos[,i] <- gsub("\\>\\s+$", "", snpVCFbyPos[,i])
}

# Format gene column and order VCF by chromosome and position
snpVCFbyPos$Gene <- gsub("\\>\\s+\\<", "; ", snpVCFbyPos$Gene)
snpVCFbyPos <- snpVCFbyPos[order(snpVCFbyPos$CHROM, snpVCFbyPos$POS),]


### Intogen Annotation
## Annotate with intogen_consequences.tsv (mainly for gene names and mutation type)
# Format names and ALLELE column for matching
names(intogenSNPs)[2] <- "CHROM"
names(intogenSNPs)[4] <- "POS"
intogenSNPs$Ref <- sapply(strsplit(as.character(intogenSNPs$ALLELE), "/"), "[", 1)
intogenSNPs$Alt <- sapply(strsplit(as.character(intogenSNPs$ALLELE), "/"), "[", 2)


# Merge by Chromosome and Position
intogenSNPsAnno <- merge(snpVCFbyPos, intogenSNPs, by=c("CHROM", "POS"))

# Subset to double check alleles and select only GENE_ID from intogenSNPs
intogenSNPsAnno <- subset(intogenSNPsAnno, ((REF==Ref) & ((ALT1==Alt) | (ALT2==Alt) | (ALT3==Alt))), select = c(CHROM:GENE5, GENE_ID))

## Annotate with intogen_variant_genes.tsv based on Ensembl Gene ID (mainly for driver and impact)
# Remove first five columns
intogenGenes <- intogenGenes[,6:16]

# Merge by Ensembl Gene ID
intogenGenesAnno <- merge(intogenSNPsAnno, intogenGenes, by="GENE_ID")

# Create driver "yes or no" column
intogenGenesAnno$Intogen_driver <- ifelse(intogenGenesAnno$INTOGEN_DRIVER == 0, "No",
                           ifelse(intogenGenesAnno$INTOGEN_DRIVER == 1, paste("Yes", "(", intogenGenesAnno$Gene, ")", sep=""), NA))

# Replace gene names that are NA with Ensembl ID
for(i in (1:nrow(intogenGenesAnno))){
  if(is.na(intogenGenesAnno$SYMBOL[i])){
    intogenGenesAnno$SYMBOL[i] <- intogenGenesAnno$GENE_ID[i]
  }
}

# Paste columns to create "intogen_info" column
intogenGenesAnno$Intogen_info <- paste(intogenGenesAnno$SYMBOL, ": ", "Variant Impact: ", intogenGenesAnno$VAR_IMPACT_DESC, ", Impact Score: ", intogenGenesAnno$VAR_IMPACT, ";", sep="")

# Aggregate by chromosome and position, keeping only intogen_info and driver columns
intogenGenesAnno <- aggregate(intogenGenesAnno[,c("rsID", "Gene", "Effect", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE", "ALT1", "ALT2", "ALT3", "allele1", "allele2", "GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "Intogen_driver", "Intogen_info")], intogenGenesAnno[,c("CHROM", "POS")], function(x) paste(unique(x[!is.na(x)]), " ", sep="", collapse=""))

## Format for future
# Order by chromosome and position
finalIntogenSNPs <- intogenGenesAnno[order(intogenGenesAnno$CHROM, intogenGenesAnno$POS),]

# Remove "No" from Intogen_driver column if yes is already there
for(i in (1:nrow(finalIntogenSNPs))){
  if(grepl("Yes", finalIntogenSNPs$Intogen_driver[i])){
    finalIntogenSNPs$Intogen_driver[i] <- gsub("No", "", finalIntogenSNPs$Intogen_driver[i])
  }
}

# Remove hanging spaces created by aggregate function
for(i in (1:ncol(finalIntogenSNPs))){
  finalIntogenSNPs[,i] <- gsub("\\>\\s+$", "", finalIntogenSNPs[,i])
}


### COSMIC Annotation
#Create new dataframe with relevant information from Cosmic and trim gene names with _EN...
CosmicVariants <- subset(CosmicCompleteExport, Mutation.GRCh37.genome.position != "", select = c(Mutation.GRCh37.genome.position, Pubmed_PMID, Gene.name, Mutation.CDS, ID_sample, Sample.name, Primary.site, Site.subtype, Primary.histology, Histology.subtype, Sample.source, Comments))
totalSamples <- length(unique(CosmicVariants[,"ID_sample"]))

#Remove chr from chromosome names in vcf. Convert position format for merge. Merge. Remove extra column used for merging. Double check mutation using subset
snpVCFbyPos$CosmicFormat <- paste(snpVCFbyPos$CHROM, ":", snpVCFbyPos$POS, "-", snpVCFbyPos$POS, sep="")
snpsInCosmic <- merge(snpVCFbyPos, CosmicVariants, by.x = "CosmicFormat", by.y = "Mutation.GRCh37.genome.position")
snpVCFbyPos$CosmicFormat <- NULL
snpsInCosmic$CosmicFormat <- NULL  
snpsInCosmic$altAllele1 <- sapply(strsplit(as.character(snpsInCosmic$Mutation.CDS), ">"), "[", 2)
snpsInCosmic$altAllele2 <- ifelse((snpsInCosmic$altAllele1 == "A"), "T",
                                  ifelse((snpsInCosmic$altAllele1 == "T"), "A",
                                         ifelse((snpsInCosmic$altAllele1 == "G"), "C",
                                                ifelse((snpsInCosmic$altAllele1 == "C"), "G", NA))))
snpsInCosmic <- subset(snpsInCosmic, (is.na(snpsInCosmic$altAllele1)) | (snpsInCosmic$ALT1 == snpsInCosmic$altAllele1) | (snpsInCosmic$ALT2 == snpsInCosmic$altAllele1) | (snpsInCosmic$ALT3 == snpsInCosmic$altAllele1) | (snpsInCosmic$ALT1 == snpsInCosmic$altAllele2) | (snpsInCosmic$ALT2 == snpsInCosmic$altAllele2) | (snpsInCosmic$ALT3 == snpsInCosmic$altAllele2))
snpsInCosmic$altAllele1 <- NULL
snpsInCosmic$altAllele2 <- NULL

snpsInCosmic$Histology.subtype <- paste(snpsInCosmic$Histology.subtype, snpsInCosmic$Comments, sep=" ")
snpsInCosmic$Comments <- NULL

#Aggregate duplicate positions and Format. If there were no matches in COSMIC, create fake dataframe. Also account for dataframe with one row because apply makes it into a list. 
if(nrow(snpsInCosmic) == 0){
  finalCosmicSNPs$CHROM[1] <- "fake"
  finalCosmicSNPs$POS[1] <- "fake"
  finalCosmicSNPs$rsID[1] <- "fake"
  finalCosmicSNPs$Gene[1] <- "fake"
  finalCosmicSNPs$Effect[1] <- "fake"
  finalCosmicSNPs$COSMIC_reference[1] <- "fake"
  finalCosmicSNPs$COSMIC_mutation[1] <- "fake"
  finalCosmicSNPs$COSMIC_sample[1] <- "fake"
  finalCosmicSNPs$COSMIC_histology[1] <- "fake"
  finalCosmicSNPs <- as.data.frame(finalCosmicSNPs)
}else{
  snpsInCosmicAgged <- aggregate(snpsInCosmic[,c(3:5,23:32)], snpsInCosmic[,1:2], FUN = function(x) c(count = length(unique(x[!is.na(x)])), paste(rle(sort(x))$values, "(", rle(sort(x))$lengths, "x", ")", "  ", sep="", collapse="")))
  snpsInCosmicAgged <- apply(snpsInCosmicAgged, 2, unlist)
  if(!is.na(nrow(snpsInCosmicAgged))){
    snpsInCosmicAgged <- as.data.frame(snpsInCosmicAgged[,c(1,2,4,6,8,10,12,14,15,16,18,20,22,24,26,28)])
    snpsInCosmicAgged$ID_sample.count <- (as.numeric(snpsInCosmicAgged$ID_sample.count)/totalSamples)*100
    snpsInCosmicAgged <- snpsInCosmicAgged[order(snpsInCosmicAgged$CHROM, snpsInCosmicAgged$POS),]
    snpsInCosmicAgged <- apply(snpsInCosmicAgged, 2, function(x) gsub("\\(1x\\)|\\(x\\)", "", x))
    snpsInCosmicAgged <- apply(snpsInCosmicAgged, 2, function(x) gsub("\\<NS\\>", "", x))
    snpsInCosmicAgged <- as.data.frame(snpsInCosmicAgged)
  }else{
    snpsInCosmicAgged <- as.data.frame(t(snpsInCosmicAgged[c(1,2,4,6,8,10,12,14,15,16,18,20,22,24,26,28)]))
    snpsInCosmicAgged$ID_sample.count <- (as.numeric(snpsInCosmicAgged$ID_sample.count)/totalSamples)*100
    snpsInCosmicAgged <- snpsInCosmicAgged[order(snpsInCosmicAgged$CHROM, snpsInCosmicAgged$POS),]
    snpsInCosmicAgged <- apply(snpsInCosmicAgged, 2, function(x) gsub("\\(1x\\)|\\(x\\)", "", x))
    snpsInCosmicAgged <- as.data.frame(t(snpsInCosmicAgged))
    snpsInCosmicAgged <- apply(snpsInCosmicAgged, 2, function(x) gsub("\\<NS\\>", "", x))
    snpsInCosmicAgged <- as.data.frame(t(snpsInCosmicAgged))
    snpsInCosmicAgged <- as.data.frame(snpsInCosmicAgged)
  }

# Remove repeat numbers from rsID, Gene, and Effect columns
  snpsInCosmicAgged$rsID. <- gsub("\\(\\w+\\)", "", snpsInCosmicAgged$rsID.)
  snpsInCosmicAgged$Gene. <- gsub("\\(\\w+\\)", "", snpsInCosmicAgged$Gene.)
  snpsInCosmicAgged$Effect. <- gsub("\\(\\w+\\)", "", snpsInCosmicAgged$Effect.)

# Remove repeat numbers for spaces
  for(i in (1:ncol(snpsInCosmicAgged))){
    snpsInCosmicAgged[,i] <- gsub("\\s+\\(\\w+\\)|^\\(\\w+\\)\\s+", "", snpsInCosmicAgged[,i])
  }

# Remove spaces appended by aggregate
  for(i in (1:ncol(snpsInCosmicAgged))){
    snpsInCosmicAgged[,i] <- gsub("\\>\\s+$", "", snpsInCosmicAgged[,i])
  }

# Format
  setnames(snpsInCosmicAgged, "Pubmed_PMID.", "COSMIC_reference")
  snpsInCosmicAgged$COSMIC_mutation <- paste(snpsInCosmicAgged$Gene.name, ": ", snpsInCosmicAgged$Mutation.CDS, sep="")
  snpsInCosmicAgged$COSMIC_sample <- paste("Sample: ", snpsInCosmicAgged$Sample.name, " ", snpsInCosmicAgged$Sample.source, " ID:", snpsInCosmicAgged$ID_sample., ", frequency in COSMIC: ", snpsInCosmicAgged$ID_sample.count, sep="")
  snpsInCosmicAgged$COSMIC_histology <- paste(snpsInCosmicAgged$Primary.site, snpsInCosmicAgged$Site.subtype, snpsInCosmicAgged$Primary.histology, snpsInCosmicAgged$Histology.subtype, sep=" ")
  snpsInCosmicAgged$COSMIC_sample <- gsub("\\s+", " ", snpsInCosmicAgged$COSMIC_sample)
  snpsInCosmicAgged$COSMIC_histology <- gsub("\\s+", " ", snpsInCosmicAgged$COSMIC_histology)
  setnames(snpsInCosmicAgged, "rsID.", "rsID")
  setnames(snpsInCosmicAgged, "Gene.", "Gene")
  setnames(snpsInCosmicAgged, "Effect.", "Effect")
  
  finalCosmicSNPs <- snpsInCosmicAgged[,c("CHROM", "POS", "rsID", "Gene", "Effect", "COSMIC_reference", "COSMIC_mutation", "COSMIC_sample", "COSMIC_histology")]
}


### ClinVar Annotation
# Format and merge on rsID and Start
vcfForClinvar <- snpVCFbyPos
names(vcfForClinvar)[1] <- "Chromosome"
names(vcfForClinvar)[2] <- "Start"
names(vcfForClinvar)[3] <- "rsID"
names(clinVarSummary)[7] <- "rsID"
vcfForClinvar[,3] <- gsub("rs", "", vcfForClinvar[,3])
clinVarRS <- merge(vcfForClinvar, clinVarSummary, by="rsID")
clinVarStart <- merge(vcfForClinvar, clinVarSummary, by=c("Chromosome", "Start"))

# Format and merge on Stop
names(vcfForClinvar)[2] <- "Stop"
names(clinVarSummary)[7] <- "rsID"
clinVarStop <- merge(vcfForClinvar, clinVarSummary, by=c("Chromosome", "Stop"))

# Remove nonSNPs
clinVarStart <- subset(clinVarStart, Start == Stop)
clinVarStop <- subset(clinVarStop, Start == Stop)

# Format for combining dataframes (columns: rsID, CHROM, POS, Gene, Effect, Type, Name, ClinicalSignificance, PhenotypeIDs, Origin, ReviewStatus, NumberSubmitters, OtherIDs, X.AlleleID, REF:GENE5)
names(clinVarRS)[c(2,3)] <- c("CHROM", "POS")
clinVarRS$Chromosome.y <- NULL
clinVarRS$Start.y <- NULL
clinVarRS$Stop <- NULL
clinVarRS <- clinVarRS[,c(1,2,3,4,5,24,25,28,32,33,36,39,42,23,6:22)]

names(clinVarStart)[c(1,2,3)] <- c("CHROM", "POS", "rsID")
clinVarStart$rsID.y <- NULL
clinVarStart$Stop <- NULL
clinVarStart <- clinVarStart[,c(3,1,2,4,5,24,25,28,32,33,36,39,42,23,6:22)]

names(clinVarStop)[c(1,2,3)] <- c("CHROM", "POS", "rsID")
clinVarStop$rsID.y <- NULL
clinVarStop$Start <- NULL
clinVarStop <- clinVarStop[,c(3,1,2,4,5,24,25,28,32,33,36,39,42,23,6:22)]

# Combine dataframes and Keep unique entries
clinVarAllMethods <- rbind(clinVarRS, clinVarStart, clinVarStop)
clinVarFinal <- unique(clinVarAllMethods)

# Double check correct allele by matching mutation base from vcf and ClinVar's Name field
if(nrow(clinVarFinal) == 0){
  finalClinVarSNPs$CHROM[1] <- "fake"
  finalClinVarSNPs$POS[1] <- "fake"
  finalClinVarSNPs$rsID[1] <- "fake"
  finalClinVarSNPs$Gene[1] <- "fake"
  finalClinVarSNPs$Effect[1] <- "fake"
  finalClinVarSNPs$ClinVar_Rating[1] <- "fake"
  finalClinVarSNPs$ClinVar_ClinicalSignificance[1] <- "fake"
  finalClinVarSNPs$ClinVar_PhenotypeIDs[1] <- "fake"
  finalClinVarSNPs$ClinVar_Origin[1] <- "fake"
  finalClinVarSNPs$ClinVar_Mutation[1] <- "fake"
  finalClinVarSNPs$ClinVar_ReviewStatus[1] <- "fake"
  finalClinVarSNPs$ClinVar_NumberSubmitters[1] <- "fake"
  finalClinVarSNPs$ClinVar_AlleleID[1] <- "fake"
  finalClinVarSNPs$ClinVar_OtherIDs[1] <- "fake"
  finalClinVarSNPs <- as.data.frame(finalClinVarSNPs)
}else{
  clinVarFinal$clinVarAllele1 <- NA
  for(i in 1:ncol(clinVarFinal)){
    if(grepl(">", clinVarFinal$Name[i])){
    clinVarFinal$clinVarAllele1[i] <- gsub(">|\\s.*$", "", regmatches(clinVarFinal$Name[i], regexpr(">.*$", clinVarFinal$Name[i])))
    }
  }
  clinVarFinal$clinVarAllele2 <- ifelse((clinVarFinal$clinVarAllele1 == "A"), "T",
                                    ifelse((clinVarFinal$clinVarAllele1 == "T"), "A",
                                           ifelse((clinVarFinal$clinVarAllele1 == "G"), "C",
                                                  ifelse((clinVarFinal$clinVarAllele1 == "C"), "G", NA))))
  
  clinVarFinal <- subset(clinVarFinal, (is.na(clinVarFinal$clinVarAllele1)) | (clinVarFinal$ALT1 == clinVarFinal$clinVarAllele1) | (clinVarFinal$ALT2 == clinVarFinal$clinVarAllele1) | (clinVarFinal$ALT3 == clinVarFinal$clinVarAllele1) | (clinVarFinal$ALT1 == clinVarFinal$clinVarAllele2) | (clinVarFinal$ALT2 == clinVarFinal$clinVarAllele2) | (clinVarFinal$ALT3 == clinVarFinal$clinVarAllele2))
  clinVarFinal$clinVarAllele1 <- NULL
  clinVarFinal$clinVarAllele2 <- NULL
  
  # Create rating column from number of submitters
  clinVarFinal$Rating <- ifelse(grepl("not", clinVarFinal$ReviewStatus), 0,
                                ifelse(grepl("conflicting", clinVarFinal$ReviewStatus), 0,
                                       ifelse(grepl("single", clinVarFinal$ReviewStatus), 1,
                                              ifelse(grepl("multiple", clinVarFinal$ReviewStatus), 2,
                                                     ifelse(grepl("expert", clinVarFinal$ReviewStatus), 3,
                                                            ifelse(grepl("professional", clinVarFinal$ReviewStatus), 4, NA  ))))))
  
  # Format (columns: CHROM, POS, rsID, Gene, Effect, Rating, ClinicalSignificance, PhenotypeIDs, Origin, Name, Type, ReviewStatus, NumberSubmitters, X.AlleleID, OtherIDs
  finalClinVarSNPs <- clinVarFinal[,c(2,3,1,4,5,32,8,9,10,7,6,11,12,14,13)]
  setnames(finalClinVarSNPs, "X.AlleleID", "AlleleID")
  for(i in 1:nrow(finalClinVarSNPs)){
    if(finalClinVarSNPs$rsID[i] != -1){
      finalClinVarSNPs$rsID[i] <- paste("rs", finalClinVarSNPs$rsID[i], sep="")
    }
  }
  finalClinVarSNPs$Mutation <- paste(finalClinVarSNPs$Type, finalClinVarSNPs$Name, sep=" ")
  finalClinVarSNPs$Name <- NULL
  finalClinVarSNPs$Type <- NULL
  names(finalClinVarSNPs)[c(6:14)] <- paste("ClinVar_", names(finalClinVarSNPs)[c(6:14)], sep="")
  finalClinVarSNPs <- finalClinVarSNPs[,c(1:9, 14, 10:13)]
}


### pharmGkb Annotation
# Merge based on rsID
pharmVar <- merge(snpVCFbyPos, pharmgkbRSID, by.x = "rsID", by.y = "Variant")

# Format pharmGkb individual SNP data
names(pharmgkbAllele) <- c("rsID", "genotype", "description", "notes1", "notes2")
pharmgkbAllele$allele1 <- substr(pharmgkbAllele$genotype, 1, 1)
pharmgkbAllele$allele2 <- substr(pharmgkbAllele$genotype, 2, 2)
  
# Add individual SNP annotation from pharmGkb
pharmVar[,28] <- list("")
pharmVar[,29] <- list("")
pharmVar[,30] <- list("")
  
names(pharmVar)[28] <- "description"
names(pharmVar)[29] <- "notes1"
names(pharmVar)[30] <- "notes2"
  
if(nrow(pharmVar)!=0) {
  for(i in 1:(nrow(pharmVar))){
    for(j in 1:(nrow(pharmgkbAllele))){
      if((pharmVar[i,1] == pharmgkbAllele[j,1]) & (((pharmVar[i,16] == pharmgkbAllele[j,6]) & (pharmVar[i,17] == pharmgkbAllele[j,7])) | ((pharmVar[i,16] == pharmgkbAllele[j,7]) & (pharmVar[i,17] == pharmgkbAllele[j,6])))){
        pharmVar[i,28] <- paste(pharmVar[i,28], pharmgkbAllele[j,3], sep = " ")
        pharmVar[i,29] <- paste(pharmVar[i,29], pharmgkbAllele[j,4], sep = " ")
        pharmVar[i,30] <- paste(pharmVar[i,30], pharmgkbAllele[j,5], sep = " ")
      }
    }
  }
}
  
# Format
finalPharmgkbSNPs <- pharmVar[,c(2,3,1,4,5,24,26,25,27,23,28:30)]
setnames(finalPharmgkbSNPs, "Gene.x", "Gene")
setnames(finalPharmgkbSNPs, "Type", "Reaction")
setnames(finalPharmgkbSNPs, "Strength.of.evidence..level.", "Evidence.Level")
names(finalPharmgkbSNPs)[c(6:13)] <- paste("PharmGkb_", names(finalPharmgkbSNPs)[c(6:13)], sep="")
setnames(finalPharmgkbSNPs, "PharmGkb_Gene.y", "PharmGkb_Gene")

if(nrow(finalPharmgkbSNPs)==0){
  for(i in 1:ncol(finalPharmgkbSNPs)){
    finalPharmgkbSNPs[1,i] <- "fake"
  }
}


### DrugBank Annotation
# Merge on snp rsID
snpsInDrugBank <- merge(snpVCFbyPos, drugBankTable, by.x="rsID", by.y="snp")

# Pull out allele (three possibilities: G>A; A allele; A Allele) and subset on it
snpsInDrugBank$altAllele1 <- sapply(strsplit(as.character(snpsInDrugBank$allele.change), "> "), "[", 2)
snpsInDrugBank$altAllele2 <- sapply(strsplit(as.character(snpsInDrugBank$allele.change), "\\sall|\\sAll"), "[", 1)

# Double check that alternative alleles match
snpsInDrugBank1 <- subset(snpsInDrugBank, (snpsInDrugBank$ALT1 == snpsInDrugBank$altAllele1) | (snpsInDrugBank$ALT2 == snpsInDrugBank$altAllele1) | (snpsInDrugBank$ALT3 == snpsInDrugBank$altAllele1))
snpsInDrugBank2 <- subset(snpsInDrugBank, (snpsInDrugBank$ALT1 == snpsInDrugBank$altAllele2) | (snpsInDrugBank$ALT2 == snpsInDrugBank$altAllele2) | (snpsInDrugBank$ALT3 == snpsInDrugBank$altAllele2))
snpsInDrugBank <- rbind(snpsInDrugBank1, snpsInDrugBank2)

snpsInDrugBank$altAllele1 <- NULL
snpsInDrugBank$altAllele2 <- NULL

# Format final table (column order: CHROM, POS, rsID, Gene, Effect, drug.name, reaction, reference, drug.number, allele.name, allele.change, gene.name, gene.symbol, uniprot)
finalDrugBankSNPs <- snpsInDrugBank[,c(2,3,1,4,5,27,25,26,28,23,24,29:31)]
setnames(finalDrugBankSNPs, "drug.name", "drug")
names(finalDrugBankSNPs)[c(6:14)] <- paste("DrugBank_", names(finalDrugBankSNPs)[c(6:14)], sep="")

if(nrow(finalDrugBankSNPs)==0){
  for(i in 1:ncol(finalDrugBankSNPs)){
    finalDrugBankSNPs[1,i] <- "fake"
  }
}


### Combine and Format Annotations
library(plyr)

# Cast as data.frames then rbind.fill
finalIntogenSNPs <- as.data.frame(finalIntogenSNPs)
finalCosmicSNPs <- as.data.frame(finalCosmicSNPs)
finalClinVarSNPs <- as.data.frame(finalClinVarSNPs)
finalPharmgkbSNPs <- as.data.frame(finalPharmgkbSNPs)
finalDrugBankSNPs <- as.data.frame(finalDrugBankSNPs)

allAnnos <- rbind.fill(finalIntogenSNPs, finalCosmicSNPs, finalClinVarSNPs, finalPharmgkbSNPs, finalDrugBankSNPs)

# Aggregate on chromosome and position
allAnnos$Effect <- gsub("\\>;\\s+$", ";", allAnnos$Effect)
allAnnos <- aggregate(allAnnos[,c(3:54)], allAnnos[,c(1,2)], function(x) paste(unique(x[!is.na(x)]), " ", sep="", collapse=""))

# Remove spaces appended by aggregate function
for(i in (1:ncol(allAnnos))){
  allAnnos[,i] <- gsub("\\>\\s+$", "", allAnnos[,i])
}

allAnnos <- allAnnos[-which(allAnnos$CHROM=="fake"),]


### Label Differentially Expressed (Yes or No)
# Import list of differentially expressed genes
#diffGenes <- names(diffExpNBinom$FCdiff)

## Label differentially expressed genes
# Create DIFF column
#allAnnos$DIFF <- list("")

# Remove spaces in gene columns
#allAnnos$GENE1 <- gsub("\\s+", "", allAnnos$GENE1)
#allAnnos$GENE2 <- gsub("\\s+", "", allAnnos$GENE2)
#allAnnos$GENE3 <- gsub("\\s+", "", allAnnos$GENE3)
#allAnnos$GENE4 <- gsub("\\s+", "", allAnnos$GENE4)
#allAnnos$GENE5 <- gsub("\\s+", "", allAnnos$GENE5)

#for(i in 1:(nrow(allAnnos))){
#  if((allAnnos$GENE1[i] %in% diffGenes)){
#    allAnnos$DIFF[i] <- paste(allAnnos$DIFF[i], "Yes", "(", allAnnos$GENE1[i], ")", sep="")
#  }else if(allAnnos$GENE2[i] %in% diffGenes){
#    allAnnos$DIFF[i] <- paste(allAnnos$DIFF[i], "Yes", "(", allAnnos$GENE2[i], ")", sep="")
#  }else if(allAnnos$GENE3[i] %in% diffGenes){
#    allAnnos$DIFF[i] <- paste(allAnnos$DIFF[i], "Yes", "(", allAnnos$GENE3[i], ")", sep="")
#  }else if(allAnnos$GENE4[i] %in% diffGenes){
#    allAnnos$DIFF[i] <- paste(allAnnos$DIFF[i], "Yes", "(", allAnnos$GENE4[i], ")", sep="")
#  }else if(allAnnos$GENE5[i] %in% diffGenes){
#    allAnnos$DIFF[i] <- paste(allAnnos$DIFF[i], "Yes", "(", allAnnos$GENE5[i], ")", sep="")
#  }else{
#    allAnnos$DIFF[i] <- "No"
#  }
#}

#allAnnos$DIFF <- unlist(allAnnos$DIFF)


### CADD Annotation
# make temporary files for input and output
tempRegions <- tempfile("cadd_regions_", fileext=".txt")
tempCADD <- tempfile("cadd_scores_", fileext=".vcf")

# create tabix region input file
positions <- paste(snpVCFbyPos$CHROM, ":", snpVCFbyPos$POS, "-", snpVCFbyPos$POS, sep="")
write.table(positions, tempRegions, row.names=FALSE, col.names=FALSE)

# extract relevant lines from CADD database
system(paste("xargs -a ", tempRegions ," -I {} tabix ", config$CADD, " {} > ", tempCADD, sep=""))
caddVCF <- read.table(tempCADD, sep = "\t", stringsAsFactors=FALSE)

# merge to add CADD scores to allAnnos table
names(caddVCF) <- c("CHROM", "POS", "REF", "ALT", "CADD_raw", "CADD_PHRED")
allAnnos <- merge(allAnnos, caddVCF, by=c("CHROM", "POS", "REF", "ALT"))


### Create final table
#finalSNPs <- as.data.frame(allAnnos[,c("CHROM", "POS", "rsID", "Gene", "DIFF", "Intogen_driver", "CADD_PHRED", "Effect", "ClinVar_Rating", "ClinVar_ClinicalSignificance", "COSMIC_histology", "COSMIC_reference", "PharmGkb_Reaction", "PharmGkb_Drugs", "DrugBank_drug", "DrugBank_reaction", "DrugBank_reference")])
finalSNPs <- as.data.frame(allAnnos[,c("CHROM", "POS", "rsID", "Gene", "Intogen_driver", "Intogen_info", "CADD_PHRED","Effect", "ClinVar_Rating", "ClinVar_ClinicalSignificance", "COSMIC_histology", "COSMIC_reference", "PharmGkb_Reaction", "PharmGkb_Drugs", "DrugBank_drug", "DrugBank_reaction", "DrugBank_reference")])

finalSNPs <- finalSNPs[order(finalSNPs$CHROM, finalSNPs$POS),]

setnames(finalSNPs, "Effect", "SnpEff_Effect")
#setnames(finalSNPs, "DIFF", "Diff_Exp")

# export as #csv
write.table(finalSNPs, paste(config$REPORT_RESULTS, '/', patientID, '/snps.csv', sep=''), row.names=FALSE, sep='\t')

### enter html links
finalSNPs <- transform(finalSNPs,
                       rsID = paste('<a href = ', 
                                    shQuote(paste('http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?',
                                                  gsub('rs', 'rs=', finalSNPs$rsID), 
                                                  sep='')), 
                                    'target="_blank" >', 
                                    finalSNPs$rsID, 
                                    '</a>'))
#finalSNPs <- transform(finalSNPs, 
#                       Gene = 
#)
finalSNPs <- transform(finalSNPs,
                       COSMIC_reference = paste('<a href = ',
                                                shQuote(paste('http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=',
                                                              finalSNPs$COSMIC_reference, 
                                                              sep='')), 
                                                'target="_blank" >', 
                                                finalSNPs$COSMIC_reference, 
                                                '</a>'))
#finalSNPs <- transform(finalSNPs, 
#                       DrugBank_reference = 
#)


### Knitr
#kable(snp_eff, format='html', table.attr = 'id=\"SNPs_SNPeff_table\"', row.names=F)
#kable(annovar, format='html', table.attr = 'id=\"SNPs_ANNOVAR_table\"', row.names=F)
kable(finalSNPs, format='html', table.attr = 'id=\"SNPs_table\"', row.names=F)


}
snpAnnotation_DNA()
