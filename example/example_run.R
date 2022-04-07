
library(signature.tools.lib)
#library(MMRDetect)

source("~/Documents/Github/MMRDetect/R/Gen_catalogues.R")
load("~/Documents/Github/MMRDetect/data/indelsig_template.rda")

source("~/Documents/Github/MMRDetect/R/MMRDetect.compute.variables.R")
load("~/Documents/Github/MMRDetect/data/PancanSig.rda")
load("~/Documents/Github/MMRDetect/data/MMRKO_subsig.rda")
load("~/Documents/Github/MMRDetect/data/MMRKO_indelsig.rda")

source("~/Documents/Github/MMRDetect/R/MMRDetect.classify.R")
load("~/Documents/Github/MMRDetect/data/MMRDclassifier.rda")

setwd('/Volumes/')

args = commandArgs(trailingOnly=TRUE)

sample_name <- args[1]
tissue <- args[2]
snv_vcf_fileLoc  <- args[3]
indel_vcf_fileLoc <- args[4]

# Rearrange SNV catalog
snv_catalogue <- signature.tools.lib::vcfToSNVcatalogue(vcfFilename = snv_vcf_fileLoc, genome.v = "hg38")

snv_catalogue_mutType <- snv_catalogue$catalogue
names(snv_catalogue_mutType)[1] <- sample_name
snv_catalogue_mutType$MutationType <- row.names(snv_catalogue_mutType)

# Generate MSI-indel catalogue
all_indels <- read.table(indel_vcf_fileLoc,sep = "\t", header = F, as.is = T,)[,c(1,2,4,5)]
names(all_indels) <- c("Chrom", "Pos", "Ref", "Alt")
all_indels$Sample <- sample_name

indel_classified <- MSI_indel_classifier(indels = all_indels, genome.v = "hg38")

indel_catalogues <- gen_indelmuttype_MMRD(indel_classified)

# Apply MMRDetect
snv_catalogue_mutType$dummy <- 1
indel_catalogues$dummy <- 1

muts_variables <- MMRDetect.compute.variables(sub_cat=snv_catalogue_mutType, 
                                              indel_cat=indel_catalogues, 
                                              tissue_type=tissue )

muts_variables <- muts_variables[muts_variables$Sample != "dummy",]

muts_variables$MMR_sum <- muts_variables$MMR_sum / sum(snv_catalogue$catalogue)
muts_variables$RepIndel_num <- muts_variables$RepIndel_num / sum(indel_catalogues[,3])

MMRDetect_classified <- MMRDetect.classify(muts_variables)

MMRDetect_classified$MMRDcall[MMRDetect_classified$glm_prob > 0.90] <- "Non-MSI-H"
MMRDetect_classified$MMRDcall[MMRDetect_classified$glm_prob <= 0.90] <- "MSI-H"

write.table(
MMRDetect_classified,
file = paste("~/Documents/Github/MMRDetect/example/",sample_name,".MMRD.txt",sep=""),
append = F, quote = FALSE, sep = "\t", 
eol = "\n", na = "NA",dec = ".", row.names = FALSE, 
col.names = TRUE
)

