library(tidyverse)

# ENVIRONMENT
setwd("/lisc/scratch/botany/frschmidt/mk_test/")

# FUNCTIONS

# VARIABLES

vcf_path <- "raw_data/snpeff_annotation/pair_cal_spn/vieref_cal_spn.ann.miss-syno.vcf.gz"
vcf_file <- readLines(vcf_path)

vcf_SNPs <- data.frame(vcf_file[grep(pattern="#CHROM",vcf_file):length(vcf_file)])
vcf_SNPs <- data.frame(do.call('rbind', strsplit(as.character(vcf_SNPs[,1]),'\t',fixed=TRUE)))
colnames(vcf_SNPs) <- as.character(unlist(vcf_SNPs[1,]))
vcf_SNPs <- vcf_SNPs[-1,]
colnames(vcf_SNPs)[1] <- "scaffold"

# read in allele frequencies (the output from vcftools is a bit unpractical, so we split the columns)
outgroup_freq_path <- "raw_data/pop_allel_freq/pair_cal_spn/calciphila.frq"
outgroup_freq <- read.delim(outgroup_freq_path, header=FALSE, skip=1)
colnames(outgroup_freq) <- c("scaffold","pos","n_alleles","n_chromosomes","ref","alt")
outgroup_freq <- separate(outgroup_freq,ref,c("ref_allele","ref_freq"), ":")
outgroup_freq <- separate(outgroup_freq,alt,c("alt_allele","alt_freq"), ":")
outgroup_freq$ref_freq<-as.numeric(outgroup_freq$ref_freq)
outgroup_freq$alt_freq<-as.numeric(outgroup_freq$alt_freq)

target_freq_path <- "raw_data/pop_allel_freq/pair_cal_spn/sppicnga.frq"
target_freq <- read.delim(target_freq_path, header=FALSE, skip=1)
colnames(target_freq) <- c("scaffold","pos","n_alleles","n_chromosomes","ref","alt")
target_freq <- separate(target_freq,ref,c("ref_allele","ref_freq"), ":") # separate columns into new columns based on delimiter ":"
target_freq <- separate(target_freq,alt,c("alt_allele","alt_freq"), ":")
target_freq$ref_freq <- as.numeric(target_freq$ref_freq)
target_freq$alt_freq <- as.numeric(target_freq$alt_freq)

# attach allele frequencies to SNP file
vcf_SNPs$outgroup_ref_freq <- outgroup_freq$ref_freq
vcf_SNPs$outgroup_alt_freq <- outgroup_freq$alt_freq

vcf_SNPs$target_ref_freq <- target_freq$ref_freq
vcf_SNPs$target_alt_freq <- target_freq$alt_freq

# only keep SNPs that are either synonymous or nonsynonymous.
#vcf_SNPs<-vcf_SNPs[grep(pattern="synonymous_variant|missense_variant",vcf_SNPs$INFO),]
vcf_SNPs$annotation <- separate(vcf_SNPs,INFO,c("A","B","C","D","E","F","G","H","J","K","L","Z","X","N","V","Q","W","R"), "\\|")[,"B"] # random number of columns with random column names, needs to be more than the maximum amount of columns possible, you only need the second column
vcf_SNPs$gene <- separate(vcf_SNPs,INFO,c("A","B","C","D","E","F","G","H","J","K","L","Z","X","N","V","Q","W","R"), "\\|")[,"D"]
vcf_SNPs$gene <- as.factor(vcf_SNPs$gene)
vcf_SNPs$scaffold <- factor(vcf_SNPs$scaffold) # reset factor levels

## count pN, pS, dN, dS

# create columns with 0s
vcf_SNPs$pN<-0
vcf_SNPs$pS<-0
vcf_SNPs$dN<-0
vcf_SNPs$dS<-0


find_div_pol_sites <- function(out_ref, out_alt, tar_ref, tar_alt, site_ann) {
  # result vector where each position corresponds to
  # 1: dN, 2: dS, 3: pN, 4: pS
  res_vec <- rep(0, 4) 
  if (is.na(out_ref) || is.na(out_alt) || is.na(tar_ref) || is.na(tar_alt)) {
    return(res_vec)
  }
  if ((out_ref == 0 || out_alt == 0) && (tar_ref == 0 || tar_alt == 0)) { 
    if ((out_ref == 0 && tar_alt == 0) || (out_alt == 0 && tar_ref == 0)) { 
      if ("missense_variant" %in% site_ann) { res_vec[1] <- 1 }
      else if ("synonymous_variant" %in% site_ann) { res_vec[2] <- 1 }
    } 
  }
  if ((out_ref == 0 || out_alt == 0) && (tar_ref != 0 && tar_alt != 0)) {
    if ("missense_variant" %in% site_ann) { res_vec[3] <- 1 }
    else if ("synonymous_variant" %in% site_ann) { res_vec[4] <- 1 }
  }
  return(res_vec)
}

# Testing function
popu <- c(1, 0, 1, 0) # e.g. REF: C, OUT: C, TAR: C
popu <- c(0.6, 0.4, 0, 1) # e.g. REF: A, OUT: a, TAR: G
popu <- c(0, 1, 0.6, 0.4) # e.g. REF: G, OUT: T, TAR: t
popu <- c(1, 0, 0.6, 0.4) # e.g. REF: T, OUT: T, TAR: g
popu <- c(0, 1, 1, 0) # e.g. REF: C, OUT: A, TAR: C
popu <- c(0, 1, 0, 1) # e.g. REF: G, OUT: T, TAR: T
popu <- c(0, 1, 1, NA)
find_div_pol_sites(popu[1], popu[2], popu[3], popu[4], "missense_variant")

vcf_SNPs$sites <- matrix(0, ncol = 4, nrow = nrow(vcf_SNPs))

# check the category a SNP falls into
for(i in 1:nrow(vcf_SNPs)) {
  vcf_SNPs$sites[i, ] <- find_div_pol_sites(
    vcf_SNPs$outgroup_ref_freq[i],
    vcf_SNPs$outgroup_alt_freq[i],
    vcf_SNPs$target_ref_freq[i],
    vcf_SNPs$target_alt_freq[i],
    vcf_SNPs$annotation[i]
  )
}

save(vcf_SNPs, file = "results/mk_test/vcf_dataframe.Rdata")

# make empty data.frame
MKT<-data.frame(gene=factor(),pN=numeric(), pS=numeric(), dN=numeric(), dS=numeric())

#sum pNs, pSs, dNs, and dSs for each contig across SNPs
for( i in 1:length(levels(vcf_SNPs$gene))) {
  temp <- vcf_SNPs[which(vcf_SNPs$gene == levels(vcf_SNPs$gene)[i]), ]
  MKT<-rbind(MKT,data.frame(
    gene = as.character(temp$gene[1]),
    dN = sum(temp$sites[, 1]),
    dS = sum(temp$sites[, 2]),
    pN = sum(temp$sites[, 3]),
    pS = sum(temp$sites[, 4])))
}

# perform MKT test
MKT$pN.pS=MKT$pN/MKT$pS # calculate ration of nonsynonymous to synonymous polymorphisms
MKT$dN.dS=MKT$dN/MKT$dS # calculate ration of nonsynonymous to synonymous substitutions
MKT$fisher.test.P<-99  # create new column for p-values
for(i in 1:nrow(MKT)){
  MKT$fisher.test.P[i]<-fisher.test(matrix(as.numeric(MKT[i,c(3,2,5,4)]), ncol=2))$p.value # calculate fisher exact test and copy p-value for every contig
  if((MKT$pN[i] == 0 && MKT$dN[i] == 0) || (MKT$pS[i] == 0 && MKT$dS[i] == 0) || (MKT$pS[i] == 0 && MKT$pN[i] == 0) || (MKT$dS[i] == 0 && MKT$dN[i] == 0)) { MKT$fisher.test.P[i]<-NA} # this lines assigns an NA to all p-values that are meaningless, because the contingency table was incomplete
  if(sum(as.numeric(MKT[i,c(3,2,5,4)])) < 3) { MKT$fisher.test.P[i]<-NA } # only use cases where total number of SNPs is higher than or equal to 3
}

save(MKT, file = "results/mk_test/mkt.Rdata")
load(file = "results/mk_test/mkt.Rdata")

MKT_noNAs<-MKT[which(MKT$fisher.test.P != "NA"),] # remove all cases where fisher exact test is meaningless


# multiple hypothesis testing. It seems to be uncommon with MKT; probably because fisher exact produces low P-values only with much higher counts than are common for SNP data
MKT_noNAs$fisher.test.P<-p.adjust(MKT_noNAs$fisher.test.P, method = "BH") # correct p-values using Benjamini & Hochberg (1985) FDR