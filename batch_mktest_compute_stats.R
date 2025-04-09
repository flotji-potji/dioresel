#!/usr/bin/env R
#
### ENVIRONMENT
library(tidyverse)

options(error = quote({
  dump.frames(to.file=T, dumpto='last.dump')
  load('last.dump.rda')
  print(last.dump)
  q()
}))

### FUNCTIONS
get_vcf_df <- function(vcf.path) {
  vcf.file <- readLines(vcf.path)
  vcf.df <- data.frame(
    vcf.file[grep(pattern = "#CHROM", vcf.file):length(vcf.file)]
  )
  vcf.df <- data.frame(
    do.call('rbind', strsplit(as.character(vcf.df[,1]),'\t',fixed=TRUE))
  )
  colnames(vcf.df) <- as.character(unlist(vcf.df[1, ]))
  vcf.df <- vcf.df[-1,]
  colnames(vcf.df)[1] <- "scaffold"
  vcf.df$scaffold <- factor(vcf.df$scaffold) # reset factor levels
  info_fields <- separate(vcf.df, INFO, as.character(1:18), "\\|")
  vcf.df$annotation <- info_fields$"2"
  vcf.df$gene <- info_fields$"4"
  vcf.df$gene <- as.factor(vcf.df$gene)
  return(vcf.df)
}

get_pop_freq_df <- function(pop_freq.path) {
  pop_freq.df <- read.delim(pop_freq.path, header=FALSE, skip=1)
  colnames(pop_freq.df) <- c("scaffold", "pos", "n_alleles",
                              "n_chromosomes","ref","alt")
  pop_freq.df <- separate(pop_freq.df, ref, c("ref_allele","ref_freq"), ":")
  pop_freq.df <- separate(pop_freq.df, alt, c("alt_allele","alt_freq"), ":")
  pop_freq.df$ref_freq <- as.numeric(pop_freq.df$ref_freq)
  pop_freq.df$alt_freq <- as.numeric(pop_freq.df$alt_freq)
  return(pop_freq.df)
}

find_div_pol_sites <- function(out_ref, out_alt, tar_ref, tar_alt, site_ann) {
  # result vector representing following conditions
  # 1: dN-outgroup, 2: dS-outgroup, 3: dN-target, 4: dS-target
  # 5: dN-poly-out, 6: dS-poly-out, 7: dN-poly-target, 8: dS-poly-target
  # 9: pN-fix-out, 10: pS-fix-out, 11: pN-poly-out, 12: pS-poly-out
  # whereas d*-outgroup = divergence fixed in outgroup alt allele
  # d*-target = divergence fixed in target alt allele
  # d*-poly-out = divergence fixed in ref target but polymorphic in outgroup
  # d*-poly-tar = divergence fixed in alt target but polymorphic in outgroup
  # p*-fix-out = polymorphism in target and fixed in outgroup
  # p*-poly-out = polymorphisms in target and polymorphisms in outgroup
  res_vec <- rep(0, 12)
  # ensure that all frequency inputs are not NA
  if (is.na(out_ref) || is.na(out_alt) || is.na(tar_ref) || is.na(tar_alt)) {
    return(list(
      dN.out = res_vec[1], dS.out = res_vec[2],
      dN.tar = res_vec[3], dS.tar = res_vec[4],
      dN.poly.out = res_vec[5], dS.poly.out = res_vec[6],
      dN.poly.tar = res_vec[7], dS.poly.tar = res_vec[8],
      pN.fix.out = res_vec[9], pS.fix.out = res_vec[10],
      pN.poly.out = res_vec[11], pS.poly.out = res_vec[12]
    ))
  }
  # d*-outgroup and d*-target
  if ((tar_ref == 1 && tar_alt == 0) || (tar_ref == 0 && tar_alt == 1)) {
    # d*-outgroup
    if (out_alt == 1 && tar_alt == 0) {
      if ("missense_variant" == site_ann) {
        res_vec[1] <- 1
      } else if ("synonymous_variant" == site_ann) {
        res_vec[2] <- 1
      }
    }
    # d*-target
    if (out_alt == 0 && tar_alt == 1) {
      if ("missense_variant" == site_ann) {
        res_vec[3] <- 1
      } else if ("synonymous_variant" == site_ann) {
        res_vec[4] <- 1
      }
    }
    # d*-poly-out and d*-poly-tar
    # d*-poly-out
    if ((out_ref > 0 && out_ref < 1) && (out_alt > 0 && out_alt < 1)) {
      if (tar_ref == 1) {
        if ("missense_variant" == site_ann) {
          res_vec[5] <- 1
        } else if ("synonymous_variant" == site_ann) {
          res_vec[6] <- 1
        }
      } else { # d*-poly-tar
        if ("missense_variant" == site_ann) {
          res_vec[7] <- 1
        } else if ("synonymous_variant" == site_ann) {
          res_vec[8] <- 1
        }
      }
    }
  # p*-fix-out and p*-poly-out
  } else if ((tar_alt > 0 && tar_alt < 1) && (tar_ref > 0 && tar_ref < 1)) {
    # p*-fix-out
    if ((out_alt == 1 && out_ref == 0) || (out_alt == 0 || out_ref == 1)) {
      if ("missense_variant" == site_ann) {
        res_vec[9] <- 1
      } else if ("synonymous_variant" == site_ann) {
        res_vec[10] <- 1
      }
    } else { # p*-poly-out
      if ("missense_variant" == site_ann) {
        res_vec[11] <- 1
      } else if ("synonymous_variant" == site_ann) {
        res_vec[12] <- 1
      }
    }
  }
  # return values contained in result vector as named list
  return(list(
    dN.out = res_vec[1], dS.out = res_vec[2],
    dN.tar = res_vec[3], dS.tar = res_vec[4],
    dN.poly.out = res_vec[5], dS.poly.out = res_vec[6],
    dN.poly.tar = res_vec[7], dS.poly.tar = res_vec[8],
    pN.fix.out = res_vec[9], pS.fix.out = res_vec[10],
    pN.poly.out = res_vec[11], pS.poly.out = res_vec[12]
  ))
}

### VARIABLES
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Requires to supply one text file")
}

vcf.path <- args[1]
outgroup_freq.path <- args[2]
target_freq.path <- args[3]

output_dir <- "raw_data/mkt_summary/"
vcf.df.output <- paste0(output_dir, "vcf.df.Rdata")
gene.summary.output <- paste0(output_dir, "gene.summary.Rdata")

### EXECUTION
setwd("/lisc/scratch/botany/frschmidt/mk_test")
# extract vcf data frame and allele frequencies
print("[*] Start loading vcf dataframe...")
vcf.df <- get_vcf_df(vcf.path)
print("[+] Finished loading")
print("[*] Start loading outgroup dataframe...")
outgroup_freq.df <- get_pop_freq_df(outgroup_freq.path)
print("[+] Finished loading")
print("[*] Start loading target dataframe...")
target_freq.df <- get_pop_freq_df(target_freq.path)
print("[+] Finished loading")
# attach allele frequencies to SNP file
vcf.df$outgroup_ref_freq <- outgroup_freq.df$ref_freq
vcf.df$outgroup_alt_freq <- outgroup_freq.df$alt_freq
vcf.df$target_ref_freq <- target_freq.df$ref_freq
vcf.df$target_alt_freq <- target_freq.df$alt_freq

temp.vcf.df <- data.frame(dN.out = 0,
dS.out = 0,
dN.tar = 0,
dS.tar = 0,
dN.poly.out = 0,
dS.poly.out = 0,
dN.poly.tar = 0,
dS.poly.tar = 0,
pN.fix.out = 0,
pS.fix.out = 0,
pN.poly.out = 0,
pS.poly.out = 0)
#save(vcf.df, file="vcf.Rdata")
# check the category a SNP falls into
print("[*] Start identifying SNP sites...")
for (i in 1:(nrow(vcf.df))) {
  if (i %% round((nrow(vcf.df) / 8)) == 0) {
    sprintf("%s percent done @ %s", i / nrow(vcf.df) * 100, date())
  }
  sprintf("%s %s %s %s", vcf.df$outgroup_ref_freq[i],
    vcf.df$outgroup_alt_freq[i],
    vcf.df$target_ref_freq[i],
    vcf.df$target_alt_freq[i])
  temp.vcf.df[i, ] <- find_div_pol_sites(
    vcf.df$outgroup_ref_freq[i],
    vcf.df$outgroup_alt_freq[i],
    vcf.df$target_ref_freq[i],
    vcf.df$target_alt_freq[i],
    vcf.df$annotation[i]
  )
}
vcf.df <- cbind(vcf.df, temp.vcf.df)

print("[+] Finished process")

ifelse(!dir.exists(file.path(output_dir)),
        dir.create(file.path(output_dir)),
        "Directory Exists")
save(vcf.df, file = vcf.df.output)

gene.summary <- data.frame(
  gene=factor(),
  dN.out=numeric(),
  dS.out=numeric(),
  dN.tar=numeric(),
  dS.tar=numeric(),
  dN.poly.out=numeric(),
  dS.poly.out=numeric(),
  dN.poly.tar=numeric(),
  dS.poly.tar=numeric(),
  pN.fix.out =numeric(),
  pS.fix.out =numeric(),
  pN.poly.out=numeric(),
  pS.poly.out=numeric()
)

print("[*] Started calculating gene summaries...")
for (i in 1:(length(levels(vcf.df$gene)))) {
  temp <- vcf.df[which(vcf.df$gene == levels(vcf.df$gene)[i]), ]
  gene.summary <- rbind(gene.summary, data.frame(
    gene = as.character(temp$gene[1]),
    dN.out  = sum(temp$dN.out ),
    dS.out  = sum(temp$dS.out ),
    dN.tar  = sum(temp$dN.tar ),
    dS.tar  = sum(temp$dS.tar ),
    dN.poly.out = sum(temp$dN.poly.out),
    dS.poly.out = sum(temp$dS.poly.out),
    dN.poly.tar = sum(temp$dN.poly.tar),
    dS.poly.tar = sum(temp$dS.poly.tar),
    pN.fix.out  = sum(temp$pN.fix.out ),
    pS.fix.out  = sum(temp$pS.fix.out ),
    pN.poly.out = sum(temp$pN.poly.out),
    pS.poly.out = sum(temp$pS.poly.out)
  ))
}
print("[+] Finished process")

save(gene.summary, file = gene.summary.output) 
