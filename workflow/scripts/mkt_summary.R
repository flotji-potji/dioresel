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
    res_vec <- rep(0, 4)
    # ensure that all frequency inputs are not NA
    if (is.na(out_ref) || is.na(out_alt) || is.na(tar_ref) || is.na(tar_alt)) {
        return(list(
        dN = res_vec[1],
        dS = res_vec[2],
        pN = res_vec[3],
        pS = res_vec[4]
        ))
    }
    # d*-outgroup and d*-target
    if ((tar_ref == 1 && tar_alt == 0) || (tar_ref == 0 && tar_alt == 1)) {
        # d*-outgroup
        if ((out_alt == 1 && tar_alt == 0) || (out_alt == 0 && tar_alt == 1)) {
            if ("missense_variant" == site_ann) {
                res_vec[1] <- 1
            } else if ("synonymous_variant" == site_ann) {
                res_vec[2] <- 1
            }
        }
    }
    if ((tar_alt > 0 && tar_alt < 1) || (tar_ref > 0 && tar_ref < 1)) {
        if ("missense_variant" == site_ann) {
            res_vec[3] <- 1
        } else if ("synonymous_variant" == site_ann) {
            res_vec[4] <- 1
        }
    }
    return(list(
        dN = res_vec[1],
        dS = res_vec[2],
        pN = res_vec[3],
        pS = res_vec[4]
    ))
}

### VARIABLES
vcf.path <- snakemake@input$vcf
outgroup_freq.path <- snakemake@input$sample1
target_freq.path <- snakemake@input$sample2

vcf.df.output <- snakemake@output$vcf_mkt
gene.summary.output <- snakemake@output$gene_summary

output_prefix <- snakemake@params$output_prefix

### EXECUTION
# extract vcf data frame and allele frequencies
vcf.df <- get_vcf_df(vcf.path)
outgroup_freq.df <- get_pop_freq_df(outgroup_freq.path)
target_freq.df <- get_pop_freq_df(target_freq.path)
# attach allele frequencies to SNP file
vcf.df$outgroup_ref_freq <- outgroup_freq.df$ref_freq
vcf.df$outgroup_alt_freq <- outgroup_freq.df$alt_freq
vcf.df$target_ref_freq <- target_freq.df$ref_freq
vcf.df$target_alt_freq <- target_freq.df$alt_freq

temp.vcf.df <- data.frame(
dN = 0,
dS = 0,
pN = 0,
pS = 0)
# check the category a SNP falls into
for (i in 1:(nrow(vcf.df))) {
  temp.vcf.df[i, ] <- find_div_pol_sites(
    vcf.df$outgroup_ref_freq[i],
    vcf.df$outgroup_alt_freq[i],
    vcf.df$target_ref_freq[i],
    vcf.df$target_alt_freq[i],
    vcf.df$annotation[i]
  )
}
vcf.df <- cbind(vcf.df, temp.vcf.df)

gene.summary <- data.frame(
  gene=factor(),
  dN=numeric(),
  dS=numeric(),
  pN=numeric(),
  pS=numeric()
)

for (i in 1:(length(levels(vcf.df$gene)))) {
  temp <- vcf.df[which(vcf.df$gene == levels(vcf.df$gene)[i]), ]
  gene.summary <- rbind(gene.summary, data.frame(
    gene = as.character(temp$gene[1]),
    dN = sum(temp$dN),
    dS = sum(temp$dS),
    pN = sum(temp$pN),
    pS = sum(temp$pS)
  ))
}
print("[+] Finished process")

# save gene summary data frame and vcf as R objects
ifelse(!dir.exists(file.path(output_prefix)),
        dir.create(file.path(output_prefix)),
        "Directory Exists")

save(vcf.df, file = vcf.df.output)

save(gene.summary, file = gene.summary.output)
