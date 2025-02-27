library(tidyverse)

# ENVIRONMENT
setwd('/lisc/scratch/botany/frschmidt/mk_test/')

# FUNCTIONS
build_region_df <- function(dir, file) {
  region.df <- read.delim(paste0(dir, file), 
                          header = FALSE, sep = ":", 
                          col.names = c("Scaffold", "Range"))
  region.df <- region.df %>%
    separate(col = Range, into = c("RangeStart", "RangeEnd"), sep = "-")
  region.df <- transform(region.df, RangeStart = as.numeric(RangeStart))
  region.df <- transform(region.df, RangeEnd = as.numeric(RangeEnd))
  region.df$Length <- region.df$RangeEnd - region.df$RangeStart
  return(region.df)
}

plot_cat_bar <- function(df.list) {
  p <- df.list[[1]] %>%
    ggplot(aes(x = Scaffold), environment = environment()) + 
    geom_bar(aes(y = (..count..)/sum(..count..))) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "", y = "Frequency")
  p + ggtitle(names(df.list))
}

# VARIABLES
transdecoder.dir <- 'results/transdecoder/'
transdecoder.complete.orf.dir <- 'results/complete_orf_transdecoder/'
transdecoder.subsets <- c("vie_gene.fa.transdecoder.missing.regions", 
                          "vie_gene.fa.transdecoder.chim.regions",
                          "vie_gene.fa.transdecoder.uniq.regions")

# EXECUTION
# extract scaffold, coordinate pairs into seperate data frames
transdecoder.subsets.list <- sapply(transdecoder.subsets, 
                                    build_region_df, dir = transdecoder.dir,
                                    simplify = FALSE, USE.NAMES = TRUE)

# PLOTTING
# plot categorical histograms (barplots) of scaffolds
plot_cat_bar(transdecoder.subsets.list[transdecoder.subsets[1]])
plot_cat_bar(transdecoder.subsets.list[transdecoder.subsets[2]])
plot_cat_bar(transdecoder.subsets.list[transdecoder.subsets[3]])

# plot histograms of length of transcripts
hist(transdecoder.subsets.list$vie_gene.fa.transdecoder.uniq.regions$Length, breaks = seq(0, max(transdecoder.subsets.list$vie_gene.fa.transdecoder.uniq.regions$Length)+100, by=100))
hist(transdecoder.subsets.list$vie_gene.fa.transdecoder.chim.regions$Length, breaks = seq(0, max(transdecoder.subsets.list$vie_gene.fa.transdecoder.chim.regions$Length)+100, by=100))
hist(transdecoder.subsets.list$vie_gene.fa.transdecoder.missing.regions$Length, breaks = seq(0, max(transdecoder.subsets.list$vie_gene.fa.transdecoder.missing.regions$Length)+100, by=100))
  

transdecoder <- read.delim(paste0(transdecoder.dir, "vie_gene.fa.transdecoder.csv"),
           header = FALSE, sep = ";", 
           col.names = c("ID", "ID2", "ORFtype"))
transdecoder <- subset(transdecoder, select = -ID2)

transdecoder <- transdecoder %>%
  separate(col = ID, into = c("ID0", "ID"), sep = "=")
transdecoder <- subset(transdecoder, select = -ID0)

transdecoder$Scaffold <- transdecoder$ID

transdecoder <- transdecoder %>%
  separate(col = Scaffold, into = c("Scaffold", "ID0"), sep = "\\.")
transdecoder <- subset(transdecoder, select = -ID0)

transdecoder <- transdecoder %>%
  separate(col = Scaffold, into = c("Scaffold", "Range"), sep = ":")
transdecoder <- transdecoder %>%
  separate(col = Range, into = c("RangeStart", "RangeEnd"), sep = "-")
transdecoder <- transform(transdecoder, RangeStart = as.numeric(RangeStart))
transdecoder <- transform(transdecoder, RangeEnd = as.numeric(RangeEnd))
transdecoder$Length <- transdecoder$RangeEnd - transdecoder$RangeStart

transdecoder <- transdecoder %>%
  separate(col = ORFtype, into = c("ORFtype", "Score"), sep = ",")
transdecoder <- transdecoder %>%
  separate(col = ORFtype, into = c("ID0", "ORFtype"), sep = ":")
transdecoder <- subset(transdecoder, select = -ID0)
transdecoder <- transdecoder %>%
  separate(col = Score, into = c("ID0", "Score"), sep = "=")
transdecoder <- subset(transdecoder, select = -ID0)

transdecoder %>% 



