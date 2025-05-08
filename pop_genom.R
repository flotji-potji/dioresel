library(tidyverse)
library(ggpubr)

##############################  LD-Decay  ######################################
my_bins <- "/lisc/scratch/botany/frschmidt/ld2/umbrosa.ld_decay_bins"

ld_bins <- read_tsv(my_bins)

ggplot(ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2)) #+ xlim(c(0, 150000))


##############################  FST  ###########################################
fst <- read_tsv("/lisc/scratch/botany/frschmidt/data/chr_ptg000001l_fst.weir.fst")

ggplot(fst, aes(POS, WEIR_AND_COCKERHAM_FST)) + geom_point()

threshold <- quantile(fst$WEIR_AND_COCKERHAM_FST, 0.975, na.rm = TRUE)

fst <- fst %>%
  mutate(outlier = ifelse(WEIR_AND_COCKERHAM_FST > threshold, "outlier", "background"))

fst %>% group_by(outlier) %>% tally()

ggplot(fst, aes(POS, WEIR_AND_COCKERHAM_FST, colour = outlier)) + geom_point()


#############################  VCF Filtering  ##################################
setwd("/lisc/scratch/frschmidt_botany/frschmidt/filtering/")
var_qual <- read_delim("chr_ptg000001l_random_sample.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
ggplot(var_qual, aes(qual)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()

var_depth <- read_delim("chr_ptg000001l_random_sample.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()
summary(var_depth$mean_depth)

var_miss <- read_delim("chr_ptg000001l_random_sample.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()
summary(var_miss$fmiss)

var_freq <- read_delim("chr_ptg000001l_random_sample.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()
summary(var_freq$maf)

ind_depth <- read_delim("chr_ptg000001l_random_sample.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
ggplot(ind_depth, aes(depth)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()

ind_miss  <- read_delim("chr_ptg000001l_random_sample.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
ggplot(ind_miss, aes(fmiss)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()

ind_het <- read_delim("chr_ptg000001l_random_sample.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()


##########################  PCA  ###############################################

setwd("/lisc/scratch/botany/frschmidt/pca/")
pca <- read_table2("ptg000001l.eigenvec", col_names = FALSE)
eigenval <- scan("ptg000001l.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))
spp[grep("flavo", pca$ind)] <- "flavocarpa"
spp[grep("fasti", pca$ind)] <- "flavocarpa"
spp[grep("umbro", pca$ind)] <- "umbrosa"

other <- rep("rest", length(pca$ind))
ultramafic <- c("Pig5223", "DB1117", "BT128", "BT127", "BT126", "BT130")
fastidiosa <- c("sa1042", "sa1045")
not_umbrosa <- c("NOT")
other[grep(paste(ultramafic, collapse = "|"), pca$ind)] <- "ultramafic"
other[grep(paste(fastidiosa, collapse = "|"), pca$ind)] <- "fastidiosa"
other[grep(paste(not_umbrosa, collapse = "|"), pca$ind)] <- "not_umbrosa"

pca <- as.tibble(data.frame(pca, spp, other))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
ggplot(pve, aes(PC, pve)) + 
  geom_bar(stat = "identity") + 
  ylab("Percentage variance explained") + 
  theme_light()

cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = other)) + 
      geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

c <- ggplot(pca, aes(PC3, PC4, col = spp, shape = other)) + 
  geom_point(size = 3)
c <- c + scale_colour_manual(values = c("red", "blue"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + 
  ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)"))

ggarrange(b, c, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")

View(pca)










