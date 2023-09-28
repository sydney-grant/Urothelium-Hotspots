




nmibc <- read.csv("E:\\Bladder Hotspot Paper\\NMIBC_Dataset.csv")
mibc <- read.csv("E:\\Bladder Hotspot Paper\\BLCA_Dataset.csv")

impact_regions <- read.csv("E:\\Bladder Hotspot Paper\\mskimpact_regions_new.csv")

colnames(impact_regions)[[1]] <- "chr"
chr_list <- unique(impact_regions$chr)


nmibc_new <- data.frame()
mibc_new <- data.frame()

for (c in chr_list){
  pos_list <- c()
  regions_sub <- subset(impact_regions, chr == c)
  for (i in 1:nrow(regions_sub)){
    pos_list <- c(pos_list, regions_sub$start[[i]]:regions_sub$stop[[i]])
  }
  nmibc_sub <- subset(nmibc, Chromosome == c & Start_Position %in% pos_list)
  nmibc_new <- rbind(nmibc_new, nmibc_sub)
  mibc_sub <- subset(mibc, Chromosome == c & Start_Position %in% pos_list)
  mibc_new <- rbind(mibc_new, mibc_sub)
}


# identify top 10% of mutated regions panel length

mc <- nmibc_new[,c(1,5,6)]
colnames(mc) <- c("gene", "chr", "pos")

chr_list <- unique(mc$chr)

all_pos <- c()
for (c in chr_list){
  mc_sub <- subset(mc, chr == c)
  all_vec <- c()
  for (p in mc_sub$pos){
    vec <- (p - 50):(p + 50)
    all_vec <- c(all_vec, vec)
  }
  all_vec <- unique(all_vec)
  all_pos <- c(all_pos, all_vec)
}


nmibc_top_10 <- length(all_pos) * 0.1


### NMIBC 10552


mc <- mibc_new[,c(1,5,6)]
colnames(mc) <- c("gene", "chr", "pos")

chr_list <- unique(mc$chr)

all_pos <- c()
for (c in chr_list){
  mc_sub <- subset(mc, chr == c)
  all_vec <- c()
  for (p in mc_sub$pos){
    vec <- (p - 50):(p + 50)
    all_vec <- c(all_vec, vec)
  }
  all_vec <- unique(all_vec)
  all_pos <- c(all_pos, all_vec)
}


mibc_top_10 <- length(all_pos) * 0.1

### MIBC 36917


### average number of mutations per sample

mibc_samples <- unique(mibc_new$Tumor_Sample_Barcode)

nmibc_samples <- unique(nmibc_new$Tumor_Sample_Barcode)

mibc_count_list <- c()
for (m in mibc_samples){
  mibc_count_list <- c(mibc_count_list, nrow(subset(mibc_new, Tumor_Sample_Barcode == m)))
}
nmibc_count_list <- c()
for (n in nmibc_samples){
  nmibc_count_list <- c(nmibc_count_list, nrow(subset(nmibc_new, Tumor_Sample_Barcode == n)))
}

mibc_avg <- mean(mibc_count_list)

nmibc_avg <- mean(nmibc_count_list)


wilcox.test(mibc_count_list, nmibc_count_list, paired = FALSE)







