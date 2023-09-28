

#### measure length for hotspot panel size

mi_data <- read.csv("E:\\Bladder Hotspot Paper\\BLCA_Dataset.csv")
nmi_data <- read.csv("E:\\Bladder Hotspot Paper\\NMIBC_Dataset.csv")

mi_hotspots <- read.csv("E:\\Bladder Hotspot Paper\\Final Analysis\\MIBC_Bins_Final.csv")
nmi_hotspots <- read.csv("E:\\Bladder Hotspot Paper\\Final Analysis\\NMIBC_Bins_Final.csv")



## nmi

vec <- seq(from = 0, to = 10500, by = 100)

nmi_samples <- unique(nmi_data$Tumor_Sample_Barcode)

nmi_df <- data.frame()

for (i in vec){
  hotspot_sub <- subset(nmi_hotspots, Cummulative_Bin_Length >= i & Cummulative_Bin_Length < (i + 101))
  count_list <- c()
    for (sample in nmi_samples){
      data_sub <- subset(nmi_data, Tumor_Sample_Barcode == sample)
      count <- 0
      for (j in 1:nrow(hotspot_sub)){
        data_sub2 <- subset(data_sub, Chromosome == hotspot_sub$chromosome[[j]] &
                            Start_Position %in% hotspot_sub$lowerbound[[j]]:hotspot_sub$upperbound[[j]])
        count <- count + nrow(data_sub2)
      }
      count_list <- c(count_list, count)
    }
  row <- data.frame("Count" = unlist(count_list), "Hotspots" = paste(i, "-", (i+100), "bp"))
  nmi_df <- rbind(nmi_df, row)
}


nmi_df$Hotspots <- factor(nmi_df$Hotspots, levels = unique(nmi_df$Hotspots))

anov <- aov(Count ~ Hotspots, data = nmi_df)
summary(anov)

TukeyHSD(anov)

nmi_anov <- as.data.frame(TukeyHSD(anov)$Hotspots)


colnames(nmi_anov)[[4]] <- "pval"
nmi_anov <- subset(nmi_anov, pval < 0.1)




vec <- seq(from = 0, to = 20000, by = 100)

mi_hotspots <- mi_hotspots[1:200,]

mi_samples <- unique(mi_data$Tumor_Sample_Barcode)

mi_df <- data.frame()

for (i in vec){
  hotspot_sub <- subset(mi_hotspots, Cummulative_Bin_Length >= i & Cummulative_Bin_Length < (i + 101))
  count_list <- c()
  for (sample in mi_samples){
    data_sub <- subset(mi_data, Tumor_Sample_Barcode == sample)
    count <- 0
    for (j in 1:nrow(hotspot_sub)){
      data_sub2 <- subset(data_sub, Chromosome == hotspot_sub$chromosome[[j]] &
                            Start_Position %in% hotspot_sub$lowerbound[[j]]:hotspot_sub$upperbound[[j]])
      count <- count + nrow(data_sub2)
    }
    count_list <- c(count_list, count)
  }
  row <- data.frame("Count" = unlist(count_list), "Hotspots" = paste(i, "-", (i+100), "bp"))
  mi_df <- rbind(mi_df, row)
}


mi_df$Hotspots <- factor(mi_df$Hotspots, levels = unique(mi_df$Hotspots))


anov.mi2 <- aov(Count ~ Hotspots, data = mi_df)
summary(anov.mi2)

mi_anov <- as.data.frame(TukeyHSD(anov.mi2)$Hotspots)


colnames(mi_anov)[[4]] <- "pval"
mi_anov <- subset(mi_anov, pval < 0.1)



