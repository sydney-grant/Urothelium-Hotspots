
library(ggplot2)
library(ggpubr)
library(ggsignif)

bladder <- read.csv("F:\\Bladder Hotspot Paper\\Final Analysis\\MIBC_Bins_Final.csv")

nmi_bladder <- read.csv("F:\\Bladder Hotspot Paper\\Final Analysis\\NMIBC_Bins_Final.csv")

bladder$percent <- (bladder$count / 409) * 100


vec <- seq(from = 0, to = 37000, by = 500)

group_list <- c()
g_list <- c()
for (v in 1:(length(vec)-1)){
  data_sub <- subset(bladder, Cummulative_Bin_Length <= vec[[v+1]] & Cummulative_Bin_Length > vec[[v]])
  lab <- paste(vec[[v]], "-",vec[[v+1]], "bp" )
  group_list <- c(group_list, rep(lab, nrow(data_sub)))
  g_list <- c(g_list, lab)
}

bladder$label <- unlist(group_list)
bladder$label <- factor(bladder$label, levels = g_list)
bladder$number <- 1:369

plot <- ggplot(bladder, aes(y = percent, x = number)) +
  geom_point(data = subset(bladder, number <= 13),
             color = "#FC4E07", fill = "#FC4E07", alpha = 1, shape = 21, size = 2) +
  geom_point(data = subset(bladder, number > 13),
             color = "gray", fill = "gray", alpha = 1, shape = 21, size = 2) +
  geom_vline(xintercept = 14, linetype = "longdash") +
  theme_classic() +
  ylab("Percentage of Samples \n with Hotspot Mutation") +
  xlab("Ranked Region Number") +
  ggtitle("Muscle Invasive Bladder Cancer") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title = element_text(size = 4),
        axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 12)) +
  ylim(0,20) +
  theme(legend.position = "none")

plot

nmi_bladder$number <- 1:105

nmi_bladder$percent <- (nmi_bladder$count / 103) * 100

plot2 <- ggplot(nmi_bladder, aes(y = percent, x = number)) +
  geom_point(data = subset(bladder, number <= 4),
             color = "#00AFBB", fill = "#00AFBB", alpha = 1, shape = 21, size = 2) +
  geom_point(data = subset(bladder, number > 4),
             color = "gray", fill = "gray", alpha = 1, shape = 21, size = 2) +
  geom_vline(xintercept = 5, linetype = "longdash") +
  theme_classic() +
  ggtitle("Non-Muscle Invasive Bladder Cancer") +
  ylab("Percentage of Samples \n with Hotspot Mutation") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title = element_text(size = 4),
        axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 12)) +
  xlab("Ranked Region Number") +
  ylim(0,20) +
  theme(legend.position = "none")

plot2

##############################################################
##############################################################
bladder <- read.csv("F:\\Bladder Hotspot Paper\\Final Analysis\\mi_significant_hotspots.csv")

nmi_bladder <- read.csv("F:\\Bladder Hotspot Paper\\Final Analysis\\nmi_significant_hotspots.csv")

mi_freq <- as.data.frame(table(bladder$gene))
mi_freq <- mi_freq[order(-mi_freq$Freq),]
mi_freq$Var1 <- factor(mi_freq$Var1, levels = c(mi_freq$Var1))

nmi_freq <- as.data.frame(table(nmi_bladder$gene))
nmi_freq <- nmi_freq[order(-nmi_freq$Freq),]
nmi_freq$Var1 <- factor(nmi_freq$Var1, levels = c(nmi_freq$Var1))


bladder_plot <- ggplot(mi_freq, aes(x=Var1, y=Freq)) +
  geom_segment( aes(x=Var1, xend=Var1, y=0, yend=Freq), color="#FC4E07", size = 1) +
  geom_point( color="#FC4E07", size=4) +
  ylab("Number of Hotspots") +
  xlab("Gene") +
  ggtitle("Muscle Invasive Bladder Cancer") +
  theme_classic()+
  ylim(0,6)+
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12, angle = 90),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title = element_text(hjust = 0.5, size = 12),
        plot.title= element_text(hjust = 0.5, size = 12))
bladder_plot

nmi_bladder_plot <- ggplot(nmi_freq, aes(x=Var1, y=Freq)) +
  geom_segment( aes(x=Var1, xend=Var1, y=0, yend=Freq), color="#00AFBB", size = 1) +
  geom_point( color="#00AFBB", size=4) +
  ylab("Number of Hotspots") +
  xlab("Gene") +
  theme_classic()+
  ylim(0,6)+
  ggtitle("Non-Muscle Invasive Bladder Cancer") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12, angle = 90),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title = element_text(hjust = 0.5, size = 12),
        plot.title= element_text(hjust = 0.5, size = 12))
nmi_bladder_plot

###########################################################
###########################################################

b_data <- read.csv("F:\\Bladder Hotspot Paper\\Li_UCC_SNV.csv")
colnames(b_data)[[1]] <- "Tumor_Sample_Barcode"

chr_list <- c()
for (i in 1:nrow(b_data)){
  chr <- substr(b_data$Chromosome[[i]], 4, nchar(b_data$Chromosome[[i]]))
  chr_list <- c(chr_list, chr)
}
b_data$Chromosome <- unlist(chr_list)
colnames(b_data)[[3]] <- "Start_Position"

samples <- unique(b_data$Tumor_Sample_Barcode)

count_list <- c()
for (s in samples){
  data.sub <- subset(b_data, Tumor_Sample_Barcode == s)
  count <- 0
  for (i in 1:nrow(bladder)){
    data.sub2 <- subset(data.sub, Chromosome == bladder$chromosome[[i]])
    data.sub3 <- subset(data.sub2, Start_Position %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
    count <- count + nrow(data.sub3)
  }
  count_list <- c(count_list, count)
}

count_df_high <- data.frame("Count" = unlist(count_list), "Sample" = unlist(samples))

bladder_low <- read.csv("F:\\Bladder Hotspot Paper\\Low_BLCA_Bins.csv")
bladder_low <- bladder_low[1:13,]

count_list <- c()
for (s in samples){
  data.sub <- subset(b_data, Tumor_Sample_Barcode == s)
  count <- 0
  for (i in 1:nrow(bladder_low)){
    data.sub2 <- subset(data.sub, Chromosome == bladder_low$chromosome[[i]])
    data.sub3 <- subset(data.sub2, Start_Position %in% bladder_low$lowerbound[[i]]:bladder_low$upperbound[[i]])
    count <- count + nrow(data.sub3)
  }
  count_list <- c(count_list, count)
}

count_df_low <- data.frame("Count" = unlist(count_list), "Sample" = unlist(samples))

wtest <- wilcox.test(count_df_high$Count, count_df_low$Count)


fige_df <- data.frame("Frequency" = c(count_df_high$Count,count_df_low$Count),
                      "Panel" = c(rep("Hotspots", nrow(count_df_high)), rep("Non-Hotspots", nrow(count_df_low))))



fig1e <- ggplot(fige_df, aes(y = Frequency, x = Panel)) +
  geom_violin(aes(fill = Panel), trim = FALSE) +
  ylab("Mutations per Sample") +
  scale_fill_manual(values = c("#bf3b06", "#FC4E07")) +
  theme_classic()+
  xlab(NULL) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        legend.position = "none", plot.title= element_text(size = 11)) +
  geom_signif(comparisons=list(c("Hotspots", "Non-Hotspots")), annotations="***",
              y_position = 4, tip_length = 0, vjust=0.5)
fig1e

dat <- read.csv("F:\\Bladder Hotspot Paper\\BLCA_Dataset.csv")
nmi <- read.csv("F:\\Bladder Hotspot Paper\\NMIBC_Dataset.csv")

####### FGFR3

count_list.nmi <- c()

for (i in unique(nmi$Tumor_Sample_Barcode)){
  data.sub1 <- subset(nmi, Tumor_Sample_Barcode == i)
  data.sub <- subset(data.sub1, Chromosome == 4 & Start_Position %in% 1803519:1803618)
  count_list.nmi <- c(count_list.nmi, nrow(data.sub))
}

count_list.mi <- c()

for (i in unique(dat$Tumor_Sample_Barcode)){
  data.sub1 <- subset(dat, Tumor_Sample_Barcode == i)
  data.sub <- subset(data.sub1, Chromosome == 4 & Start_Position %in% 1803519:1803618)
  count_list.mi <- c(count_list.mi, nrow(data.sub))
}

wilcox.test(count_list.mi, count_list.nmi, paired = FALSE)

fc.fgfr3 <- mean(count_list.nmi) / mean(count_list.mi)
fgfr3.pval <- wilcox.test(count_list.mi, count_list.nmi, paired = FALSE)$p.value

plotfgfr3.ann = expression("FGFR3 chr: 4 bp: 1803519\u22121803618")
fgfr3.df <- data.frame("Count" = c(count_list.mi, count_list.nmi), "Type" = c(rep("MIBC", length(count_list.mi)), rep("NMIBC", length(count_list.nmi))))
fgfr3.df$Type <- factor(fgfr3.df$Type, levels = c("NMIBC", "MIBC"))
fgfr3.plot <- ggplot(fgfr3.df, aes(y = Count, x = Type)) +
  geom_violin(aes(fill = Type), trim = FALSE) +
  ylab("Mutations per Sample") +
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_classic()+
  xlab(NULL) +
  ggtitle(plotfgfr3.ann) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        legend.position = "none", plot.title= element_text(size = 12)) +
  geom_signif(comparisons=list(c("NMIBC", "MIBC")), annotations="***",
              y_position = 3, tip_length = 0, vjust=0.5)
fgfr3.plot


####### PIK3CA (NMIBC) ********** p < 0.05 higher in NMIBC

count_list.nmi <- c()

for (i in unique(nmi$Tumor_Sample_Barcode)){
  data.sub1 <- subset(nmi, Tumor_Sample_Barcode == i)
  data.sub <- subset(data.sub1, Chromosome == 3 & Start_Position %in% 178936038:178936137)
  count_list.nmi <- c(count_list.nmi, nrow(data.sub))
}

count_list.mi <- c()

for (i in unique(dat$Tumor_Sample_Barcode)){
  data.sub1 <- subset(dat, Tumor_Sample_Barcode == i)
  data.sub <- subset(data.sub1, Chromosome == 3 & Start_Position %in% 178936038:178936137)
  count_list.mi <- c(count_list.mi, nrow(data.sub))
}

wilcox.test(count_list.mi, count_list.nmi, paired = FALSE)

fc.pik3ca <- mean(count_list.nmi) / mean(count_list.mi)

pik3ca.pval <- wilcox.test(count_list.mi, count_list.nmi, paired = FALSE)$p.value

plotpik3ca.ann = expression("PIK3CA chr: 3 bp: 178936038\u2212178936137")
pik3ca.df <- data.frame("Count" = c(count_list.mi, count_list.nmi), "Type" = c(rep("MIBC", length(count_list.mi)), rep("NMIBC", length(count_list.nmi))))
pik3ca.df$Type <- factor(pik3ca.df$Type, levels = c("NMIBC", "MIBC"))
pik3ca.plot <- ggplot(pik3ca.df, aes(y = Count, x = Type)) +
  geom_violin(aes(fill = Type), trim = FALSE) +
  ylab("Mutations per Sample") +
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_classic()+
  xlab(NULL) +
  ggtitle(plotpik3ca.ann) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        legend.position = "none", plot.title= element_text(size = 12)) +
  geom_signif(comparisons=list(c("MIBC", "NMIBC")), annotations="*",
              y_position = 3, tip_length = 0, vjust=0.5)
pik3ca.plot

####### TP53(MIBC) ********** p < 0.05 higher in MIBC

count_list.nmi <- c()

for (i in unique(nmi$Tumor_Sample_Barcode)){
  data.sub1 <- subset(nmi, Tumor_Sample_Barcode == i)
  data.sub <- subset(data.sub1, Chromosome == 17 & Start_Position %in% 7577494:7577593)
  count_list.nmi <- c(count_list.nmi, nrow(data.sub))
}

count_list.mi <- c()

for (i in unique(dat$Tumor_Sample_Barcode)){
  data.sub1 <- subset(dat, Tumor_Sample_Barcode == i)
  data.sub <- subset(data.sub1, Chromosome == 17 & Start_Position %in% 7577494:7577593)
  count_list.mi <- c(count_list.mi, nrow(data.sub))
}

tp53.pval <- wilcox.test(count_list.mi, count_list.nmi, paired = FALSE)$p.value

fc.tp53 <- mean(count_list.mi) / mean(count_list.nmi)

plottp53.ann = expression("TP53 chr: 17 bp: 7577494\u22127577593")
tp53.df <- data.frame("Count" = c(count_list.mi, count_list.nmi), "Type" = c(rep("MIBC", length(count_list.mi)), rep("NMIBC", length(count_list.nmi))))
tp53.df$Type <- factor(tp53.df$Type, levels = c("NMIBC", "MIBC"))
tp53.plot <- ggplot(tp53.df, aes(y = Count, x = Type)) +
  geom_violin(aes(fill = Type), trim = FALSE) +
  ylab("Mutations per Sample") +
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_classic()+
  xlab(NULL) +
  ggtitle(plottp53.ann) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        legend.position = "none", plot.title= element_text(size = 12)) +
  geom_signif(comparisons=list(c("MIBC", "NMIBC")), annotations="*",
              y_position = 3, tip_length = 0, vjust=0.5)
tp53.plot


A <- ggarrange(plot2, plot, ncol = 2, labels = c("A"))
B <- ggarrange(nmi_bladder_plot, bladder_plot, ncol = 2, labels = c("B"))
CD <- ggarrange(fig1e, tp53.plot, ncol = 2, labels = c("C", "D"))
E <- ggarrange(fgfr3.plot, pik3ca.plot, ncol = 2, labels = c("E", "F"))

FIG1 <- ggarrange(A, B, CD, E, nrow = 4)
FIG1




