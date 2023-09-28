library(ggplot2)
library(ggrepel)

nu <- read.csv("F:\\Bladder Hotspot Paper\\Lawson_Bladder_WES_Dataset.csv")
nu_si <- read.csv("F:\\Bladder Hotspot Paper\\Lawson_SI.csv")


nu <- subset(nu, histological_feature == "Urothelium")

id_list <- c()
for (i in 1:nrow(nu)){
  id <- substr(nu[i,1], 1, 7)
  id_list <- c(id_list, id)
}


nu <- nu[,c(1,4,5)]
nu$patient <- unlist(id_list)
colnames(nu) <- c("samples", "chr", "pos", "patient")

colnames(nu_si)[[1]] <- "patient"

nu <- merge(nu, nu_si, by = "patient")

nu_highrisk <- subset(nu, Status == "Cancer Patient")
nu_lowrisk <- subset(nu, Status == "Healthy Donor")

val <- read.csv("F:\\Bladder Hotspot Paper\\Li_NU_SNV.csv")
val_si <- read.csv("F:\\Bladder Hotspot Paper\\Li_Si.csv")

colnames(val)[[1]] <- "Tumor_Sample_Barcode"

chr_list <- c()
for (i in 1:nrow(val)){
  chr <- substr(val$Chromosome[[i]], 4, nchar(val$Chromosome[[i]]))
  chr_list <- c(chr_list, chr)
}
val$Chromosome <- unlist(chr_list)
colnames(val)[[3]] <- "Start_Position"

pt_list <- c()
for (i in 1:nrow(val)){
  p <- val$Tumor_Sample_Barcode[[i]]
  p_split <- strsplit(p, "U")
  pt <- unlist(p_split)[[1]]
  pt_list <- c(pt_list, pt)
}

colnames(val_si)[[1]] <- "Patient"
val$Patient <- unlist(pt_list)

val <- merge(val, val_si, by = "Patient")

val <- val[,c(1:4,10:12)]
colnames(val) <- c("patient", "samples", "chr", "pos", "gender", "age", "smoking")

smoking_list <- c()
for (i in 1:length(val$smoking)){
  if (val$smoking[[i]] == "No"){
    smoking_list <- c(smoking_list, "Non-smoker")
  }
  if (val$smoking[[i]] == "Yes"){
    smoking_list <- c(smoking_list, "Smoker")
  }

}
val$smoking <- unlist(smoking_list)

nu_highrisk <- nu_highrisk[,1:7]
colnames(nu_highrisk) <- c("patient", "samples", "chr", "pos", "gender", "age", "smoking")

nu_lowrisk <- nu_lowrisk[,1:7]
colnames(nu_lowrisk) <- c("patient", "samples", "chr", "pos", "gender", "age", "smoking")



nu_highrisk <- rbind(nu_highrisk, val)

nu_all <- rbind(nu_highrisk, nu_lowrisk)

hr_samples <- unique(nu_highrisk$patient)
lr_samples <- unique(nu_lowrisk$patient)

#####################################################

bladder1 <- read.csv("F:\\Bladder Hotspot Paper\\Final Analysis\\MIBC_Bins_Final.csv")

nmi_bladder <- read.csv("F:\\Bladder Hotspot Paper\\Final Analysis\\NMIBC_Bins_Final.csv")

####################################################

dat <- read.csv("F:\\Bladder Hotspot Paper\\BLCA_Dataset.csv")
nmi <- read.csv("F:\\Bladder Hotspot Paper\\NMIBC_Dataset.csv")


##### add genes to hotspots
gene_list <- c()
for (i in 1:nrow(bladder1)){
  data.sub1 <- subset(dat, Chromosome == bladder1$chromosome[[i]])
  data.sub2 <- subset(data.sub1, Start_Position %in% bladder1$lowerbound[[i]]:bladder1$upperbound[[i]])
  gene_list <- c(gene_list, unique(data.sub2$ï..Hugo_Symbol)[[1]])
}
bladder1$gene <- unlist(gene_list)

gene_list <- c()
for (i in 1:nrow(nmi_bladder)){
  data.sub1 <- subset(nmi, Chromosome == nmi_bladder$chromosome[[i]])
  data.sub2 <- subset(data.sub1, Start_Position %in% nmi_bladder$lowerbound[[i]]:nmi_bladder$upperbound[[i]])
  gene_list <- c(gene_list, unique(data.sub2$Hugo_Symbol)[[1]])
}
nmi_bladder$gene <- unlist(gene_list)

bc_all <- rbind(dat[,c(5,6,17)], nmi[,c(5,6,17)])

## correlation all healthy vs all bc

nu_hotspots_all <- c()
blca_hotspots_all <- c()

bladder <- rbind(bladder1, nmi_bladder)


for (i in 1:nrow(bladder)){
  data1_sub <- subset(nu_all, chr == bladder$chromosome[[i]] &
                        pos %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count1 <- nrow(data1_sub) / length(unique(nu_all$samples))
  nu_hotspots_all <- c(nu_hotspots_all, count1)
  data2_sub <- subset(bc_all, Chromosome == bladder$chromosome[[i]] &
                        Start_Position %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count2 <- nrow(data2_sub) / length(unique(bc_all$Tumor_Sample_Barcode))
  blca_hotspots_all <- c(blca_hotspots_all, count2)
}

cor.test(nu_hotspots_all, blca_hotspots_all, method = "pearson")


## correlation nmi vs mi in mi hotspots

nmi_mi <- c()
mi_mi <- c()

bladder <- bladder1


for (i in 1:nrow(bladder)){
  data1_sub <- subset(dat, Chromosome == bladder$chromosome[[i]] &
                        Start_Position %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count2 <- nrow(data2_sub) / length(unique(dat$Tumor_Sample_Barcode))
  mi_mi <- c(mi_mi, count2)

  data2_sub <- subset(nmi, Chromosome == bladder$chromosome[[i]] &
                        Start_Position %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count2 <- nrow(data2_sub) / length(unique(nmi$Tumor_Sample_Barcode))
  nmi_mi <- c(nmi_mi, count2)
}

cor.test(nmi_mi, mi_mi, method = "pearson")


## correlation nmi vs mi in nmi hotspots

nmi_nmi <- c()
mi_nmi <- c()

bladder <- nmi_bladder


for (i in 1:nrow(bladder)){
  data1_sub <- subset(nmi, Chromosome == bladder$chromosome[[i]] &
                        Start_Position %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count1 <- nrow(data1_sub) / length(unique(nmi$Tumor_Sample_Barcode))
  nmi_nmi <- c(nmi_nmi, count1)
  data2_sub <- subset(dat, Chromosome == bladder$chromosome[[i]] &
                        Start_Position %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count2 <- nrow(data2_sub) / length(unique(dat$Tumor_Sample_Barcode))
  mi_nmi <- c(mi_nmi, count2)
}

cor.test(nmi_nmi, mi_nmi, method = "pearson")


## correlation hr lr vs mibc

hr_mi <- c()
lr_mi <- c()

bladder <- bladder1


for (i in 1:nrow(bladder)){
  data1_sub <- subset(nu_highrisk, chr == bladder$chromosome[[i]] &
                        pos %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count1 <- nrow(data1_sub) / length(unique(nu_highrisk$samples))
  hr_mi <- c(hr_mi, count1)
  data2_sub <- subset(nu_lowrisk, chr == bladder$chromosome[[i]] &
                        pos %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count2 <- nrow(data2_sub) / length(unique(nu_lowrisk$samples))
  lr_mi <- c(lr_mi, count2)
}

cor.test(hr_mi, mi_mi, method = "pearson")

cor.test(lr_mi, mi_mi, method = "pearson")


## correlation hr lr vs nmibc

hr_nmi <- c()
lr_nmi <- c()

bladder <- nmi_bladder


for (i in 1:nrow(bladder)){
  data1_sub <- subset(nu_highrisk, chr == bladder$chromosome[[i]] &
                        pos %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count1 <- nrow(data1_sub) / length(unique(nu_highrisk$samples))
  hr_nmi <- c(hr_nmi, count1)
  data2_sub <- subset(nu_lowrisk, chr == bladder$chromosome[[i]] &
                        pos %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
  count2 <- nrow(data2_sub) / length(unique(nu_lowrisk$samples))
  lr_nmi <- c(lr_nmi, count2)
}

cor.test(hr_nmi, nmi_nmi, method = "pearson")

cor.test(lr_nmi, nmi_nmi, method = "pearson")



all_hotspots <- data.frame("NU" = unlist(nu_hotspots_all), "BC" = unlist(blca_hotspots_all),
                           "gene" = c(bladder1$gene, nmi_bladder$gene))

mi_hotspots <- data.frame("MIBC" = unlist(mi_mi), "HR" = unlist(hr_mi),
                          "LR" = unlist(lr_mi), "NMIBC" = unlist(nmi_mi), "gene" = bladder1$gene)

nmi_hotspots <- data.frame("MIBC" = unlist(mi_nmi), "HR" = unlist(hr_nmi),
                          "LR" = unlist(lr_nmi), "NMIBC" = unlist(nmi_nmi), "gene" = nmi_bladder$gene)

plot1.ann = expression(paste("R = 0.385, ", italic(p), " < 2.2 x 10"^"\u221216", sep = ""))
plot1 <- ggplot(all_hotspots, aes(x = NU, y = BC)) +
  geom_jitter(data = subset(all_hotspots, NU <= 0.012 & BC <= 0.1), size = 2, shape = 19, color = c("#E7B800")) +
  geom_point(data = subset(all_hotspots, NU > 0.012 | BC > 0.1), size = 2, shape = 19, color = c("#E7B800")) +
  geom_text_repel(data = subset(all_hotspots, NU > 0.01 | BC > 0.1), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("NU \n (Mutations per Region)") +
  ylab("Bladder Cancer \n (Mutations per Region)") +
  ggtitle("BC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  xlim(0,0.025) +
  ylim(0,0.28) +
  geom_text(x = 0.0125, y = 0.27, label = plot1.ann, size = 4)

plot1

plotb1.ann = expression(paste("R = 0.341, ", italic(p), " = 1.59 x 10"^"\u221211", sep = ""))
plotb1 <- ggplot(mi_hotspots, aes(x = NMIBC, y = MIBC)) +
  geom_jitter(data = subset(mi_hotspots, NMIBC <= 0.1 & MIBC <= 0.05), size = 2, shape = 19, color = c("#FC4E07")) +
  geom_point(data = subset(mi_hotspots, NMIBC > 0.1 | MIBC > 0.05), size = 2, shape = 19, color = c("#FC4E07")) +
  geom_text_repel(data = subset(mi_hotspots, NMIBC > 0.1 | MIBC > 0.05), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("NMIBC \n (Mutations per Region)") +
  ylab("MIBC \n (Mutations per Region)") +
  ggtitle("MIBC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  #xlim(0,0.42) +
  #ylim(0,0.28) +
  geom_text(x = 0.21, y = 0.08, label = plotb1.ann, size = 4)

plotb1

plotb2.ann = expression(paste("R = 0.627, ", italic(p), " = 8.07 x 10"^"\u221213", sep = ""))
plotb2 <- ggplot(nmi_hotspots, aes(x = NMIBC, y = MIBC)) +
  geom_jitter(data = subset(nmi_hotspots, NMIBC <= 0.1 & MIBC <= 0.07), size = 2, shape = 19, color = c("#00AFBB")) +
  geom_point(data = subset(nmi_hotspots, NMIBC > 0.1 | MIBC > 0.07), size = 2, shape = 19, color = c("#00AFBB")) +
  geom_text_repel(data = subset(nmi_hotspots, NMIBC > 0.1 | MIBC > 0.07), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("NMIBC \n (Mutations per Region)") +
  ylab("MIBC \n (Mutations per Region)") +
  ggtitle("NMIBC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  xlim(0,0.82) +
  ylim(0,0.25) +
  geom_text(x = 0.41, y = 0.24, label = plotb2.ann, size = 4)

plotb2

plothrmi.ann = expression(paste("R = 0.425, ", italic(p), " < 2.2 x 10"^"\u221216", sep = ""))
plot.hr.mi <- ggplot(mi_hotspots, aes(x = HR, y = MIBC)) +
  geom_jitter(data = subset(mi_hotspots, HR <= 0.01 & MIBC <= 0.025), size = 2, shape = 19, color = c("#FC4E07")) +
  geom_point(data = subset(mi_hotspots, HR > 0.01 | MIBC > 0.025), size = 2, shape = 19, color = c("#FC4E07")) +
  geom_text_repel(data = subset(mi_hotspots, HR > 0.012 | MIBC > 0.025), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("High-Risk NU \n (Mutations per Region)") +
  ylab("MIBC \n (Mutations per Region)") +
  ggtitle("MIBC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  #xlim(0,0.08) +
  #ylim(0,0.25) +
  geom_text(x = 0.01, y = 0.08, label = plothrmi.ann, size = 4)

plot.hr.mi

plotlrmi.ann = expression(paste("R = 0.140, ", italic(p), " = 0.00698", sep = ""))
plot.lr.mi <- ggplot(mi_hotspots, aes(x = LR, y = MIBC)) +
  geom_jitter(data = subset(mi_hotspots, LR <= 0.01 & MIBC <= 0.05), size = 2, shape = 19, color = c("#FC4E07")) +
  geom_point(data = subset(mi_hotspots, LR > 0.01 | MIBC > 0.05), size = 2, shape = 19, color = c("#FC4E07")) +
  geom_text_repel(data = subset(mi_hotspots, LR > 0.01 | MIBC > 0.05), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("Low-Risk NU \n (Mutations per Region)") +
  ylab("MIBC \n (Mutations per Region)") +
  ggtitle("MIBC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  #xlim(0,0.025) +
  #ylim(0,0.25) +
  geom_text(x = 0.0125, y = 0.095, label = plotlrmi.ann, size = 4)

plot.lr.mi

plothrnmi.ann = expression(paste("R = 0.320, ", italic(p), " = 0.0009", sep = ""))
plot.hr.nmi <- ggplot(nmi_hotspots, aes(x = HR, y = NMIBC)) +
  geom_jitter(data = subset(nmi_hotspots, HR <= 0.01 & NMIBC <= 0.2), size = 2, shape = 19, color = c("#00AFBB")) +
  geom_point(data = subset(nmi_hotspots, HR > 0.01 | NMIBC > 0.2), size = 2, shape = 19, color = c("#00AFBB")) +
  geom_text_repel(data = subset(nmi_hotspots, HR > 0.01 | NMIBC > 0.2), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("High-Risk NU \n (Mutations per Region)") +
  ylab("NMIBC \n (Mutations per Region)") +
  ggtitle("NMIBC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  xlim(0,0.025) +
  ylim(0,0.83) +
  geom_text(x = 0.0175, y = 0.8, label = plothrnmi.ann, size = 4)

plot.hr.nmi

require(ggrepel)

plotlrnmi.ann = expression(paste("R = 0.082, ", italic(p), " = 0.404", sep = ""))
plot.lr.nmi <- ggplot(nmi_hotspots, aes(x = LR, y = NMIBC)) +
  geom_jitter(data = subset(nmi_hotspots, LR <= 0.01 & NMIBC <= 0.2), size = 2, shape = 19, color = c("#00AFBB")) +
  geom_point(data = subset(nmi_hotspots, LR > 0.01 | NMIBC > 0.2), size = 2, shape = 19, color = c("#00AFBB")) +
  geom_text_repel(data = subset(nmi_hotspots, LR > 0.01 | NMIBC > 0.2), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("Low-Risk NU \n (Mutations per Region)") +
  ylab("NMIBC \n (Mutations per Region)") +
  ggtitle("NMIBC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  xlim(-0.001,0.02) +
  ylim(-0.001,0.83) +
  geom_text(x = 0.01, y = 0.82, label = plotlrnmi.ann, size = 4)

plot.lr.nmi

cor.test(lr_mi, hr_mi, method = "pearson")

plotriskcomp.mi.ann = expression(paste("R = 0.0314, ", italic(p), " = 0.547", sep = ""))
plot.riskcomp.mi <- ggplot(mi_hotspots, aes(x = LR, y = HR)) +
  geom_jitter(data = subset(mi_hotspots, LR <= 0.01 & HR <= 0.01), size = 2, shape = 19, color = c("#FC4E07")) +
  geom_point(data = subset(mi_hotspots, LR > 0.01 | HR > 0.01), size = 2, shape = 19, color = c("#FC4E07")) +
  geom_text_repel(data = subset(mi_hotspots, LR > 0.01 | HR > 0.012), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("Low-Risk NU \n (Mutations per Region)") +
  ylab("High-Risk NU \n (Mutations per Region)") +
  ggtitle("MIBC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  #xlim(-0.001,0.025) +
  #ylim(-0.001,0.083) +
  geom_text(x = 0.0175, y = 0.015, label = plotriskcomp.mi.ann, size = 4)

plot.riskcomp.mi

cor.test(lr_nmi, hr_nmi, method = "pearson")

plotriskcomp.nmi.ann = expression(paste("R = \u22120.020, ", italic(p), " = 0.837", sep = ""))
plot.riskcomp.nmi <- ggplot(nmi_hotspots, aes(x = LR, y = HR)) +
  geom_jitter(data = subset(nmi_hotspots, LR <= 0.005 & HR <= 0.01), size = 2, shape = 19, color = c("#00AFBB")) +
  geom_point(data = subset(nmi_hotspots, LR > 0.005 | HR > 0.01), size = 2, shape = 19, color = c("#00AFBB")) +
  geom_text_repel(data = subset(nmi_hotspots, LR > 0.01 | HR > 0.012), aes(label = gene),
                  size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), force = 20) +
  xlab("Low-Risk NU \n (Mutations per Region)") +
  ylab("High-Risk NU \n (Mutations per Region)") +
  ggtitle("NMIBC10") +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  theme_classic() +
  xlim(-0.001,0.02) +
  ylim(-0.001,0.02) +
  geom_text(x = 0.012, y = 0.018, label = plotriskcomp.nmi.ann, size = 4)

plot.riskcomp.nmi

count_df1 <- read.csv("F:\\Bladder Hotspot Paper\\LRvsHR_MIBC_Bins.csv")

count_df21Type <- factor(count_df1$Type, levels = c("Low-Risk", "High-Risk"))

wc_test <- wilcox.test(hr_count, lr_count, paired = FALSE)

fig2f.1 <- ggplot(count_df1, aes(y = Count, x = Type)) +
  geom_violin(aes(fill = Type), trim = FALSE) +
  ylab("Mutations per Sample") +
  scale_fill_manual(values = c("coral4", "coral3")) +
  xlab(NULL) +
  ggtitle("MBIC Hotspots") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 10),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12),
        legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) +
  geom_signif(comparisons=list(c("High-Risk", "Low-Risk")), annotations="***",
              y_position = 14, tip_length = 0, vjust=0.5)

count_df2 <- read.csv("F:\\Bladder Hotspot Paper\\LRvsHR_NMIBC_Bins.csv")

count_df2$Type <- factor(count_df2$Type, levels = c("Low-Risk", "High-Risk"))

wc_test <- wilcox.test(hr_count, lr_count, paired = FALSE)

fig2f.2 <- ggplot(count_df2, aes(y = Count, x = Type)) +
  geom_violin(aes(fill = Type), trim = FALSE) +
  ylab("Mutations per Sample") +
  scale_fill_manual(values = c("plum4", "plum3")) +
  xlab(NULL) +
  theme_classic() +
  ggtitle("NMIBC Hotspots") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 10),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12),
        legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14))


#B <- ggarrange(plotb1, plotb2, ncol = 2)
A.B <- ggarrange(plotb1, plotb2,plot1,  ncol = 3, labels = c("A", " ", "B"))
C <- ggarrange(plot.hr.mi, plot.lr.mi, ncol = 2, labels = c("C"))
D <- ggarrange(plot.hr.nmi, plot.lr.nmi, ncol = 2, labels = c("D"))
E <- ggarrange(plot.riskcomp.mi, plot.riskcomp.nmi, ncol = 2, labels = c("E", " "))
#F <- ggarrange(fig2f.1, fig2f.2, nrow = 2, labels = c("F"))

FIG2 <- ggarrange(A.B, C, D, E, nrow = 4)

FIG2



