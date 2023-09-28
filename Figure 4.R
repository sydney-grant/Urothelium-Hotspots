


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

bladder <- read.csv("F:\\Bladder Hotspot Paper\\Final Analysis\\mi_significant_hotspots.csv")

nmi_bladder <- read.csv("F:\\Bladder Hotspot Paper\\Final Analysis\\nmi_significant_hotspots.csv")

####################################################


mi_hr_count <- c()
nmi_hr_count <- c()
both_hr_count <- c()
for (s in hr_samples){
  data.sub <- subset(nu_highrisk, patient == s)
  mi_count <- 0
  nmi_count <- 0
  both_count <- 0
  for (i in 1:nrow(bladder)){
    data.sub2 <- subset(data.sub, chr == bladder$chromosome[[i]])
    data.sub3 <- subset(data.sub, pos %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
    mi_count <- mi_count + nrow(data.sub3)
    both_count <- both_count + nrow(data.sub3)
  }

  for (i in 1:nrow(nmi_bladder)){
    data.sub2 <- subset(data.sub, chr == nmi_bladder$chromosome[[i]])
    data.sub3 <- subset(data.sub, pos %in% nmi_bladder$lowerbound[[i]]:nmi_bladder$upperbound[[i]])
    nmi_count <- nmi_count + nrow(data.sub3)
    both_count <- both_count + nrow(data.sub3)
  }

  mi_count <- mi_count / length(unique(data.sub$samples))
  mi_hr_count <- c(mi_hr_count, mi_count)
  nmi_count <- nmi_count / length(unique(data.sub$samples))
  nmi_hr_count <- c(nmi_hr_count, nmi_count)
  both_count <- both_count / length(unique(data.sub$samples))
  both_hr_count <- c(both_hr_count, both_count)
}


mi_lr_count <- c()
nmi_lr_count <- c()
both_lr_count <- c()
for (s in lr_samples){
  data.sub <- subset(nu_lowrisk, patient == s)
  mi_count <- 0
  nmi_count <- 0
  both_count <- 0
  for (i in 1:nrow(bladder)){
    data.sub2 <- subset(data.sub, chr == bladder$chromosome[[i]])
    data.sub3 <- subset(data.sub, pos %in% bladder$lowerbound[[i]]:bladder$upperbound[[i]])
    mi_count <- mi_count + nrow(data.sub3)
    both_count <- both_count + nrow(data.sub3)
  }

  for (i in 1:nrow(nmi_bladder)){
    data.sub2 <- subset(data.sub, chr == nmi_bladder$chromosome[[i]])
    data.sub3 <- subset(data.sub, pos %in% nmi_bladder$lowerbound[[i]]:nmi_bladder$upperbound[[i]])
    nmi_count <- nmi_count + nrow(data.sub3)
    both_count <- both_count + nrow(data.sub3)
  }

  mi_count <- mi_count / length(unique(data.sub$samples))
  mi_lr_count <- c(mi_lr_count, mi_count)
  nmi_count <- nmi_count / length(unique(data.sub$samples))
  nmi_lr_count <- c(nmi_lr_count, nmi_count)
  both_count <- both_count / length(unique(data.sub$samples))
  both_lr_count <- c(both_lr_count, both_count)
}

nu_1 <- data.frame("NMIBC_Count" = c(nmi_hr_count, nmi_lr_count),
                   "MIBC_Count" = c(mi_hr_count, mi_lr_count), "All_Count" = c(both_hr_count, both_lr_count),
                      "patient" = c(hr_samples, lr_samples),
                   "Risk_Type" = c(rep("HR", length(hr_samples)), rep("LR", length(lr_samples))))
nu_simp <- unique(nu_all[,c(1,5:7)])
df <- merge(nu_1, nu_simp, by = "patient")



#all <- df[,-1]

#all <- df[,-c(1:4)]
# set the seed, and put aside a test set

library(randomForest)
library(pROC)
library(caret)
library(leaps)
library(glmnet)

############### logistic regression

set.seed(6)

df.log <- df
df.log$Risk_Type <- as.numeric(as.factor(df.log$Risk_Type)) - 1
#df.log$smoking <- as.numeric(as.factor(df.log$smoking)) - 1
df.log$gender <- as.numeric(as.factor(df.log$gender)) - 1
df.log$age <- sqrt(max(df.log$age + 1) - df.log$age)

non <- c()
ex <- c()
current <- c()
for (i in 1:nrow(df)){
  if (df$smoking[[i]] == "Non-smoker"){
    non <- c(non, 1)
    ex <- c(ex, 0)
    current <- c(current, 0)
  }
  if (df$smoking[[i]] == "Ex-smoker"){
    non <- c(non, 0)
    ex <- c(ex, 1)
    current <- c(current, 0)
  }
  if (df$smoking[[i]] == "Smoker"){
    non <- c(non, 0)
    ex <- c(ex, 0)
    current <- c(current, 1)
  }
}

df.log$current_smoker <- unlist(current)
df.log$ex_smoker <- unlist(ex)
df.log$non_smoker <- unlist(non)


df.log <- df.log[,c(2:4,6,7,9:11)]
df.log$Risk_Type <- as.numeric(factor(df$Risk_Type)) - 1

mibc_only <- as.data.frame(scale(df.log[,c(2,4:8)]))
mibc_only$Risk_Type <- df.log$Risk_Type
nmibc_only <- as.data.frame(scale(df.log[,c(1,4:8)]))
nmibc_only$Risk_Type <- df.log$Risk_Type
all_bc <- as.data.frame(scale(df.log[,c(3:8)]))
all_bc$Risk_Type <- df.log$Risk_Type
clin_only <- as.data.frame(scale(df.log[,c(4:8)]))
clin_only$Risk_Type <- df.log$Risk_Type
mut_only <- as.data.frame(scale(df.log[,c(1:3,8)]))
mut_only$Risk_Type <- df.log$Risk_Type
both <- as.data.frame(scale(df.log[,c(1:2,4:8)]))
both$Risk_Type <- df.log$Risk_Type

test_indis <- sample(1:nrow(df.log), .30*nrow(df.log))
test.mibc <- mibc_only[test_indis, ]
training.mibc <- mibc_only[-test_indis, ]
set.seed(6)

fit <- glm(Risk_Type ~ ., data = training.mibc)
y_hat <- predict(fit, newdata = test.mibc)
y_true <- test.mibc$Risk_Type

lr.mibc_auc <- auc(y_true, y_hat)
lr.mibc_rocobj <- roc(y_true, y_hat)

test.nmibc <- nmibc_only[test_indis, ]
training.nmibc <- nmibc_only[-test_indis, ]


fit <- glm(Risk_Type ~ ., data = training.nmibc)
y_hat <- predict(fit, newdata = test.nmibc)
y_true <- test.nmibc$Risk_Type

lr.nmibc_auc <- auc(y_true, y_hat)
lr.nmibc_rocobj <- roc(y_true, y_hat)


test.clin <- clin_only[test_indis, ]
training.clin <- clin_only[-test_indis, ]


fit <- glm(Risk_Type ~ ., data = training.clin)
y_hat <- predict(fit, newdata = test.clin)
y_true <- test.clin$Risk_Type

lr.clin_auc <- auc(y_true, y_hat)
lr.clin_rocobj <- roc(y_true, y_hat)

test <- all_bc[test_indis, ]
training <- all_bc[-test_indis, ]


fit <- glm(Risk_Type ~ ., data = training)
y_hat <- predict(fit, newdata = test)
y_true <- test$Risk_Type

lr.all_auc <- auc(y_true, y_hat)
lr.all_rocobj <- roc(y_true, y_hat)

test.mut <- mut_only[test_indis, ]
training.mut <- mut_only[-test_indis, ]


fit <- glm(Risk_Type ~ ., data = training.mut)
y_hat <- predict(fit, newdata = test.mut)
y_true <- test.mut$Risk_Type

lr.mut_auc <- auc(y_true, y_hat)
lr.mut_rocobj <- roc(y_true, y_hat)

test.both <- both[test_indis, ]
training.both <- both[-test_indis, ]


fit <- glm(Risk_Type ~ ., data = training.both)
y_hat <- predict(fit, newdata = test.both)
y_true <- test.both$Risk_Type

lr.both_auc <- auc(y_true, y_hat)
lr.both_rocobj <- roc(y_true, y_hat)

lr.mibc_auc
lr.nmibc_auc
lr.clin_auc
lr.all_auc
lr.mut_auc
lr.both_auc

############# neural network

library(neuralnet)
set.seed(6)
nn <- neuralnet(Risk_Type~., data = training.mibc)
y_hat <- predict(nn, newdata = test.mibc)
y_true <- test.mibc$Risk_Type

nn.mibc_auc <- auc(y_true, as.vector(y_hat))
nn.mibc_rocobj <- roc(y_true, y_hat)


nn <- neuralnet(Risk_Type~., data = training.nmibc)
y_hat <- predict(nn, newdata = test.nmibc)
y_true <- test.nmibc$Risk_Type

nn.nmibc_auc <- auc(y_true, y_hat)
nn.nmibc_rocobj <- roc(y_true, y_hat)


nn <- neuralnet(Risk_Type~., data = training.clin)
y_hat <- predict(nn, newdata = test.clin)
y_true <- test.clin$Risk_Type

nn.clin_auc <- auc(y_true, y_hat)
nn.clin_rocobj <- roc(y_true, y_hat)

nn <- neuralnet(Risk_Type~., data = training)
y_hat <- predict(nn, newdata = test)
y_true <- test$Risk_Type

nn.all_auc <- auc(y_true, y_hat)
nn.all_rocobj <- roc(y_true, y_hat)

nn <- neuralnet(Risk_Type~., data = training.mut)
y_hat <- predict(nn, newdata = test.mut)
y_true <- test.mut$Risk_Type

nn.mut_auc <- auc(y_true, y_hat)
nn.mut_rocobj <- roc(y_true, y_hat)


nn <- neuralnet(Risk_Type~., data = training.both)
y_hat <- predict(nn, newdata = test.both)
y_true <- test.both$Risk_Type

nn.both_auc <- auc(y_true, y_hat)
nn.both_rocobj <- roc(y_true, y_hat)

nn.mibc_auc
nn.nmibc_auc
nn.clin_auc
nn.all_auc
nn.mut_auc
nn.both_auc

################## random forest
set.seed(6)
rf.fit <- randomForest(Risk_Type~., data = training.mibc, n.tree = 10000)
y_hat <- predict(rf.fit, newdata = test.mibc, type = "response")
y_true <- test.mibc$Risk_Type
mibc.imp <- importance(rf.fit)
rf.mibc_auc <- auc(y_true, y_hat)
rf.mibc_rocobj <- roc(y_true, y_hat)

set.seed(6)
rf.fit <- randomForest(Risk_Type~., data = training.nmibc, n.tree = 10000)
y_hat <- predict(rf.fit, newdata = test.nmibc, type = "response")
y_true <- test.nmibc$Risk_Type
nmibc.imp <- importance(rf.fit)
rf.nmibc_auc <- auc(y_true, y_hat)
rf.nmibc_rocobj <- roc(y_true, y_hat)

set.seed(6)
rf.fit <- randomForest(Risk_Type~., data = training.clin, n.tree = 10000)
y_hat <- predict(rf.fit, newdata = test.clin, type = "response")
y_true <- test.clin$Risk_Type
clin.imp <- importance(rf.fit)
rf.clin_auc <- auc(y_true, y_hat)
rf.clin_rocobj <- roc(y_true, y_hat)

set.seed(6)
rf.fit <- randomForest(Risk_Type~., data = training, n.tree = 10000)
y_hat <- predict(rf.fit, newdata = test, type = "response")
y_true <- test$Risk_Type
varImpPlot(rf.fit)
rf.all_auc <- auc(y_true, y_hat)
rf.all_rocobj <- roc(y_true, y_hat)

set.seed(6)
rf.fit <- randomForest(Risk_Type~., data = training.mut, n.tree = 10000)
y_hat <- predict(rf.fit, newdata = test.mut, type = "response")
y_true <- test.mut$Risk_Type
varImpPlot(rf.fit)
rf.mut_auc <- auc(y_true, y_hat)
rf.mut_rocobj <- roc(y_true, y_hat)

set.seed(6)
rf.fit <- randomForest(Risk_Type~., data = training.both, n.tree = 10000)
y_hat <- predict(rf.fit, newdata = test.both, type = "response")
y_true <- test.both$Risk_Type
both.imp <- importance(rf.fit)
rf.both_auc <- auc(y_true, y_hat)
rf.both_rocobj <- roc(y_true, y_hat)

rf.mibc_auc
rf.nmibc_auc
rf.clin_auc
rf.all_auc
rf.mut_auc
rf.both_auc


#### graphing


# mibc, nmibc, clin, both


mibc.roc_df <- data.frame("Specificity" = c(rev(lr.mibc_rocobj$specificities), rev(nn.mibc_rocobj$specificities), rev(rf.mibc_rocobj$specificities)),
                     "Sensitivity" = c(rev(lr.mibc_rocobj$sensitivities), rev(nn.mibc_rocobj$sensitivities), rev(rf.mibc_rocobj$sensitivities)),
                     "Model" = c(rep("Logistic Regression", length(lr.mibc_rocobj$sensitivities)), rep("Neural Network", length(nn.mibc_rocobj$sensitivities)),
                                 rep("Random Forest", length(rf.mibc_rocobj$sensitivities))))


mibc.roc_plot <- ggplot(mibc.roc_df, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(aes(color = Model), linewidth = 1) +
  scale_x_reverse() +
  theme_minimal() +
  scale_color_manual(values = c("#E7B800", "#FC4E07", "#00AFBB")) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  geom_segment(
    aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed", size = 1.25) +
  ggtitle("MIBC Hotspots and Personal Risk Factors") +
  annotate("text", x = 0.25, y = 0.375, col = "#E7B800", label = "AUC = 0.8929", size = 4) +
  annotate("text", x = 0.25, y = 0.3, col = "#FC4E07", label = "AUC = 0.8643", size = 4) +
  annotate("text", x = 0.25, y = 0.225, col = "#00AFBB", label = "AUC = 0.9286", size = 4)

nmibc.roc_df <- data.frame("Specificity" = c(rev(lr.nmibc_rocobj$specificities), rev(nn.nmibc_rocobj$specificities), rev(rf.nmibc_rocobj$specificities)),
                          "Sensitivity" = c(rev(lr.nmibc_rocobj$sensitivities), rev(nn.nmibc_rocobj$sensitivities), rev(rf.nmibc_rocobj$sensitivities)),
                          "Model" = c(rep("Logistic Regression", length(lr.nmibc_rocobj$sensitivities)), rep("Neural Network", length(nn.nmibc_rocobj$sensitivities)),
                                      rep("Random Forest", length(rf.nmibc_rocobj$sensitivities))))


nmibc.roc_plot <- ggplot(nmibc.roc_df, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(aes(color = Model), linewidth = 1) +
  scale_x_reverse() +
  theme_minimal() +
  scale_color_manual(values = c("#E7B800", "#FC4E07", "#00AFBB")) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  geom_segment(
    aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed", size = 1.25) +
  ggtitle("NMIBC Hotspots and Personal Risk Factors") +
  annotate("text", x = 0.25, y = 0.375, col = "#E7B800", label = "AUC = 0.8786", size = 4) +
  annotate("text", x = 0.25, y = 0.3, col = "#FC4E07", label = "AUC = 0.8643", size = 4) +
  annotate("text", x = 0.25, y = 0.225, col = "#00AFBB", label = "AUC = 0.9000", size = 4)


clin.roc_df <- data.frame("Specificity" = c(rev(lr.clin_rocobj$specificities), rev(nn.clin_rocobj$specificities), rev(rf.clin_rocobj$specificities)),
                          "Sensitivity" = c(rev(lr.clin_rocobj$sensitivities), rev(nn.clin_rocobj$sensitivities), rev(rf.clin_rocobj$sensitivities)),
                          "Model" = c(rep("Logistic Regression", length(lr.clin_rocobj$sensitivities)), rep("Neural Network", length(nn.clin_rocobj$sensitivities)),
                                      rep("Random Forest", length(rf.clin_rocobj$sensitivities))))


clin.roc_plot <- ggplot(clin.roc_df, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(aes(color = Model), linewidth = 1) +
  scale_x_reverse() +
  theme_minimal() +
  scale_color_manual(values = c("#E7B800", "#FC4E07", "#00AFBB")) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  geom_segment(
    aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed", size = 1.25) +
  ggtitle("Personal Risk Factors Only") +
  annotate("text", x = 0.25, y = 0.375, col = "#E7B800", label = "AUC = 0.8786", size = 4) +
  annotate("text", x = 0.25, y = 0.3, col = "#FC4E07", label = "AUC = 0.8429", size = 4) +
  annotate("text", x = 0.25, y = 0.225, col = "#00AFBB", label = "AUC = 0.8964", size = 4)



mut.roc_df <- data.frame("Specificity" = c(rev(lr.mut_rocobj$specificities), rev(nn.mut_rocobj$specificities), rev(rf.mut_rocobj$specificities)),
                          "Sensitivity" = c(rev(lr.mut_rocobj$sensitivities), rev(nn.mut_rocobj$sensitivities), rev(rf.mut_rocobj$sensitivities)),
                          "Model" = c(rep("Logistic Regression", length(lr.mut_rocobj$sensitivities)), rep("Neural Network", length(nn.mut_rocobj$sensitivities)),
                                      rep("Random Forest", length(rf.mut_rocobj$sensitivities))))


mut.roc_plot <- ggplot(mut.roc_df, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(aes(color = Model), linewidth = 1) +
  scale_x_reverse() +
  theme_minimal() +
  scale_color_manual(values = c("#E7B800", "#FC4E07", "#00AFBB")) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  geom_segment(
    aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed", size = 1.25) +
  ggtitle("MIBC and NMIBC Hotspots") +
  annotate("text", x = 0.25, y = 0.375, col = "#E7B800", label = "AUC = 0.6964", size = 4) +
  annotate("text", x = 0.25, y = 0.3, col = "#FC4E07", label = "AUC = 0.7036", size = 4) +
  annotate("text", x = 0.25, y = 0.225, col = "#00AFBB", label = "AUC = 0.7429", size = 4)











##########################
##########################
##########################
##########################





all_cor <- df[,-1]
colnames(all_cor) <- c("NMIBC_Count", "MIBC_Count", "Both", "Risk_Type", "Gender", "Age", "Smoking_Status")
all_cor$Smoking_Status <- as.numeric(as.factor(all_cor$Smoking_Status))
all_cor$Gender <- as.numeric(as.factor(all_cor$Gender))
all_cor$Risk_Type <- as.numeric(as.factor(all_cor$Risk_Type))


library(corrplot)
correlations <- cor(all_cor[,-3])
cor.plot <- corrplot(correlations, tl.cex = 1.25)


mibc.vip_df <- data.frame("Variable" = c("Current Smoker","Non-Smoker",
                                    " MIBC Hotspots",  "Ex-Smoker", "Gender",
                                     "Age"),
                     "MeanDecreaseGini" = c(0.2695625, 0.5272787, 0.5874977, 0.8169978, 0.8606566, 1.8862795))

mibc.vip_df$Variable <- factor(mibc.vip_df$Variable, levels = mibc.vip_df$Variable)

mibc.vip_plot <- ggplot(mibc.vip_df, aes(x = MeanDecreaseGini, y = Variable)) +
  geom_point(color="#00AFBB", fill = "#00AFBB", size=4, shape = 23) +
  ylab(NULL) +
  xlab("Mean Decrease Gini") +
  theme_classic()+
  ggtitle("MIBC Hotspots and Personal Risk Factors") +
  xlim(0,2) +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12, angle = 90),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12, angle = 25), title = element_text(size = 12),
        plot.title= element_text(size = 12))
mibc.vip_plot

 
nmibc.vip_df <- data.frame("Variable" = c("NMIBC Hotspots","Current Smoker","Non-Smoker","Gender", "Ex-Smoker",
                                           
                                         "Age"),
                          "MeanDecreaseGini" = c(0.03294403, 0.34170427, 0.50785108, 0.73431767, 0.88472890, 1.92890221))

nmibc.vip_df$Variable <- factor(nmibc.vip_df$Variable, levels = nmibc.vip_df$Variable)

nmibc.vip_plot <- ggplot(nmibc.vip_df, aes(x = MeanDecreaseGini, y = Variable)) +
  geom_point(color="#00AFBB", fill = "#00AFBB", size=4, shape = 23) +
  ylab(NULL) +
  xlab("Mean Decrease Gini") +
  theme_classic()+
  xlim(0,2) +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12, angle = 90),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12, angle = 25), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  ggtitle("NMIBC Hotspots and Personal Risk Factors") 
nmibc.vip_plot

clin.vip_df <- data.frame("Variable" = c("    Current Smoker","Gender", "Non-Smoker", "Ex-Smoker",
                                          "Age"),
                           "MeanDecreaseGini" = c(0.2545370, 0.4657496, 0.4926756, 0.6720723, 1.0346603))

clin.vip_df$Variable <- factor(clin.vip_df$Variable, levels = clin.vip_df$Variable)

clin.vip_plot <- ggplot(clin.vip_df, aes(x = MeanDecreaseGini, y = Variable)) +
  geom_point(color="#00AFBB", fill = "#00AFBB", size=4, shape = 23) +
  ylab(NULL) +
  xlab("Mean Decrease Gini") +
  theme_classic()+
  xlim(0,2) +
  theme(axis.text.x = element_text(vjust = 0.25, hjust=0.25, size = 12, angle = 90),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12, angle = 25), title = element_text(size = 12),
        plot.title= element_text(size = 12)) +
  ggtitle("Personal Risk Factors Only") 
clin.vip_plot


library(ggpubr)

b <- ggarrange(mibc.roc_plot, nmibc.roc_plot, clin.roc_plot, nrow = 3, ncol = 1, labels = c("A"))

c <- ggarrange(mibc.vip_plot, nmibc.vip_plot, clin.vip_plot, nrow = 3, ncol = 1, labels = c("B"))

fig3 <- ggarrange(b, c, ncol = 2)

fig3

