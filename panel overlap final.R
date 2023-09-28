

# check MIBC10 and NMIBC10 panel overlap



mi <- read.csv("E:\\Bladder Hotspot Paper\\Final Analysis\\MIBC_Bins_Final.csv")


nmi <- read.csv("E:\\Bladder Hotspot Paper\\Final Analysis\\NMIBC_Bins_Final.csv")





both <- rbind(mi, nmi)


mi_unique <- data.frame()
nmi_unique <- data.frame()
shared <- data.frame()

for (i in 1:nrow(both)){
  print(i)
  mi_sub <- subset(mi, chromosome == both$chromosome[[i]])
  mi_pos <- c()
  if (nrow(mi_sub) > 0){
  for (n in 1:nrow(mi_sub)){
    mi_pos <- c(mi_pos, mi_sub$lowerbound[[n]]:mi_sub$upperbound[[n]])
  }
  }

  nmi_sub <- subset(nmi, chromosome == both$chromosome[[i]])
  nmi_pos <- c()
  if (nrow(nmi_sub) > 0){
  for (m in 1:nrow(nmi_sub)){
    nmi_pos <- c(nmi_pos, nmi_sub$lowerbound[[m]]:nmi_sub$upperbound[[m]])
  }
  }


  if (length(intersect(mi_pos, both$lowerbound[[i]]:both$upperbound[[i]])) >= 1 & length(intersect(nmi_pos, both$lowerbound[[i]]:both$upperbound[[i]])) >= 1){
    shared <- rbind(shared, both[i,])
  }
  if (length(intersect(mi_pos, both$lowerbound[[i]]:both$upperbound[[i]])) < 1 & length(intersect(nmi_pos, both$lowerbound[[i]]:both$upperbound[[i]])) >= 1){
    nmi_unique <- rbind(nmi_unique, both[i,])
  }
  if (length(intersect(mi_pos, both$lowerbound[[i]]:both$upperbound[[i]])) >= 1 & length(intersect(nmi_pos, both$lowerbound[[i]]:both$upperbound[[i]])) < 1){
    mi_unique <- rbind(mi_unique, both[i,])
  }
}

unique_hotspots <- data.frame()
for (i in 1:nrow(shared)){
  shared.new <- shared
  sub <- subset(shared.new, chromosome == shared$chromosome[[i]])
  sub2 <- subset(sub, lowerbound %in% shared$lowerbound[[i]]:shared$upperbound[[i]] | upperbound %in% shared$lowerbound[[i]]:shared$upperbound[[i]])
  new_hotspot <- data.frame("lowerbound" = min(sub2$lowerbound), "upperbound" = max(sub2$upperbound), "chromosome" = shared$chromosome[[i]])
  unique_hotspots <- rbind(unique_hotspots, new_hotspot)
}


unique_hotspots <- unique(unique_hotspots)


nmi_unique_hotspots <- data.frame()
for (i in 1:nrow(nmi_unique)){
  shared.new <- nmi_unique
  sub <- subset(shared.new, chromosome == nmi_unique$chromosome[[i]])
  sub2 <- subset(sub, lowerbound %in% nmi_unique$lowerbound[[i]]:nmi_unique$upperbound[[i]] | upperbound %in% nmi_unique$lowerbound[[i]]:nmi_unique$upperbound[[i]])
  new_hotspot <- data.frame("lowerbound" = min(sub2$lowerbound), "upperbound" = max(sub2$upperbound), "chromosome" = nmi_unique$chromosome[[i]])
  nmi_unique_hotspots <- rbind(nmi_unique_hotspots, new_hotspot)
}


nmi_unique_hotspots <- unique(nmi_unique_hotspots)




mi_unique_hotspots <- data.frame()
for (i in 1:nrow(mi_unique)){
  shared.new <- mi_unique
  sub <- subset(shared.new, chromosome == mi_unique$chromosome[[i]])
  sub2 <- subset(sub, lowerbound %in% mi_unique$lowerbound[[i]]:mi_unique$upperbound[[i]] | upperbound %in% mi_unique$lowerbound[[i]]:mi_unique$upperbound[[i]])
  new_hotspot <- data.frame("lowerbound" = min(sub2$lowerbound), "upperbound" = max(sub2$upperbound), "chromosome" = mi_unique$chromosome[[i]])
  mi_unique_hotspots <- rbind(mi_unique_hotspots, new_hotspot)
}


mi_unique_hotspots <- unique(mi_unique_hotspots)









