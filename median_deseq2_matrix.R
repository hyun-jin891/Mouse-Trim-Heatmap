mat <- readRDS("/home/xkdl27/mouse_DB/heat_map/normalized_matrix.rds")

sampleNum <- c(8, 13, 11, 1, 4, 11, 4, 3, 12, 2, 5, 9, 9, 9, 5, 2, 7, 8, 7, 9, 7, 17, 4, 7, 13, 4, 1,
               10, 9, 2, 7, 13, 9, 12, 6, 7, 6, 5, 9, 9, 6, 10, 10, 1, 10, 4, 6, 10, 4, 2, 2, 6, 
               8, 25, 11, 4, 9, 10, 8, 10, 7, 2, 2, 3, 7)

curIndex <- 1
print(curIndex + sampleNum[1] - 1)
inx <- curIndex + sampleNum[1] - 1
median_mat <- mat[, curIndex: inx]

temp <- matrix(0, length(median_mat[,1]), 1)

for (i in 1:(length(rownames(median_mat)))){
  median_temp <- median(as.numeric(median_mat[i,]))
  temp[i] <- median_temp
}

name_mat <- readRDS("/home/xkdl27/mouse_DB/heat_map/mouse_merged_TRIM_matrix_final.rds")

sampleName <- colnames(name_mat)

colnames(temp) <- sampleName[1]

median_mat <- temp
rownames(median_mat) <- rownames(mat)

curIndex <- curIndex + sampleNum[1]

for (i in 2:length(sampleNum)){
  inx <- curIndex + sampleNum[i] - 1
  median_mat_temp <- mat[, curIndex : inx, drop=FALSE]
  temp <- matrix(0, length(median_mat_temp[,1]), 1)
  
  for (j in 1:length(rownames(median_mat_temp))){
    median_temp <- median(as.numeric(median_mat_temp[j,]))
    temp[j] <- median_temp
  }
  
  colnames(temp) <- sampleName[i]
  median_mat <- cbind(median_mat, temp)
  rownames(median_mat) <- rownames(mat)
  
  curIndex <- curIndex + sampleNum[i]
  
}


saveRDS(median_mat, file = "/home/xkdl27/mouse_DB/heat_map/median_expression_matrix.rds")

