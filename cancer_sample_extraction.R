####### Bladder cancer


raw_BDCA <- read.table("mouse_DB/cancer_expression_DB/Bladder_cancer_expression_matrix.tsv", sep="\t", header=FALSE)


raw_BDCA <- as.matrix(raw_BDCA[,1:5])

rownames(raw_BDCA) <- raw_BDCA[, 1]

raw_BDCA <- raw_BDCA[, -1]
colnames(raw_BDCA) <- raw_BDCA[1,]
raw_BDCA <- raw_BDCA[-1,]

raw_BDCA <- matrix(apply(raw_BDCA, 2, as.integer), nrow=nrow(raw_BDCA), dimnames=dimnames(raw_BDCA))

rowname <- rownames(raw_BDCA)


saveRDS(raw_BDCA, file="/home/xkdl27/mouse_DB/new_normalization_output/BDCA.rds")



####### BRCA

raw_BRCA <- read.table("mouse_DB/cancer_expression_DB/BRCA_expression_matrix.tsv", sep="\t", header=FALSE)


colnames(raw_BRCA) <- raw_BRCA[1,]
raw_BRCA <- raw_BRCA[-1,]
rowname <- raw_BRCA[, 1]

raw_BRCA <- as.matrix(raw_BRCA[, colnames(raw_BRCA) %in% c("GSM2144727", "GSM2144729", "GSM2144735", "GSM2144736", "GSM2144745", "GSM3502135")])

rownames(raw_BRCA) <- rowname


raw_BRCA <- matrix(apply(raw_BRCA, 2, as.integer), nrow=nrow(raw_BRCA), dimnames=dimnames(raw_BRCA))

saveRDS(raw_BRCA, file="/home/xkdl27/mouse_DB/new_normalization_output/BRCA.rds")



####### COAD

raw_COAD <- read.table("mouse_DB/cancer_expression_DB/Colon_Adenocarcinoma_expression_matrix.tsv", sep="\t", header=FALSE)



colnames(raw_COAD) <- raw_COAD[1,]
raw_COAD <- raw_COAD[-1,]
rowname <- raw_COAD[, 1]

raw_COAD <- as.matrix(raw_COAD[, colnames(raw_COAD) %in% c("GSM4336388", "GSM5057361", "GSM5057362", "GSM5057363")])

rownames(raw_COAD) <- rowname


raw_COAD <- matrix(apply(raw_COAD, 2, as.integer), nrow=nrow(raw_COAD), dimnames=dimnames(raw_COAD))

saveRDS(raw_COAD, file="/home/xkdl27/mouse_DB/new_normalization_output/COAD.rds")

####### GBM

raw_GBM <- read.table("mouse_DB/cancer_expression_DB/GBM_expression_matrix.tsv", sep="\t", header=FALSE)

colnames(raw_GBM) <- raw_GBM[1,]
raw_GBM <- raw_GBM[-1,]
rowname <- raw_GBM[, 1]

raw_GBM <- as.matrix(raw_GBM[, colnames(raw_GBM) %in% c("GSM3502792", "GSM3502800")])

rownames(raw_GBM) <- rowname
raw_GBM <- matrix(apply(raw_GBM, 2, as.integer), nrow=nrow(raw_GBM), dimnames=dimnames(raw_GBM))
saveRDS(raw_GBM, file="/home/xkdl27/mouse_DB/new_normalization_output/GBM.rds")

####### LUAD

raw_LUAD <- read.table("/home/xkdl27/mouse_DB/cancer_expression_DB/LUAD_expression_matrix.tsv", sep="\t", header=FALSE)

colnames(raw_LUAD) <- raw_LUAD[1,]
raw_LUAD <- raw_LUAD[-1,]
rowname <- raw_LUAD[, 1]
raw_LUAD <- raw_LUAD[,-1]


raw_LUAD <- as.matrix(raw_LUAD)

rownames(raw_LUAD) <- rowname


raw_LUAD <- matrix(apply(raw_LUAD, 2, as.integer), nrow=nrow(raw_LUAD), dimnames=dimnames(raw_LUAD))


saveRDS(raw_LUAD, file="/home/xkdl27/mouse_DB/new_normalization_output/LUAD.rds")

####### PDAC

raw_PDAC <- read.table("mouse_DB/cancer_expression_DB/PDAC_expression_matrix.tsv", sep="\t", header=FALSE)

colnames(raw_PDAC) <- raw_PDAC[1,]
raw_PDAC <- raw_PDAC[-1,]
rowname <- raw_PDAC[, 1]
raw_PDAC <- raw_PDAC[,-1]


raw_PDAC <- as.matrix(raw_PDAC[, colnames(raw_PDAC) %in% c("GSM2409486", "GSM2409487", "GSM2414654", "GSM2414655", "GSM2414657", "GSM2414658", "GSM2414660", "GSM2414662",	"GSM2414663", "GSM2414664")])

rownames(raw_PDAC) <- rowname


raw_PDAC <- matrix(apply(raw_PDAC, 2, as.integer), nrow=nrow(raw_PDAC), dimnames=dimnames(raw_PDAC))



saveRDS(raw_PDAC, file="/home/xkdl27/mouse_DB/new_normalization_output/PDAC.rds")

####### LUSC

raw_LUSC <- read.table("mouse_DB/cancer_expression_DB/lung_squamous_cell_carcinoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_LUSC)
colnames(raw_LUSC) <- raw_LUSC[1,]
raw_LUSC <- raw_LUSC[-1,]
rowname <- raw_LUSC[, 1]
raw_LUSC <- raw_LUSC[,-1]


raw_LUSC <- as.matrix(raw_LUSC[, colnames(raw_LUSC) %in% c("GSM3510091", "GSM3510092", "GSM3510095", "GSM3510096", "GSM3510097", "GSM4495271", "GSM4495272", "GSM4495273", "GSM4495274", "GSM4495275",  "GSM4495276")])

rownames(raw_LUSC) <- rowname


raw_LUSC <- matrix(apply(raw_LUSC, 2, as.integer), nrow=nrow(raw_LUSC), dimnames=dimnames(raw_LUSC))

raw_LUSC <- normalized_expression_matrix2(raw_LUSC)

new_raw_LUSC <- matrix(0, length(rowname), 1)
colnames(new_raw_LUSC) <- "LUSC"
rownames(new_raw_LUSC) <- rowname

for (i in 1:length(rowname)){
  new_raw_LUSC[i, 1] <- mean(raw_LUSC[i,])
}

saveRDS(raw_LUSC, file="/home/xkdl27/mouse_DB/new_normalization_output/LUSC.rds")

####### OV

raw_OV <- read.table("mouse_DB/cancer_expression_DB/ovarian_cancer_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_OV)
colnames(raw_OV) <- raw_OV[1,]
raw_OV <- raw_OV[-1,]
rowname <- raw_OV[, 1]
raw_OV <- raw_OV[,-1]


raw_OV <- as.matrix(raw_OV[, colnames(raw_OV) %in% c("GSM3693233", "GSM3693234", "GSM3693235", "GSM3693236", "GSM3693237", "GSM3693238", "GSM3693240", "GSM3693241")])

rownames(raw_OV) <- rowname


raw_OV <- matrix(apply(raw_OV, 2, as.integer), nrow=nrow(raw_OV), dimnames=dimnames(raw_OV))

raw_OV <- normalized_expression_matrix2(raw_OV)

new_raw_OV <- matrix(0, length(rowname), 1)
colnames(new_raw_OV) <- "OV"
rownames(new_raw_OV) <- rowname

for (i in 1:length(rowname)){
  new_raw_OV[i, 1] <- mean(raw_OV[i,])
}
saveRDS(raw_OV, file="/home/xkdl27/mouse_DB/new_normalization_output/OV.rds")

######## prostate cancer



raw_pcancer <- read.table("mouse_DB/cancer_expression_DB/prostate_cancer_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_pcancer)
colnames(raw_pcancer) <- raw_pcancer[1,]
raw_pcancer <- raw_pcancer[-1,]
rowname <- raw_pcancer[, 1]
raw_pcancer <- raw_pcancer[,-1]


raw_pcancer <- as.matrix(raw_pcancer[, colnames(raw_pcancer) %in% c("GSM3070579", "GSM3070580", "GSM3070581", "GSM3070582", "GSM3070583", "GSM3070584", "GSM3592433")])

rownames(raw_pcancer) <- rowname


raw_pcancer <- matrix(apply(raw_pcancer, 2, as.integer), nrow=nrow(raw_pcancer), dimnames=dimnames(raw_pcancer))

raw_pcancer <- normalized_expression_matrix2(raw_pcancer)

new_raw_pcancer <- matrix(0, length(rowname), 1)
colnames(new_raw_pcancer) <- "PCa"
rownames(new_raw_pcancer) <- rowname

for (i in 1:length(rowname)){
  new_raw_pcancer[i, 1] <- mean(raw_pcancer[i,])
}



saveRDS(raw_pcancer, file="/home/xkdl27/mouse_DB/new_normalization_output/prostate_cancer.rds")


######## THCA



raw_THCA <- read.table("mouse_DB/cancer_expression_DB/thyroid_carcinoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_THCA)
colnames(raw_THCA) <- raw_THCA[1,]
raw_THCA <- raw_THCA[-1,]
rowname <- raw_THCA[, 1]
raw_THCA <- raw_THCA[,-1]


raw_THCA <- as.matrix(raw_THCA[, colnames(raw_THCA) %in% c("GSM3302918", "GSM3302919", "GSM3302920", "GSM3302921", "GSM3302922", "GSM3302923", "GSM3302924" )])

rownames(raw_THCA) <- rowname


raw_THCA <- matrix(apply(raw_THCA, 2, as.integer), nrow=nrow(raw_THCA), dimnames=dimnames(raw_THCA))

raw_THCA <- normalized_expression_matrix2(raw_THCA)

new_raw_THCA <- matrix(0, length(rowname), 1)
colnames(new_raw_THCA) <- "THCA"
rownames(new_raw_THCA) <- rowname

for (i in 1:length(rowname)){
  new_raw_THCA[i, 1] <- mean(raw_THCA[i,])
}
saveRDS(raw_THCA, file="/home/xkdl27/mouse_DB/new_normalization_output/THCA.rds")

######## RCC


raw_RCC <- read.table("mouse_DB/cancer_expression_DB/Renal_Cell_Carcinoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_RCC)
colnames(raw_RCC) <- raw_RCC[1,]
raw_RCC <- raw_RCC[-1,]
rowname <- raw_RCC[, 1]
raw_RCC <- raw_RCC[,-1]


raw_RCC <- as.matrix(raw_RCC[, colnames(raw_RCC) %in% c("GSM4933264", "GSM4933265" )])

rownames(raw_RCC) <- rowname


raw_RCC <- matrix(apply(raw_RCC, 2, as.integer), nrow=nrow(raw_RCC), dimnames=dimnames(raw_RCC))

raw_RCC <- normalized_expression_matrix2(raw_RCC)

new_raw_RCC <- matrix(0, length(rowname), 1)
colnames(new_raw_RCC) <- "RCC"
rownames(new_raw_RCC) <- rowname

for (i in 1:length(rowname)){
  new_raw_RCC[i, 1] <- mean(raw_RCC[i,])
}

saveRDS(raw_RCC, file="/home/xkdl27/mouse_DB/new_normalization_output/RCC.rds")


######## LPS


raw_LPS <- read.table("mouse_DB/cancer_expression_DB/Liposarcoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_LPS)


colnames(raw_LPS) <- raw_LPS[1,]
raw_LPS <- raw_LPS[-1,]
rowname <- raw_LPS[, 1]
raw_LPS <- raw_LPS[,-1]


raw_LPS <- as.matrix(raw_LPS[, colnames(raw_LPS) %in% c("GSM2127221", "GSM2127222", "GSM2127223", "GSM2127224", "GSM2127225", "GSM2127226", "GSM2127227", "GSM2127228")])


rownames(raw_LPS) <- rowname


raw_LPS <- matrix(apply(raw_LPS, 2, as.integer), nrow=nrow(raw_LPS), dimnames=dimnames(raw_LPS))

raw_LPS <- normalized_expression_matrix2(raw_LPS)

new_raw_LPS <- matrix(0, length(rowname), 1)
colnames(new_raw_LPS) <- "LPS"
rownames(new_raw_LPS) <- rowname

for (i in 1:length(rowname)){
  new_raw_LPS[i, 1] <- mean(raw_LPS[i,])
}

saveRDS(raw_LPS, file="/home/xkdl27/mouse_DB/new_normalization_output/LPS.rds")

######## HCC


raw_HCC <- read.table("mouse_DB/cancer_expression_DB/hepatocellular_carcinoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_HCC)
colnames(raw_HCC) <- raw_HCC[1,]
raw_HCC <- raw_HCC[-1,]
rowname <- raw_HCC[, 1]
raw_HCC <- raw_HCC[,-1]


raw_HCC <- as.matrix(raw_HCC[, colnames(raw_HCC) %in% c("GSM3467276", "GSM3467277", "GSM3467278", "GSM3537438", "GSM3537439", "GSM3537440")])

rownames(raw_HCC) <- rowname


raw_HCC <- matrix(apply(raw_HCC, 2, as.integer), nrow=nrow(raw_HCC), dimnames=dimnames(raw_HCC))

raw_HCC <- normalized_expression_matrix2(raw_HCC)

new_raw_HCC <- matrix(0, length(rowname), 1)
colnames(new_raw_HCC) <- "HCC"
rownames(new_raw_HCC) <- rowname

for (i in 1:length(rowname)){
  new_raw_HCC[i, 1] <- mean(raw_HCC[i,])
}

saveRDS(raw_HCC, file="/home/xkdl27/mouse_DB/new_normalization_output/HCC.rds")

######## MEL


raw_MEL <- read.table("mouse_DB/cancer_expression_DB/melanoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_MEL)
colnames(raw_MEL) <- raw_MEL[1,]
raw_MEL <- raw_MEL[-1,]
rowname <- raw_MEL[, 1]
raw_MEL <- raw_MEL[,-1]


raw_MEL <- as.matrix(raw_MEL[, colnames(raw_MEL) %in% c("GSM2219662", "GSM2219663", "GSM2219664", "GSM2276599")])

rownames(raw_MEL) <- rowname


raw_MEL <- matrix(apply(raw_MEL, 2, as.integer), nrow=nrow(raw_MEL), dimnames=dimnames(raw_MEL))

raw_MEL <- normalized_expression_matrix2(raw_MEL)

new_raw_MEL <- matrix(0, length(rowname), 1)
colnames(new_raw_MEL) <- "MEL"
rownames(new_raw_MEL) <- rowname

for (i in 1:length(rowname)){
  new_raw_MEL[i, 1] <- mean(raw_MEL[i,])
}
saveRDS(raw_MEL, file="/home/xkdl27/mouse_DB/new_normalization_output/MEL.rds")

######## AML


raw_AML <- read.table("mouse_DB/cancer_expression_DB/Acute_Myeloid_Leukemia_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_AML)
colnames(raw_AML) <- raw_AML[1,]
raw_AML <- raw_AML[-1,]
rowname <- raw_AML[, 1]
raw_AML <- raw_AML[,-1]


raw_AML <- as.matrix(raw_AML[, colnames(raw_AML) %in% c("GSM1672253", "GSM1672254", "GSM1672255", "GSM1672256", "GSM1672257", "GSM1672258", "GSM1672259", "GSM1672260", "GSM1672261", "GSM1672262")])

rownames(raw_AML) <- rowname


raw_AML <- matrix(apply(raw_AML, 2, as.integer), nrow=nrow(raw_AML), dimnames=dimnames(raw_AML))

raw_AML <- normalized_expression_matrix2(raw_AML)

new_raw_AML <- matrix(0, length(rowname), 1)
colnames(new_raw_AML) <- "AML"
rownames(new_raw_AML) <- rowname

for (i in 1:length(rowname)){
  new_raw_AML[i, 1] <- mean(raw_AML[i,])
}
saveRDS(raw_AML, file="/home/xkdl27/mouse_DB/new_normalization_output/AML.rds")

######## NB


raw_NB <- read.table("mouse_DB/cancer_expression_DB/Neuroblastoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_NB)
colnames(raw_NB) <- raw_NB[1,]
raw_NB <- raw_NB[-1,]
rowname <- raw_NB[, 1]
raw_NB <- raw_NB[,-1]


raw_NB <- as.matrix(raw_NB[, colnames(raw_NB) %in% c("GSM3583214", "GSM3583215", "GSM3583216", "GSM3583217", "GSM4159921", "GSM4159922", "GSM4159923", "GSM4819873", "GSM4819874")])

rownames(raw_NB) <- rowname


raw_NB <- matrix(apply(raw_NB, 2, as.integer), nrow=nrow(raw_NB), dimnames=dimnames(raw_NB))

raw_NB <- normalized_expression_matrix2(raw_NB)

new_raw_NB <- matrix(0, length(rowname), 1)
colnames(new_raw_NB) <- "NB"
rownames(new_raw_NB) <- rowname

for (i in 1:length(rowname)){
  new_raw_NB[i, 1] <- mean(raw_NB[i,])
}

saveRDS(raw_NB, file="/home/xkdl27/mouse_DB/new_normalization_output/NB.rds")

######## CCA


raw_CCA <- read.table("mouse_DB/cancer_expression_DB/Cholangiocarcinoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_CCA)
colnames(raw_CCA) <- raw_CCA[1,]
raw_CCA <- raw_CCA[-1,]
rowname <- raw_CCA[, 1]
raw_CCA <- raw_CCA[,-1]


raw_CCA <- as.matrix(raw_CCA)

rownames(raw_CCA) <- rowname


raw_CCA <- matrix(apply(raw_CCA, 2, as.integer), nrow=nrow(raw_CCA), dimnames=dimnames(raw_CCA))

raw_CCA <- normalized_expression_matrix2(raw_CCA)

new_raw_CCA <- matrix(0, length(rowname), 1)
colnames(new_raw_CCA) <- "CCA"
rownames(new_raw_CCA) <- rowname

for (i in 1:length(rowname)){
  new_raw_CCA[i, 1] <- mean(raw_CCA[i,])
}

saveRDS(raw_CCA, file="/home/xkdl27/mouse_DB/new_normalization_output/CCA.rds")

######## ESCC


raw_ESCC <- read.table("/home/xkdl27/mouse_DB/cancer_expression_DB/Esophageal_Squamous_Cell_Carcinoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ESCC)
colnames(raw_ESCC) <- raw_ESCC[1,]
raw_ESCC <- raw_ESCC[-1,]
rowname <- raw_ESCC[, 1]
raw_ESCC <- raw_ESCC[,-1]


raw_ESCC <- as.matrix(raw_ESCC[, colnames(raw_ESCC) %in% c("GSM4272746", "GSM4272747")])

rownames(raw_ESCC) <- rowname


raw_ESCC <- matrix(apply(raw_ESCC, 2, as.integer), nrow=nrow(raw_ESCC), dimnames=dimnames(raw_ESCC))

raw_ESCC <- normalized_expression_matrix2(raw_ESCC)

new_raw_ESCC <- matrix(0, length(rowname), 1)
colnames(new_raw_ESCC) <- "ESCC"
rownames(new_raw_ESCC) <- rowname

for (i in 1:length(rowname)){
  new_raw_ESCC[i, 1] <- mean(raw_ESCC[i,])
}

saveRDS(raw_ESCC, file="/home/xkdl27/mouse_DB/new_normalization_output/ESCC.rds")

######## OS


raw_OS <- read.table("mouse_DB/cancer_expression_DB/Osteosarcoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_OS)
colnames(raw_OS) <- raw_OS[1,]
raw_OS <- raw_OS[-1,]
rowname <- raw_OS[, 1]
raw_OS <- raw_OS[,-1]


raw_OS <- as.matrix(raw_OS[, colnames(raw_OS) %in% c("GSM1013699", "GSM1013700", "GSM1013701", "GSM1013707", "GSM1013708", "GSM1013709", "GSM1013710", "GSM1422264", "GSM1422265", "GSM1422266")])

rownames(raw_OS) <- rowname


raw_OS <- matrix(apply(raw_OS, 2, as.integer), nrow=nrow(raw_OS), dimnames=dimnames(raw_OS))

raw_OS <- normalized_expression_matrix2(raw_OS)

new_raw_OS <- matrix(0, length(rowname), 1)
colnames(new_raw_OS) <- "OS"
rownames(new_raw_OS) <- rowname

for (i in 1:length(rowname)){
  new_raw_OS[i, 1] <- mean(raw_OS[i,])
}

saveRDS(raw_OS, file="/home/xkdl27/mouse_DB/new_normalization_output/OS.rds")

######## RMS


raw_RMS <- read.table("mouse_DB/cancer_expression_DB/Rhabdomyosarcoma_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_RMS)
colnames(raw_RMS) <- raw_RMS[1,]
raw_RMS <- raw_RMS[-1,]
rowname <- raw_RMS[, 1]
raw_RMS <- raw_RMS[,-1]


raw_RMS <- as.matrix(raw_RMS[, colnames(raw_RMS) %in% c("GSM5348396", "GSM5348397")])

rownames(raw_RMS) <- rowname


raw_RMS <- matrix(apply(raw_RMS, 2, as.integer), nrow=nrow(raw_RMS), dimnames=dimnames(raw_RMS))

raw_RMS <- normalized_expression_matrix2(raw_RMS)

new_raw_RMS <- matrix(0, length(rowname), 1)
colnames(new_raw_RMS) <- "RMS"
rownames(new_raw_RMS) <- rowname

for (i in 1:length(rowname)){
  new_raw_RMS[i, 1] <- mean(raw_RMS[i,])
}
saveRDS(raw_RMS, file="/home/xkdl27/mouse_DB/new_normalization_output/RMS.rds")

######## SCLC


raw_SCLC <- read.table("mouse_DB/cancer_expression_DB/Small_Cell_Lung_Cancer_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_SCLC)
colnames(raw_SCLC) <- raw_SCLC[1,]
raw_SCLC <- raw_SCLC[-1,]
rowname <- raw_SCLC[, 1]
raw_SCLC <- raw_SCLC[,-1]


raw_SCLC <- as.matrix(raw_SCLC[, colnames(raw_SCLC) %in% c("GSM2249540", "GSM2249541", "GSM2249543")])

rownames(raw_SCLC) <- rowname


raw_SCLC <- matrix(apply(raw_SCLC, 2, as.integer), nrow=nrow(raw_SCLC), dimnames=dimnames(raw_SCLC))

raw_SCLC <- normalized_expression_matrix2(raw_SCLC)

new_raw_SCLC <- matrix(0, length(rowname), 1)
colnames(new_raw_SCLC) <- "SCLC"
rownames(new_raw_SCLC) <- rowname

for (i in 1:length(rowname)){
  new_raw_SCLC[i, 1] <- mean(raw_SCLC[i,])
}

saveRDS(raw_SCLC, file="/home/xkdl27/mouse_DB/new_normalization_output/SCLC.rds")