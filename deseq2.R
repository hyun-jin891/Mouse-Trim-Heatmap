adi <- readRDS("adipose.rds")
adr <- readRDS("adrenal_gland.rds")
amy <- readRDS("amygdala.rds")
bm <- readRDS("bm_matrix.rds")
mam <- readRDS("mammary_gland_matrix.rds")
ce <- readRDS("ce_matrix.rds")
co <- readRDS("co_matrix.rds")
cer <- readRDS("cer_matrix.rds")
cho <- readRDS("choroid_plexus_matrix.rds")
col <- readRDS("colon_matrix.rds")
duo <- readRDS("duodenum_matrix.rds") 
epi <- readRDS("epididymis_matrix.rds")
eso <- readRDS("esophagus_matrix.rds")
fal <- readRDS("fallopian_tube_matrix.rds")
car <- readRDS("cartilage_matrix.rds")
he <- readRDS("heart_matrix.rds")
hip <- readRDS("hippocampus_matrix.rds")
hyp <- readRDS("hypothalamus_matrix.rds")
kid <- readRDS("kidney_matrix.rds")
li <- readRDS("liver_matrix.rds")
lu <- readRDS("lung_matrix.rds")
ly <- readRDS("lymph_node_matrix.rds")
mid <- readRDS("midbrain_matrix.rds")
ova <- readRDS("ovary_matrix.rds") 
pa <- readRDS("pancreas_matrix.rds")
pit <- readRDS("pituitary_gland_matrix.rds")
pl <- readRDS("placenta_matrix.rds") 
pro <- readRDS("prostate_matrix.rds")
rec <- readRDS("rectum_matrix.rds")
ret <- readRDS("retina_matrix.rds") 
sa <- readRDS("salivary_gland_matrix.rds")
sk <- readRDS("skeletal_muscle_matrix.rds") 
ski <- readRDS("skin_matrix.rds")
sm <- readRDS("small_intestine_matrix.rds")
smo <- readRDS("smo_matrix.rds")
spi <- readRDS("spinal_cord_matrix.rds")
spl <- readRDS("spleen_matrix.rds")
sto <- readRDS("stomach_matrix.rds") 
test <- readRDS("testis_matrix.rds")
thy <- readRDS("thymus_matrix.rds")
tra <- readRDS("trachea_matrix.rds") 
to <- readRDS("tongue_matrix.rds")
tendon <- readRDS("patellar_tendon_matrix.rds")
vag <- readRDS("vagina_matrix.rds") 
BDCA <- readRDS("BDCA.rds")
BRCA <- readRDS("BRCA.rds")
COAD <- readRDS("COAD.rds")
GBM <- readRDS("GBM.rds") 
LUAD <- readRDS("LUAD.rds")
PDAC <- readRDS("PDAC.rds")
LUSC <- readRDS("LUSC.rds")
OV <- readRDS("OV.rds") 
prostc <- readRDS("prostate_cancer.rds")
THCA <- readRDS("THCA.rds")
RCC <- readRDS("RCC.rds")
LPS <- readRDS("LPS.rds")
HCC <- readRDS("HCC.rds")
MEL <- readRDS("MEL.rds")
AML <- readRDS("AML.rds") 
NB <- readRDS("NB.rds")
CCA <- readRDS("CCA.rds")
ESCC <- readRDS("ESCC.rds")
OS <- readRDS("OS.rds") 
RMS <- readRDS("RMS.rds")
SCLC <- readRDS("SCLC.rds")

expr_matrix <- cbind(adi, adr, amy, bm, car, ce, cer, cho, co, col, duo, epi, eso, fal,
                     he, hip, hyp, kid, li, lu, ly, mam, mid, ova, pa, tendon, pit, pl, pro, rec, ret, sa, sk, ski,
                     sm, smo, spi, spl, sto, test, thy, to, tra, vag, AML, BDCA, BRCA, CCA, COAD, ESCC, GBM, 
                     HCC, LPS, LUAD, LUSC, MEL, NB, OS, OV, PDAC, prostc, RCC,
                     RMS, SCLC, THCA)


library(DESeq2)

sampleNum <- c(8, 13, 11, 1, 4, 11, 4, 3, 12, 2, 5, 9, 9, 9, 5, 2, 7, 8, 7, 9, 7, 17, 4, 7, 13, 4, 1,
               10, 9, 2, 7, 13, 9, 12, 6, 7, 6, 5, 9, 9, 6, 10, 10, 1, 10, 4, 6, 10, 4, 2, 2, 6, 
               8, 25, 11, 4, 9, 10, 8, 10, 7, 2, 2, 3, 7)


sampleName <- # tissue names



totalSample <- c()
curIndex <- 1


for (i in 1:65){
  for (j in 1:sampleNum[i]){
    totalSample <- c(totalSample, sampleName[i])
  }
  curIndex <- curIndex + sampleNum[i]
}

# metadata (which tissue is current sample ID involved in)
# Generally speaking, which treatment does current sample ID go through

col_data <- data.frame(
  totalSample 
)
rownames(col_data) <- colnames(expr_matrix)
colnames(col_data) <- "condition"


dds <- DESeqDataSetFromMatrix(countData = expr_matrix, colData = col_data, design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
normcount <- counts(dds, normalized=TRUE)

saveRDS(as.data.frame(normcount), file="/home/xkdl27/mouse_DB/heat_map/normalized_matrix.rds")

