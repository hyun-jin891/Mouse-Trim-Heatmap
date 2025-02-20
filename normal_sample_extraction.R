######## adipose(Brown)


raw_adi <- read.table("mouse_DB/normal_expression_DB/adipose_tissue_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_adi)
colnames(raw_adi) <- raw_adi[1,]
raw_adi <- raw_adi[-1,]
rowname <- raw_adi[, 1]
raw_adi <- raw_adi[,-1]


raw_adi <- as.matrix(raw_adi[, colnames(raw_adi) %in% c("GSM1605934", "GSM1605935", "GSM1605936", "GSM1605937", "GSM1605938", "GSM1605939", "GSM1605940", "GSM1605941")])

rownames(raw_adi) <- rowname


raw_adi <- matrix(apply(raw_adi, 2, as.integer), nrow=nrow(raw_adi), dimnames=dimnames(raw_adi))

raw_adi <- normalized_expression_matrix(raw_adi)

new_raw_adi <- matrix(0, length(rowname), 1)
colnames(new_raw_adi) <- "brown adipose tissue"
rownames(new_raw_adi) <- rowname

for (i in 1:length(rowname)){
  new_raw_adi[i, 1] <- mean(raw_adi[i,])
}

saveRDS(raw_adi, file="/home/xkdl27/mouse_DB/new_normalization_output/adipose.rds")


######## adrenal gland


raw_adr <- read.table("mouse_DB/normal_expression_DB/adrenal_gland_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_adr)
colnames(raw_adr) <- raw_adr[1,]
raw_adr <- raw_adr[-1,]
rowname <- raw_adr[, 1]
raw_adr <- raw_adr[,-1]


raw_adr <- as.matrix(raw_adr[, colnames(raw_adr) %in% c("GSM1321278", "GSM1321279", "GSM1321281", "GSM1321282", "GSM1321283", "GSM1321284", "GSM1321285", "GSM4662459", "GSM4662460", "GSM4662461", "GSM4662462", "GSM4662463", "GSM4662464")])

rownames(raw_adr) <- rowname


raw_adr <- matrix(apply(raw_adr, 2, as.integer), nrow=nrow(raw_adr), dimnames=dimnames(raw_adr))
raw_adr <- normalized_expression_matrix(raw_adr)
new_raw_adr <- matrix(0, length(rowname), 1)
colnames(new_raw_adr) <- "adrenal gland"
rownames(new_raw_adr) <- rowname

for (i in 1:length(rowname)){
  new_raw_adr[i, 1] <- mean(raw_adr[i,])
}


saveRDS(raw_adr, file="/home/xkdl27/mouse_DB/new_normalization_output/adrenal_gland.rds")


######## amygdala


raw_amy <- read.table("mouse_DB/normal_expression_DB/amygdala_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_amy)
colnames(raw_amy) <- raw_amy[1,]
raw_amy <- raw_amy[-1,]
rowname <- raw_amy[, 1]
raw_amy <- raw_amy[,-1]


raw_amy <- as.matrix(raw_amy[, colnames(raw_amy) %in% c("GSM1860527", "GSM2163594", "GSM2163597", "GSM2163599", "GSM2787214", "GSM4550522", "GSM4550525", "GSM4550527", "GSM4550528", "GSM4550533", "GSM4550534"
)])

rownames(raw_amy) <- rowname


raw_amy <- matrix(apply(raw_amy, 2, as.integer), nrow=nrow(raw_amy), dimnames=dimnames(raw_amy))
raw_amy <- normalized_expression_matrix(raw_amy)
new_raw_amy <- matrix(0, length(rowname), 1)
colnames(new_raw_amy) <- "amygdala"
rownames(new_raw_amy) <- rowname

for (i in 1:length(rowname)){
  new_raw_amy[i, 1] <- mean(raw_amy[i,])
}

saveRDS(raw_amy, file="/home/xkdl27/mouse_DB/new_normalization_output/amygdala.rds")

######## bone marrow

setwd("/home/xkdl27")
raw_bm <- read.table("mouse_DB/normal_expression_DB/bone_marrow_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_bm)
colnames(raw_bm) <- raw_bm[1,]
raw_bm <- raw_bm[-1,]
rowname <- raw_bm[, 1]
raw_bm <- raw_bm[,-1]


raw_bm <- as.matrix(raw_bm[, colnames(raw_bm) %in% c("GSM1006301")])
colnames(raw_bm) <- "GSM1006301"
rownames(raw_bm) <- rowname


raw_bm <- matrix(apply(raw_bm, 2, as.integer), nrow=nrow(raw_bm), dimnames=dimnames(raw_bm))
raw_bm <- normalized_expression_matrix(raw_bm)
new_raw_bm <- matrix(0, length(rowname), 1)
colnames(new_raw_bm) <- "bone marrow"
rownames(new_raw_bm) <- rowname

for (i in 1:length(rowname)){
  new_raw_bm[i, 1] <- mean(raw_bm[i,])
}

saveRDS(raw_bm, file="/home/xkdl27/mouse_DB/new_normalization_output/bm_matrix.rds")

######## mammary gland

setwd("/home/xkdl27")
raw_br <- read.table("mouse_DB/normal_expression_DB/mammary_gland_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_br)
colnames(raw_br) <- raw_br[1,]
raw_br <- raw_br[-1,]
rowname <- raw_br[, 1]
raw_br <- raw_br[,-1]


raw_br <- as.matrix(raw_br[, colnames(raw_br) %in% c("GSM925010", "GSM2510630", "GSM2510631", "GSM2510632", "GSM2510634", "GSM2510635", "GSM2510636","GSM2510637","GSM2510638","GSM2510639","GSM2510640","GSM2510641","GSM2510642","GSM2510643","GSM2510644","GSM2510645","GSM2510646")])

rownames(raw_br) <- rowname


raw_br <- matrix(apply(raw_br, 2, as.integer), nrow=nrow(raw_br), dimnames=dimnames(raw_br))
raw_br <- normalized_expression_matrix(raw_br)
new_raw_br <- matrix(0, length(rowname), 1)
colnames(new_raw_br) <- "mammary gland"
rownames(new_raw_br) <- rowname

for (i in 1:length(rowname)){
  new_raw_br[i, 1] <- mean(raw_br[i,])
}

saveRDS(raw_br, file="/home/xkdl27/mouse_DB/new_normalization_output/mammary_gland_matrix.rds")

######## cerebellum


raw_ce <- read.table("mouse_DB/normal_expression_DB/cerebellum_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ce)
colnames(raw_ce) <- raw_ce[1,]
raw_ce <- raw_ce[-1,]
rowname <- raw_ce[, 1]
raw_ce <- raw_ce[,-1]


raw_ce <- as.matrix(raw_ce[, colnames(raw_ce) %in% c("GSM1063318", "GSM1063320", "GSM1063322", "GSM1321310", "GSM1321311", "GSM1321312", "GSM1321313", "GSM1321314", "GSM1321315", "GSM1321316", "GSM1321317")])

rownames(raw_ce) <- rowname


raw_ce <- matrix(apply(raw_ce, 2, as.integer), nrow=nrow(raw_ce), dimnames=dimnames(raw_ce))
raw_ce <- normalized_expression_matrix(raw_ce)
new_raw_ce <- matrix(0, length(rowname), 1)
colnames(new_raw_ce) <- "cerebellum"
rownames(new_raw_ce) <- rowname

for (i in 1:length(rowname)){
  new_raw_ce[i, 1] <- mean(raw_ce[i,])
}

saveRDS(raw_ce, file='/home/xkdl27/mouse_DB/new_normalization_output/ce_matrix.rds')


######## cerebral cortex


raw_co <- read.table("mouse_DB/normal_expression_DB/cerebral_cortex_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_co)
colnames(raw_co) <- raw_co[1,]
raw_co <- raw_co[-1,]
rowname <- raw_co[, 1]
raw_co <- raw_co[,-1]


raw_co <- as.matrix(raw_co[, colnames(raw_co) %in% c("GSM1974083", "GSM1974084", "GSM1974085", "GSM1974086", "GSM1974087", "GSM1974088", "GSM1974089", "GSM1974090", "GSM1974091", "GSM1974092", "GSM1974093", "GSM1974094")])

rownames(raw_co) <- rowname


raw_co <- matrix(apply(raw_co, 2, as.integer), nrow=nrow(raw_co), dimnames=dimnames(raw_co))
raw_co <- normalized_expression_matrix(raw_co)
new_raw_co <- matrix(0, length(rowname), 1)
colnames(new_raw_co) <- "cerebral cortex"
rownames(new_raw_co) <- rowname

for (i in 1:length(rowname)){
  new_raw_co[i, 1] <- mean(raw_co[i,])
}

saveRDS(raw_co, file='/home/xkdl27/mouse_DB/new_normalization_output/co_matrix.rds')


######## cervix


raw_cer <- read.table("mouse_DB/normal_expression_DB/cervix_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_cer)
colnames(raw_cer) <- raw_cer[1,]
raw_cer <- raw_cer[-1,]
rowname <- raw_cer[, 1]
raw_cer <- raw_cer[,-1]


raw_cer <- as.matrix(raw_cer[, colnames(raw_cer) %in% c("GSM2284953", "GSM2284954", "GSM2598454", "GSM2598455"
)])

rownames(raw_cer) <- rowname


raw_cer <- matrix(apply(raw_cer, 2, as.integer), nrow=nrow(raw_cer), dimnames=dimnames(raw_cer))
raw_cer <- normalized_expression_matrix(raw_cer)
new_raw_cer <- matrix(0, length(rowname), 1)
colnames(new_raw_cer) <- "cervix"
rownames(new_raw_cer) <- rowname

for (i in 1:length(rowname)){
  new_raw_cer[i, 1] <- mean(raw_cer[i,])
}

saveRDS(raw_cer, file='/home/xkdl27/mouse_DB/new_normalization_output/cer_matrix.rds')

######## choroid plexus


raw_ch <- read.table("mouse_DB/normal_expression_DB/choroid_plexus_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ch)
colnames(raw_ch) <- raw_ch[1,]
raw_ch <- raw_ch[-1,]
rowname <- raw_ch[, 1]
raw_ch <- raw_ch[,-1]


raw_ch <- as.matrix(raw_ch[, colnames(raw_ch) %in% c("GSM1533011", "GSM1533012", "GSM1533013")])

rownames(raw_ch) <- rowname


raw_ch <- matrix(apply(raw_ch, 2, as.integer), nrow=nrow(raw_ch), dimnames=dimnames(raw_ch))
raw_ch <- normalized_expression_matrix(raw_ch)
new_raw_ch <- matrix(0, length(rowname), 1)
colnames(new_raw_ch) <- "choroid plexus"
rownames(new_raw_ch) <- rowname

for (i in 1:length(rowname)){
  new_raw_ch[i, 1] <- mean(raw_ch[i,])
}

saveRDS(raw_ch, file='/home/xkdl27/mouse_DB/new_normalization_output/choroid_plexus_matrix.rds')

######## colon


raw_col <- read.table("mouse_DB/normal_expression_DB/colon_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_col)
colnames(raw_col) <- raw_col[1,]
raw_col <- raw_col[-1,]
rowname <- raw_col[, 1]
raw_col <- raw_col[,-1]


raw_col <- as.matrix(raw_col[, colnames(raw_col) %in% c("GSM1020641", "GSM1083971")])

rownames(raw_col) <- rowname


raw_col <- matrix(apply(raw_col, 2, as.integer), nrow=nrow(raw_col), dimnames=dimnames(raw_col))
raw_col <- normalized_expression_matrix(raw_col)
new_raw_col <- matrix(0, length(rowname), 1)
colnames(new_raw_col) <- "colon"
rownames(new_raw_col) <- rowname

for (i in 1:length(rowname)){
  new_raw_col[i, 1] <- mean(raw_col[i,])
}


saveRDS(raw_col, file='/home/xkdl27/mouse_DB/new_normalization_output/colon_matrix.rds')


######## duodenum


raw_du <- read.table("mouse_DB/normal_expression_DB/duodenum_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_du)
colnames(raw_du) <- raw_du[1,]
raw_du <- raw_du[-1,]
rowname <- raw_du[, 1]
raw_du <- raw_du[,-1]


raw_du <- as.matrix(raw_du[, colnames(raw_du) %in% c("GSM1694276", "GSM1694277", "GSM1869947", "GSM1869948", "GSM1869949")])

rownames(raw_du) <- rowname


raw_du <- matrix(apply(raw_du, 2, as.integer), nrow=nrow(raw_du), dimnames=dimnames(raw_du))
raw_du <- normalized_expression_matrix(raw_du)
new_raw_du <- matrix(0, length(rowname), 1)
colnames(new_raw_du) <- "duodenum"
rownames(new_raw_du) <- rowname

for (i in 1:length(rowname)){
  new_raw_du[i, 1] <- mean(raw_du[i,])
}

saveRDS(raw_du, file='/home/xkdl27/mouse_DB/new_normalization_output/duodenum_matrix.rds')

######## epididymis


raw_ep <- read.table("mouse_DB/normal_expression_DB/epididymis_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ep)
colnames(raw_ep) <- raw_ep[1,]
raw_ep <- raw_ep[-1,]
rowname <- raw_ep[, 1]
raw_ep <- raw_ep[,-1]


raw_ep <- as.matrix(raw_ep[, colnames(raw_ep) %in% c("GSM1038116", "GSM4110164", "GSM4110165", "GSM4110166", "GSM4174662", "GSM4174663", "GSM4174664", "GSM4174665", "GSM4174666")])

rownames(raw_ep) <- rowname


raw_ep <- matrix(apply(raw_ep, 2, as.integer), nrow=nrow(raw_ep), dimnames=dimnames(raw_ep))
raw_ep <- normalized_expression_matrix(raw_ep)
new_raw_ep <- matrix(0, length(rowname), 1)
colnames(new_raw_ep) <- "epididymis"
rownames(new_raw_ep) <- rowname

for (i in 1:length(rowname)){
  new_raw_ep[i, 1] <- mean(raw_ep[i,])
}

saveRDS(raw_ep, file='/home/xkdl27/mouse_DB/new_normalization_output/epididymis_matrix.rds')

######## esophagus


raw_es <- read.table("mouse_DB/normal_expression_DB/esophagus_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_es)
colnames(raw_es) <- raw_es[1,]
raw_es <- raw_es[-1,]
rowname <- raw_es[, 1]
raw_es <- raw_es[,-1]


raw_es <- as.matrix(raw_es[, colnames(raw_es) %in% c("GSM1622725", "GSM1622726", "GSM1622727", "GSM1622728", "GSM4484805", "GSM4484806", "GSM4484807", "GSM4484809", "GSM4484810")])

rownames(raw_es) <- rowname


raw_es <- matrix(apply(raw_es, 2, as.integer), nrow=nrow(raw_es), dimnames=dimnames(raw_es))
raw_es <- normalized_expression_matrix(raw_es)
new_raw_es <- matrix(0, length(rowname), 1)
colnames(new_raw_es) <- "esophagus"
rownames(new_raw_es) <- rowname

for (i in 1:length(rowname)){
  new_raw_es[i, 1] <- mean(raw_es[i,])
}

saveRDS(raw_es, file='/home/xkdl27/mouse_DB/new_normalization_output/esophagus_matrix.rds')

######## fallopian tube

setwd("/home/xkdl27")
raw_fa <- read.table("mouse_DB/normal_expression_DB/fallopian_tube_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_fa)
colnames(raw_fa) <- raw_fa[1,]
raw_fa <- raw_fa[-1,]
rowname <- raw_fa[, 1]
raw_fa <- raw_fa[,-1]


raw_fa <- as.matrix(raw_fa[, colnames(raw_fa) %in% c("GSM4555255", "GSM4555256", "GSM4555257", "GSM4555264", "GSM4555265", "GSM4555266", "GSM4555273", "GSM4555274", "GSM4555275" )])

rownames(raw_fa) <- rowname


raw_fa <- matrix(apply(raw_fa, 2, as.integer), nrow=nrow(raw_fa), dimnames=dimnames(raw_fa))
raw_fa <- normalized_expression_matrix(raw_fa)
new_raw_fa <- matrix(0, length(rowname), 1)
colnames(new_raw_fa) <- "fallopian tube"
rownames(new_raw_fa) <- rowname

for (i in 1:length(rowname)){
  new_raw_fa[i, 1] <- mean(raw_fa[i,])
}

saveRDS(raw_fa, file='/home/xkdl27/mouse_DB/new_normalization_output/fallopian_tube_matrix.rds')

######## cartilage


raw_car <- read.table("mouse_DB/normal_expression_DB/cartilage_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_car)
colnames(raw_car) <- raw_car[1,]
raw_car <- raw_car[-1,]
rowname <- raw_car[, 1]
raw_car <- raw_car[,-1]


raw_car <- as.matrix(raw_car[, colnames(raw_car) %in% c("GSM1971620", "GSM2258533", "GSM2258534", "GSM2258535")])

rownames(raw_car) <- rowname


raw_car <- matrix(apply(raw_car, 2, as.integer), nrow=nrow(raw_car), dimnames=dimnames(raw_car))
raw_car <- normalized_expression_matrix(raw_car)
new_raw_car <- matrix(0, length(rowname), 1)
colnames(new_raw_car) <- "cartilage"
rownames(new_raw_car) <- rowname

for (i in 1:length(rowname)){
  new_raw_car[i, 1] <- mean(raw_car[i,])
}

saveRDS(raw_car, file='/home/xkdl27/mouse_DB/new_normalization_output/cartilage_matrix.rds')

######## heart


raw_he <- read.table("mouse_DB/normal_expression_DB/heart_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_he)
colnames(raw_he) <- raw_he[1,]
raw_he <- raw_he[-1,]
rowname <- raw_he[, 1]
raw_he <- raw_he[,-1]


raw_he <- as.matrix(raw_he[, colnames(raw_he) %in% c("GSM1012166", "GSM1015154", "GSM1087136", "GSM1087137", "GSM1087138" )])

rownames(raw_he) <- rowname


raw_he <- matrix(apply(raw_he, 2, as.integer), nrow=nrow(raw_he), dimnames=dimnames(raw_he))
raw_he <- normalized_expression_matrix(raw_he)
new_raw_he <- matrix(0, length(rowname), 1)
colnames(new_raw_he) <- "heart"
rownames(new_raw_he) <- rowname

for (i in 1:length(rowname)){
  new_raw_he[i, 1] <- mean(raw_he[i,])
}

saveRDS(raw_he, file='/home/xkdl27/mouse_DB/new_normalization_output/heart_matrix.rds')

######## hippocampus


raw_hi <- read.table("mouse_DB/normal_expression_DB/hippocampus_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_hi)
colnames(raw_hi) <- raw_hi[1,]
raw_hi <- raw_hi[-1,]
rowname <- raw_hi[, 1]
raw_hi <- raw_hi[,-1]


raw_hi <- as.matrix(raw_hi[, colnames(raw_hi) %in% c("GSM1303800", "GSM1303801" )])

rownames(raw_hi) <- rowname


raw_hi <- matrix(apply(raw_hi, 2, as.integer), nrow=nrow(raw_hi), dimnames=dimnames(raw_hi))
raw_hi <- normalized_expression_matrix(raw_hi)
new_raw_hi <- matrix(0, length(rowname), 1)
colnames(new_raw_hi) <- "hippocampus"
rownames(new_raw_hi) <- rowname

for (i in 1:length(rowname)){
  new_raw_hi[i, 1] <- mean(raw_hi[i,])
}

saveRDS(raw_hi, file='/home/xkdl27/mouse_DB/new_normalization_output/hippocampus_matrix.rds')


######## hypothalamus


raw_hy <- read.table("mouse_DB/normal_expression_DB/hypothalamus_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_hy)
colnames(raw_hy) <- raw_hy[1,]
raw_hy <- raw_hy[-1,]
rowname <- raw_hy[, 1]
raw_hy <- raw_hy[,-1]


raw_hy <- as.matrix(raw_hy[, colnames(raw_hy) %in% c("GSM1321326", "GSM1321327", "GSM1321328", "GSM1321329", "GSM1321330", "GSM1321331", "GSM1321333")])

rownames(raw_hy) <- rowname


raw_hy <- matrix(apply(raw_hy, 2, as.integer), nrow=nrow(raw_hy), dimnames=dimnames(raw_hy))
raw_hy <- normalized_expression_matrix(raw_hy)
new_raw_hy <- matrix(0, length(rowname), 1)
colnames(new_raw_hy) <- "hypothalamus"
rownames(new_raw_hy) <- rowname

for (i in 1:length(rowname)){
  new_raw_hy[i, 1] <- mean(raw_hy[i,])
}

saveRDS(raw_hy, file='/home/xkdl27/mouse_DB/new_normalization_output/hypothalamus_matrix.rds')

######## kidney


raw_ki <- read.table("mouse_DB/normal_expression_DB/kidney_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ki)
colnames(raw_ki) <- raw_ki[1,]
raw_ki <- raw_ki[-1,]
rowname <- raw_ki[, 1]
raw_ki <- raw_ki[,-1]


raw_ki <- as.matrix(raw_ki[, colnames(raw_ki) %in% c("GSM1015153", "GSM1020643", "GSM1198555", "GSM1198556", "GSM1198557", "GSM1198558", "GSM1198559", "GSM1198560")])

rownames(raw_ki) <- rowname


raw_ki <- matrix(apply(raw_ki, 2, as.integer), nrow=nrow(raw_ki), dimnames=dimnames(raw_ki))
raw_ki <- normalized_expression_matrix(raw_ki)
new_raw_ki <- matrix(0, length(rowname), 1)
colnames(new_raw_ki) <- "kidney"
rownames(new_raw_ki) <- rowname

for (i in 1:length(rowname)){
  new_raw_ki[i, 1] <- mean(raw_ki[i,])
}

saveRDS(raw_ki, file='/home/xkdl27/mouse_DB/new_normalization_output/kidney_matrix.rds')

######## liver


raw_li <- read.table("mouse_DB/normal_expression_DB/liver_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_li)
colnames(raw_li) <- raw_li[1,]
raw_li <- raw_li[-1,]
rowname <- raw_li[, 1]
raw_li <- raw_li[,-1]


raw_li <- as.matrix(raw_li[, colnames(raw_li) %in% c("GSM1000571", "GSM1002564", "GSM1015152", "GSM1020644", "GSM1020661", "GSM1101899", "GSM1101900")])

rownames(raw_li) <- rowname


raw_li <- matrix(apply(raw_li, 2, as.integer), nrow=nrow(raw_li), dimnames=dimnames(raw_li))
raw_li <- normalized_expression_matrix(raw_li)
new_raw_li <- matrix(0, length(rowname), 1)
colnames(new_raw_li) <- "liver"
rownames(new_raw_li) <- rowname

for (i in 1:length(rowname)){
  new_raw_li[i, 1] <- mean(raw_li[i,])
}

saveRDS(raw_li, file='/home/xkdl27/mouse_DB/new_normalization_output/liver_matrix.rds')



######## lung


raw_lu <- read.table("mouse_DB/normal_expression_DB/lung_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_lu)
colnames(raw_lu) <- raw_lu[1,]
raw_lu <- raw_lu[-1,]
rowname <- raw_lu[, 1]
raw_lu <- raw_lu[,-1]


raw_lu <- as.matrix(raw_lu[, colnames(raw_lu) %in% c("GSM1020653", "GSM1020662", "GSM1193000", "GSM1193001", "GSM1193935", "GSM1193936", "GSM1193937", "GSM1193938", "GSM1193939")])

rownames(raw_lu) <- rowname


raw_lu <- matrix(apply(raw_lu, 2, as.integer), nrow=nrow(raw_lu), dimnames=dimnames(raw_lu))
raw_lu <- normalized_expression_matrix(raw_lu)
new_raw_lu <- matrix(0, length(rowname), 1)
colnames(new_raw_lu) <- "lung"
rownames(new_raw_lu) <- rowname

for (i in 1:length(rowname)){
  new_raw_lu[i, 1] <- mean(raw_lu[i,])
}
saveRDS(raw_lu, file='/home/xkdl27/mouse_DB/new_normalization_output/lung_matrix.rds')

######## lymph node


raw_ly <- read.table("mouse_DB/normal_expression_DB/lymph_node_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ly)
colnames(raw_ly) <- raw_ly[1,]
raw_ly <- raw_ly[-1,]
rowname <- raw_ly[, 1]
raw_ly <- raw_ly[,-1]


raw_ly <- as.matrix(raw_ly[, colnames(raw_ly) %in% c("GSM1533008", "GSM1533051", "GSM1533052", "GSM1533053", "GSM1533060", "GSM1533061", "GSM1533062")])

rownames(raw_ly) <- rowname


raw_ly <- matrix(apply(raw_ly, 2, as.integer), nrow=nrow(raw_ly), dimnames=dimnames(raw_ly))
raw_ly <- normalized_expression_matrix(raw_ly)
new_raw_ly <- matrix(0, length(rowname), 1)
colnames(new_raw_ly) <- "lymph node"
rownames(new_raw_ly) <- rowname

for (i in 1:length(rowname)){
  new_raw_ly[i, 1] <- mean(raw_ly[i,])
}

saveRDS(raw_ly, file='/home/xkdl27/mouse_DB/new_normalization_output/lymph_node_matrix.rds')

######## midbrain


raw_mi <- read.table("mouse_DB/normal_expression_DB/midbrain_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_mi)
colnames(raw_mi) <- raw_mi[1,]
raw_mi <- raw_mi[-1,]
rowname <- raw_mi[, 1]
raw_mi <- raw_mi[,-1]


raw_mi <- as.matrix(raw_mi[, colnames(raw_mi) %in% c("GSM1573638", "GSM1573639", "GSM1573640", "GSM1573641")])

rownames(raw_mi) <- rowname


raw_mi <- matrix(apply(raw_mi, 2, as.integer), nrow=nrow(raw_mi), dimnames=dimnames(raw_mi))
raw_mi <- normalized_expression_matrix(raw_mi)
new_raw_mi <- matrix(0, length(rowname), 1)
colnames(new_raw_mi) <- "midbrain"
rownames(new_raw_mi) <- rowname

for (i in 1:length(rowname)){
  new_raw_mi[i, 1] <- mean(raw_mi[i,])
}

saveRDS(raw_mi, file='/home/xkdl27/mouse_DB/new_normalization_output/midbrain_matrix.rds')

######## ovary


raw_ov <- read.table("mouse_DB/normal_expression_DB/ovary_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ov)
colnames(raw_ov) <- raw_ov[1,]
raw_ov <- raw_ov[-1,]
rowname <- raw_ov[, 1]
raw_ov <- raw_ov[,-1]


raw_ov <- as.matrix(raw_ov[, colnames(raw_ov) %in% c("GSM1196046", "GSM1196047", "GSM1585050", "GSM1585051", "GSM2154780", "GSM2154781", "GSM2320312")])

rownames(raw_ov) <- rowname


raw_ov <- matrix(apply(raw_ov, 2, as.integer), nrow=nrow(raw_ov), dimnames=dimnames(raw_ov))
raw_ov <- normalized_expression_matrix(raw_ov)
new_raw_ov <- matrix(0, length(rowname), 1)
colnames(new_raw_ov) <- "ovary"
rownames(new_raw_ov) <- rowname

for (i in 1:length(rowname)){
  new_raw_ov[i, 1] <- mean(raw_ov[i,])
}

saveRDS(raw_ov, file='/home/xkdl27/mouse_DB/new_normalization_output/ovary_matrix.rds')

######## pancreas


raw_pa <- read.table("mouse_DB/normal_expression_DB/pancreas_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_pa)
colnames(raw_pa) <- raw_pa[1,]
raw_pa <- raw_pa[-1,]
rowname <- raw_pa[, 1]
raw_pa <- raw_pa[,-1]


raw_pa <- as.matrix(raw_pa[, colnames(raw_pa) %in% c("GSM1150322", "GSM1808836", "GSM1808837", "GSM1808838", "GSM1901172", "GSM1901173", "GSM1901174", "GSM2043596", "GSM2057880", "GSM2057881", "GSM2057882", "GSM2057883", "GSM2057884")])

rownames(raw_pa) <- rowname


raw_pa <- matrix(apply(raw_pa, 2, as.integer), nrow=nrow(raw_pa), dimnames=dimnames(raw_pa))
raw_pa <- normalized_expression_matrix(raw_pa)
new_raw_pa <- matrix(0, length(rowname), 1)
colnames(new_raw_pa) <- "pancreas"
rownames(new_raw_pa) <- rowname

for (i in 1:length(rowname)){
  new_raw_pa[i, 1] <- mean(raw_pa[i,])
}
saveRDS(raw_pa, file='/home/xkdl27/mouse_DB/new_normalization_output/pancreas_matrix.rds')

######## pituitary gland


raw_pi <- read.table("mouse_DB/normal_expression_DB/pituitary_gland_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_pi)
colnames(raw_pi) <- raw_pi[1,]
colname <- colnames(raw_pi)
raw_pi <- raw_pi[-1,]
rowname <- raw_pi[, 1]
raw_pi <- raw_pi[,-1]


raw_pi <- as.matrix(raw_pi)

rownames(raw_pi) <- rowname
colnames(raw_pi) <- colname[2]

raw_pi <- matrix(apply(raw_pi, 2, as.integer), nrow=nrow(raw_pi), dimnames=dimnames(raw_pi))
raw_pi <- normalized_expression_matrix(raw_pi)
new_raw_pi <- matrix(0, length(rowname), 1)
colnames(new_raw_pi) <- "pituitary gland"
rownames(new_raw_pi) <- rowname

for (i in 1:length(rowname)){
  new_raw_pi[i, 1] <- mean(raw_pi[i,])
}
saveRDS(raw_pi, file='/home/xkdl27/mouse_DB/new_normalization_output/pituitary_gland_matrix.rds')

######## placenta

setwd("/home/xkdl27")
raw_pl <- read.table("mouse_DB/normal_expression_DB/placenta_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_pl)
colnames(raw_pl) <- raw_pl[1,]
raw_pl <- raw_pl[-1,]
rowname <- raw_pl[, 1]
raw_pl <- raw_pl[,-1]

raw_pl <- as.matrix(raw_pl[, colnames(raw_pl) %in% c("GSM1088642", "GSM1088643", "GSM1088644", "GSM1088645", "GSM1088646", "GSM1088647", "GSM1088648", "GSM1088649", "GSM1088650", "GSM1088651")])

rownames(raw_pl) <- rowname


raw_pl <- matrix(apply(raw_pl, 2, as.integer), nrow=nrow(raw_pl), dimnames=dimnames(raw_pl))
raw_pl <- normalized_expression_matrix(raw_pl)
new_raw_pl <- matrix(0, length(rowname), 1)
colnames(new_raw_pl) <- "placenta"
rownames(new_raw_pl) <- rowname

for (i in 1:length(rowname)){
  new_raw_pl[i, 1] <- mean(raw_pl[i,])
}

saveRDS(raw_pl, file='/home/xkdl27/mouse_DB/new_normalization_output/placenta_matrix.rds')

######## prostate


raw_pr <- read.table("mouse_DB/normal_expression_DB/prostate_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_pr)
colnames(raw_pr) <- raw_pr[1,]
raw_pr <- raw_pr[-1,]
rowname <- raw_pr[, 1]
raw_pr <- raw_pr[,-1]


raw_pr <- as.matrix(raw_pr[, colnames(raw_pr) %in% c("GSM2150215", "GSM2150216", "GSM2150217", "GSM2152244", "GSM2152245", "GSM2152246", "GSM2152247", "GSM2152248", "GSM2152249")])

rownames(raw_pr) <- rowname


raw_pr <- matrix(apply(raw_pr, 2, as.integer), nrow=nrow(raw_pr), dimnames=dimnames(raw_pr))
raw_pr <- normalized_expression_matrix(raw_pr)
new_raw_pr <- matrix(0, length(rowname), 1)
colnames(new_raw_pr) <- "prostate"
rownames(new_raw_pr) <- rowname

for (i in 1:length(rowname)){
  new_raw_pr[i, 1] <- mean(raw_pr[i,])
}

saveRDS(raw_pr, file='/home/xkdl27/mouse_DB/new_normalization_output/prostate_matrix.rds')

######## rectum


raw_re <- read.table("mouse_DB/normal_expression_DB/rectum_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_re)
colnames(raw_re) <- raw_re[1,]
raw_re <- raw_re[-1,]
rowname <- raw_re[, 1]
raw_re <- raw_re[,-1]


raw_re <- as.matrix(raw_re[, colnames(raw_re) %in% c("GSM2471462", "GSM2471478" )])

rownames(raw_re) <- rowname


raw_re <- matrix(apply(raw_re, 2, as.integer), nrow=nrow(raw_re), dimnames=dimnames(raw_re))
raw_re <- normalized_expression_matrix(raw_re)
new_raw_re <- matrix(0, length(rowname), 1)
colnames(new_raw_re) <- "rectum"
rownames(new_raw_re) <- rowname

for (i in 1:length(rowname)){
  new_raw_re[i, 1] <- mean(raw_re[i,])
}

saveRDS(raw_re, file='/home/xkdl27/mouse_DB/new_normalization_output/rectum_matrix.rds')

######## retina


raw_ret <- read.table("mouse_DB/normal_expression_DB/retina_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ret)
colnames(raw_ret) <- raw_ret[1,]
raw_ret <- raw_ret[-1,]
rowname <- raw_ret[, 1]
raw_ret <- raw_ret[,-1]

print(TRUE %in% (colnames(raw_ret) %in% c("GSM1362143")))
raw_ret <- as.matrix(raw_ret[, colnames(raw_ret) %in% c("GSM1361882", "GSM1362143", "GSM1362144", "GSM1362145", "GSM1482108", "GSM1482109", "GSM1482110")])

rownames(raw_ret) <- rowname


raw_ret <- matrix(apply(raw_ret, 2, as.integer), nrow=nrow(raw_ret), dimnames=dimnames(raw_ret))
raw_ret <- normalized_expression_matrix(raw_ret)
new_raw_ret <- matrix(0, length(rowname), 1)
colnames(new_raw_ret) <- "retina"
rownames(new_raw_ret) <- rowname

for (i in 1:length(rowname)){
  new_raw_ret[i, 1] <- mean(raw_ret[i,])
}

saveRDS(raw_ret, file='/home/xkdl27/mouse_DB/new_normalization_output/retina_matrix.rds')

######## salivary gland


raw_sa <- read.table("mouse_DB/normal_expression_DB/salivary_gland_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_sa)
colnames(raw_sa) <- raw_sa[1,]
raw_sa <- raw_sa[-1,]
rowname <- raw_sa[, 1]
raw_sa <- raw_sa[,-1]


raw_sa <- as.matrix(raw_sa[, colnames(raw_sa) %in% c("GSM2142875", "GSM2142876", "GSM2142877", "GSM2142878", "GSM2142879", "GSM2142880", "GSM2142881", "GSM2142882", "GSM2142883", "GSM2142884", "GSM2142885", "GSM2142886", "GSM2142887")])

rownames(raw_sa) <- rowname


raw_sa <- matrix(apply(raw_sa, 2, as.integer), nrow=nrow(raw_sa), dimnames=dimnames(raw_sa))
raw_sa <- normalized_expression_matrix(raw_sa)
new_raw_sa <- matrix(0, length(rowname), 1)
colnames(new_raw_sa) <- "salivary gland"
rownames(new_raw_sa) <- rowname

for (i in 1:length(rowname)){
  new_raw_sa[i, 1] <- mean(raw_sa[i,])
}

saveRDS(raw_sa, file='/home/xkdl27/mouse_DB/new_normalization_output/salivary_gland_matrix.rds')

######## skeletal muscle


raw_sk <- read.table("mouse_DB/normal_expression_DB/skeletal_muscle_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_sk)
colnames(raw_sk) <- raw_sk[1,]
raw_sk <- raw_sk[-1,]
rowname <- raw_sk[, 1]
raw_sk <- raw_sk[,-1]


raw_sk <- as.matrix(raw_sk[, colnames(raw_sk) %in% c("GSM1015155", "GSM1321358", "GSM1321359", "GSM1321360", "GSM1321361", "GSM1321362", "GSM1321363", "GSM1321364", "GSM1321365")])

rownames(raw_sk) <- rowname


raw_sk <- matrix(apply(raw_sk, 2, as.integer), nrow=nrow(raw_sk), dimnames=dimnames(raw_sk))
raw_sk <- normalized_expression_matrix(raw_sk)
new_raw_sk <- matrix(0, length(rowname), 1)
colnames(new_raw_sk) <- "skeletal muscle"
rownames(new_raw_sk) <- rowname

for (i in 1:length(rowname)){
  new_raw_sk[i, 1] <- mean(raw_sk[i,])
}

saveRDS(raw_sk, file='/home/xkdl27/mouse_DB/new_normalization_output/skeletal_muscle_matrix.rds')




######## skin


raw_ski <- read.table("mouse_DB/normal_expression_DB/skin_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ski)
colnames(raw_ski) <- raw_ski[1,]
raw_ski <- raw_ski[-1,]
rowname <- raw_ski[, 1]
raw_ski <- raw_ski[,-1]


raw_ski <- as.matrix(raw_ski[, colnames(raw_ski) %in% c("GSM1141084", "GSM1141085", "GSM1141086", "GSM1141087", "GSM1141088", "GSM1141089", "GSM1141090", "GSM1141091", "GSM1141092", "GSM1141093", "GSM1141094", "GSM1141095")])

rownames(raw_ski) <- rowname


raw_ski <- matrix(apply(raw_ski, 2, as.integer), nrow=nrow(raw_ski), dimnames=dimnames(raw_ski))
raw_ski <- normalized_expression_matrix(raw_ski)
new_raw_ski <- matrix(0, length(rowname), 1)
colnames(new_raw_ski) <- "skin"
rownames(new_raw_ski) <- rowname

for (i in 1:length(rowname)){
  new_raw_ski[i, 1] <- mean(raw_ski[i,])
}



saveRDS(raw_ski, file='/home/xkdl27/mouse_DB/new_normalization_output/skin_matrix.rds')



######## small intestine


raw_sm <- read.table("mouse_DB/normal_expression_DB/small_intestine_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_sm)
colnames(raw_sm) <- raw_sm[1,]
raw_sm <- raw_sm[-1,]
rowname <- raw_sm[, 1]
raw_sm <- raw_sm[,-1]


raw_sm <- as.matrix(raw_sm[, colnames(raw_sm) %in% c("GSM1971614", "GSM2027190", "GSM2027191", "GSM2027192", "GSM2027193", "GSM2027194")])

rownames(raw_sm) <- rowname


raw_sm <- matrix(apply(raw_sm, 2, as.integer), nrow=nrow(raw_sm), dimnames=dimnames(raw_sm))
raw_sm <- normalized_expression_matrix(raw_sm)
new_raw_sm <- matrix(0, length(rowname), 1)
colnames(new_raw_sm) <- "small intestine"
rownames(new_raw_sm) <- rowname

for (i in 1:length(rowname)){
  new_raw_sm[i, 1] <- mean(raw_sm[i,])
}


saveRDS(raw_sm, file='/home/xkdl27/mouse_DB/new_normalization_output/small_intestine_matrix.rds')



######## smooth muscle


raw_smo <- read.table("mouse_DB/normal_expression_DB/smooth_muscle_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_smo)
colnames(raw_smo) <- raw_smo[1,]
raw_smo <- raw_smo[-1,]
rowname <- raw_smo[, 1]
raw_smo <- raw_smo[,-1]


raw_smo <- as.matrix(raw_smo[, colnames(raw_smo) %in% c("GSM1941383", "GSM1941384", "GSM1941385", "GSM2583647", "GSM3508445", "GSM3508446", "GSM3508447")])

rownames(raw_smo) <- rowname


raw_smo <- matrix(apply(raw_smo, 2, as.integer), nrow=nrow(raw_smo), dimnames=dimnames(raw_smo))
raw_smo <- normalized_expression_matrix(raw_smo)
new_raw_smo <- matrix(0, length(rowname), 1)
colnames(new_raw_smo) <- "smooth muscle"
rownames(new_raw_smo) <- rowname

for (i in 1:length(rowname)){
  new_raw_smo[i, 1] <- mean(raw_smo[i,])
}

saveRDS(raw_smo, file='/home/xkdl27/mouse_DB/new_normalization_output/smo_matrix.rds')

######## spinal cord


raw_sp <- read.table("mouse_DB/normal_expression_DB/spinal_cord_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_sp)
colnames(raw_sp) <- raw_sp[1,]
raw_sp <- raw_sp[-1,]
rowname <- raw_sp[, 1]
raw_sp <- raw_sp[,-1]


raw_sp <- as.matrix(raw_sp[, colnames(raw_sp) %in% c("GSM1040149", "GSM1103369", "GSM1103370", "GSM1346037", "GSM1869164", "GSM1869165")])

rownames(raw_sp) <- rowname


raw_sp <- matrix(apply(raw_sp, 2, as.integer), nrow=nrow(raw_sp), dimnames=dimnames(raw_sp))
raw_sp <- normalized_expression_matrix(raw_sp)
new_raw_sp <- matrix(0, length(rowname), 1)
colnames(new_raw_sp) <- "spinal cord"
rownames(new_raw_sp) <- rowname

for (i in 1:length(rowname)){
  new_raw_sp[i, 1] <- mean(raw_sp[i,])
}
saveRDS(raw_sp, file='/home/xkdl27/mouse_DB/new_normalization_output/spinal_cord_matrix.rds')

######## spleen



raw_spl <- read.table("mouse_DB/normal_expression_DB/spleen_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_spl)
colnames(raw_spl) <- raw_spl[1,]
raw_spl <- raw_spl[-1,]
rowname <- raw_spl[, 1]
raw_spl <- raw_spl[,-1]


raw_spl <- as.matrix(raw_spl[, colnames(raw_spl) %in% c("GSM1020647", "GSM1020664", "GSM1468523", "GSM1468524", "GSM1468525")])

rownames(raw_spl) <- rowname


raw_spl <- matrix(apply(raw_spl, 2, as.integer), nrow=nrow(raw_spl), dimnames=dimnames(raw_spl))
raw_spl <- normalized_expression_matrix(raw_spl)
new_raw_spl <- matrix(0, length(rowname), 1)
colnames(new_raw_spl) <- "spleen"
rownames(new_raw_spl) <- rowname

for (i in 1:length(rowname)){
  new_raw_spl[i, 1] <- mean(raw_spl[i,])
}
saveRDS(raw_spl, file='/home/xkdl27/mouse_DB/new_normalization_output/spleen_matrix.rds')

######## stomach


raw_sto <- read.table("mouse_DB/normal_expression_DB/stomach_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_sto)
colnames(raw_sto) <- raw_sto[1,]
raw_sto <- raw_sto[-1,]
rowname <- raw_sto[, 1]
raw_sto <- raw_sto[,-1]


raw_sto <- as.matrix(raw_sto[, colnames(raw_sto) %in% c("GSM1891196", "GSM2071342", "GSM2071343", "GSM2071376", "GSM2071377", "GSM2071427", "GSM2071428", "GSM2071624", "GSM2071625")])

rownames(raw_sto) <- rowname


raw_sto <- matrix(apply(raw_sto, 2, as.integer), nrow=nrow(raw_sto), dimnames=dimnames(raw_sto))
raw_sto <- normalized_expression_matrix(raw_sto)
new_raw_sto <- matrix(0, length(rowname), 1)
colnames(new_raw_sto) <- "stomach"
rownames(new_raw_sto) <- rowname

for (i in 1:length(rowname)){
  new_raw_sto[i, 1] <- mean(raw_sto[i,])
}

saveRDS(raw_sto, file='/home/xkdl27/mouse_DB/new_normalization_output/stomach_matrix.rds')

######## testis


raw_te <- read.table("mouse_DB/normal_expression_DB/testis_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_te)
colnames(raw_te) <- raw_te[1,]
raw_te <- raw_te[-1,]
rowname <- raw_te[, 1]
raw_te <- raw_te[,-1]


raw_te <- as.matrix(raw_te[, colnames(raw_te) %in% c("GSM1029433", "GSM1029434", "GSM1029435", "GSM1029436", "GSM1029437", "GSM1029438", "GSM1029439", "GSM1229972", "GSM1234254")])

rownames(raw_te) <- rowname


raw_te <- matrix(apply(raw_te, 2, as.integer), nrow=nrow(raw_te), dimnames=dimnames(raw_te))
raw_te <- normalized_expression_matrix(raw_te)
new_raw_te <- matrix(0, length(rowname), 1)
colnames(new_raw_te) <- "testis"
rownames(new_raw_te) <- rowname

for (i in 1:length(rowname)){
  new_raw_te[i, 1] <- mean(raw_te[i,])
}
saveRDS(raw_te, file='/home/xkdl27/mouse_DB/new_normalization_output/testis_matrix.rds')

######## thymus


raw_tm <- read.table("mouse_DB/normal_expression_DB/thymus_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_tm)
colnames(raw_tm) <- raw_tm[1,]
raw_tm <- raw_tm[-1,]
rowname <- raw_tm[, 1]
raw_tm <- raw_tm[,-1]


raw_tm <- as.matrix(raw_tm[, colnames(raw_tm) %in% c("GSM1533034", "GSM1533035", "GSM1533036", "GSM1533066", "GSM1533067", "GSM1533068")])

rownames(raw_tm) <- rowname


raw_tm <- matrix(apply(raw_tm, 2, as.integer), nrow=nrow(raw_tm), dimnames=dimnames(raw_tm))
raw_tm <- normalized_expression_matrix(raw_tm)
new_raw_tm <- matrix(0, length(rowname), 1)
colnames(new_raw_tm) <- "thymus"
rownames(new_raw_tm) <- rowname

for (i in 1:length(rowname)){
  new_raw_tm[i, 1] <- mean(raw_tm[i,])
}
saveRDS(raw_tm, file='/home/xkdl27/mouse_DB/new_normalization_output/thymus_matrix.rds')

######## trachea


raw_th <- read.table("mouse_DB/normal_expression_DB/Trachea_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_th)
colnames(raw_th) <- raw_th[1,]
raw_th <- raw_th[-1,]
rowname <- raw_th[, 1]
raw_th <- raw_th[,-1]


raw_th <- as.matrix(raw_th[, colnames(raw_th) %in% c("GSM4258341", "GSM4258342", "GSM4258343", "GSM4258344", "GSM4258345", "GSM4258346", "GSM4258347", "GSM4258348", "GSM4258349", "GSM4258350")])

rownames(raw_th) <- rowname


raw_th <- matrix(apply(raw_th, 2, as.integer), nrow=nrow(raw_th), dimnames=dimnames(raw_th))
raw_th <- normalized_expression_matrix(raw_th)
new_raw_th <- matrix(0, length(rowname), 1)
colnames(new_raw_th) <- "trachea"
rownames(new_raw_th) <- rowname

for (i in 1:length(rowname)){
  new_raw_th[i, 1] <- mean(raw_th[i,])
}
saveRDS(raw_th, file='/home/xkdl27/mouse_DB/new_normalization_output/trachea_matrix.rds')

######## tongue


raw_to <- read.table("mouse_DB/normal_expression_DB/tongue_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_to)
colnames(raw_to) <- raw_to[1,]
raw_to <- raw_to[-1,]
rowname <- raw_to[, 1]
raw_to <- raw_to[,-1]


raw_to <- as.matrix(raw_to[, colnames(raw_to) %in% c("GSM1119114", "GSM1310899", "GSM1310900", "GSM1310901", "GSM1310902", "GSM1310903", "GSM1970851", "GSM1970852", "GSM1970853", "GSM1970854")])

rownames(raw_to) <- rowname


raw_to <- matrix(apply(raw_to, 2, as.integer), nrow=nrow(raw_to), dimnames=dimnames(raw_to))
raw_to <- normalized_expression_matrix(raw_to)
new_raw_to <- matrix(0, length(rowname), 1)
colnames(new_raw_to) <- "tongue"
rownames(new_raw_to) <- rowname

for (i in 1:length(rowname)){
  new_raw_to[i, 1] <- mean(raw_to[i,])
}

saveRDS(raw_to, file='/home/xkdl27/mouse_DB/new_normalization_output/tongue_matrix.rds')


######## patellar tendon


raw_ur <- read.table("/home/xkdl27/mouse_DB/normal_expression_DB/tendon_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_ur)
colnames(raw_ur) <- raw_ur[1,]
raw_ur <- raw_ur[-1,]
rowname <- raw_ur[, 1]
raw_ur <- raw_ur[,-1]


raw_ur <- as.matrix(raw_ur[, colnames(raw_ur) %in% c("GSM4112693", "GSM4112694", "GSM4112695", "GSM4112696")])

rownames(raw_ur) <- rowname


raw_ur <- matrix(apply(raw_ur, 2, as.integer), nrow=nrow(raw_ur), dimnames=dimnames(raw_ur))
raw_ur <- normalized_expression_matrix(raw_ur)
new_raw_ur <- matrix(0, length(rowname), 1)
colnames(new_raw_ur) <- "patellar tendon"
rownames(new_raw_ur) <- rowname

for (i in 1:length(rowname)){
  new_raw_ur[i, 1] <- mean(raw_ur[i,])
}

saveRDS(raw_ur, file='/home/xkdl27/mouse_DB/new_normalization_output/patellar_tendon_matrix.rds')


######## vagina


raw_v <- read.table("mouse_DB/normal_expression_DB/vagina_expression_matrix.tsv", sep="\t", header=FALSE)
getwd()

print(raw_v)
colnames(raw_v) <- raw_v[1,]
colname <- colnames(raw_v)
raw_v <- raw_v[-1,]
rowname <- raw_v[, 1]
raw_v <- raw_v[,-1]


raw_v <- as.matrix(raw_v[, colnames(raw_v) %in% c("GSM1654317")])

rownames(raw_v) <- rowname
colnames(raw_v) <- colname[2]


raw_v <- matrix(apply(raw_v, 2, as.integer), nrow=nrow(raw_v), dimnames=dimnames(raw_v))
raw_v <- normalized_expression_matrix(raw_v)
new_raw_v <- matrix(0, length(rowname), 1)
colnames(new_raw_v) <- "vagina"
rownames(new_raw_v) <- rowname

for (i in 1:length(rowname)){
  new_raw_v[i, 1] <- mean(raw_v[i,])
}

saveRDS(raw_v, file='/home/xkdl27/mouse_DB/new_normalization_output/vagina_matrix.rds')
