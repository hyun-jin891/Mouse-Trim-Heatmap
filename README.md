# 🔥Mouse Trim Heatmap
* PDAC(Pancreatic Ductal Adenocarcinoma)에서의 TRIM 계열 유전자 관련 연구를 진행 도중 Mouse의 데이터를 써도 될까? 라는 질문에서 출발한 프로젝트
* Mouse 조직에서의 Trim 계열 유전자 발현 양상을 시각화한 heatmap을 그리기로 결정

* For our study related to PDAC(Pancreatic Ductal Adenocarcinoma), the question "Is it ok to use the data of mouse?" makes me start this project.
* I draw Trim heatmap of mouse

# 💾 DB
* 공개 DB를 사용
* ARCHS4 DB: [Link](https://maayanlab.cloud/archs4/)
* Human과 mouse 특정 조직의 count matrix(Bulk RNA-seq)를 구할 수 있음
* 조직을 검색하면, R script를 다운받을 수 있음
* R script를 실행하면, 로컬에 human 또는 mouse의 전체 DB가 설치됨 (한 번 설치되면, 해당 과정은 생략)
    * 검색했던 조직과 연관된 GEO sample ID를 추출해 count matrix를 구축하고 rds의 형태로 저장

* I use a public DB
* ARCHS4 DB: [Link](https://maayanlab.cloud/archs4/)
* From this DB, I could get the count matrix of the specific tissues of humand and mouse
* If you search the tissue you want, you can download R script
* If you run this script, the humand or mouse gene DB is installed in your environment
    * It automatically extract GEO sample ID related to your focused tissue, make the count matrix, and save this as rds file

![Image](https://github.com/user-attachments/assets/7a4e5e8a-5b2f-4a00-bef9-c7b08e259173)

# 💻 Workflow
* 44개의 mouse 정상 조직과 21개의 암 조직의 count matrix를 ARCHS4 DB를 통해 획득
* 각 조직을 대표할 수 있다고 판단되는 sample ID만을 선별해서, 각 조직별로 전처리 (normal_sample_extraction.R & cancer_sample_extraction.R)
* 모든 조직을 합쳐서 하나의 matrix로 만든 다음, DESeq2를 활용해 normalization 진행 (deseq2.R)
* 조직 당 여러 개 존재하는 sample 열을 median 계산을 통해 하나의 열로 바꿈 (median_deseq2_matrix.R)
* Trim 계열 유전자 행만 추출 -> heatmap 생성 (final_heatmap.R)

* I get the count matrix of 44 normal tissues & 21 cancer tissues through ARCHS4
* I extract some of sample IDs that can be representative of the tissue -> preprocessing per each tissue (normal_sample_extraction.R & cancer_sample_extraction.R)
* I make a single matrix through merging the count matrix of each tissue, and do normalization through DESeq2 (deseq2.R)
* The sample ID columns per tissue are converted to just one column that has median value of sample IDs of one tissue (median_deseq2_matrix.R)
* I extract the Trim gene rows -> draw heatmap (final_heatmap.R)

# ☺️ Project Result

![Image](https://github.com/user-attachments/assets/89f8faa6-5b95-4f4c-bebe-62ce4e41f72f)



# Tissues
* Normal Tissue
    * Brown Adipose
    * Adrenal Gland
    * Amygdala
    * Bone Marrow
    * Mammary Gland
    * Cerebellum
    * Cerebral Cortex
    * Cervix
    * Choroid Plexus
    * Colon
    * Duodenum
    * Epididymis
    * Esophagus
    * Fallopian Tube
    * Cartilage
    * Heart
    * Hippocampus
    * Hypothalamus
    * Kidney
    * Liver
    * Lung
    * Lymph Node
    * Midbrain
    * Ovary
    * Pancreas
    * Pituitary Gland
    * Placenta
    * Prostate
    * Rectum
    * Retina
    * Salivary Gland
    * Skeletal Muscle
    * Skin
    * Small Intestine
    * Smooth Muscle
    * Spleen
    * Spinal Cord
    * Stomach
    * Testis
    * Thymus
    * Trachea
    * Tongue
    * Patellar Tendon
    * Vagina
* Cancer
    * BDCA
    * BRCA
    * COAD
    * GBM
    * LUAD
    * PDAC
    * LUSC
    * OV
    * PCa
    * THCA
    * RCC
    * LPS
    * HCC
    * MEL
    * AML
    * NB
    * CCA
    * ESCC
    * OS
    * RMS
    * SCLC

# Trim Gene
* Trim34a
* Trim6
* Trim5
* Trim21
* Trim68
* Trim58
* Trim38
* Trim17
* Trim11
* Mefv
* Trim27
* Trim7
* Trim10
* Trim26
* Trim39
* Trim41
* Trim60
* Trim61
* Trim40
* Trim31
* Trim31
* Trim43b
* Trim14
* Triml2
* Trim50
* Trim72
* Trim69
* Trim62
* Trim35
* Trim44
* Bspry
* Trim37
* Trim55
* Trim54
* Trim63
* Trim46
* Trim36
* Mid1
* Mid2
* Trim9
* Trim67
* Trim13
* Trim59
* Trim25
* Trim47
* Trim29
* Trim16
* Trim65
* Trim8
* Trim2
* Trim3
* Trim71
* Trim45
* Trim32
* Trim56
* Pml
* Rnf207
* Trim23
* Trim33
* Trim24
* Trim28
* Trim66
* Trim42





