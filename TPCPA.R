rm(list=ls())

getwd()
setwd("/Users/naturexv/Datasets/TPCPA_pro/")

library(data.table)
library(Seurat)
library(scCustomize)
library(tidyverse)
library(impute)
library(ggplot2)
library(factoextra)

metadata <- read.table("Metadata.csv", header = T, sep = ",")
colnames(metadata)
metadata <- metadata[,c(1:23)]
metadata <- column_to_rownames(.data = metadata, var = "Sample")

data <- fread("Intensity_data.csv", header = T, sep = ",")
colnames(data)
data <- data[,c(6,11:1182)]
colnames(data)[1] <- "genesymbol"
data <- distinct(.data = data, genesymbol, .keep_all = T)
data <- column_to_rownames(.data = data, var = "genesymbol")

boxplot(data, las = 2, main = "Raw log2 Intensities across samples")

# 假设 data 是 log2 转换后的表达矩阵，行为蛋白，列为样本；先转置，样本为行
data_t <- t(data)

# 找到矩阵中最小的非 NA 值
min_val <- min(data_t, na.rm = TRUE)

# 用最小值 - 1 替换 NA
data_imputed <- data_t
data_imputed[is.na(data_imputed)] <- min_val - 1
data_imputed_t <- t(data_imputed) 

# 创建 Seurat 对象（蛋白表达数据不需要是 raw count）
seu <- CreateSeuratObject(counts = data_imputed_t, meta.data = metadata)

# 将 counts 层的值复制到 data 层
seu[["RNA"]] <- SetAssayData(seu[["RNA"]], slot = "data", new.data = GetAssayData(seu, layer = "counts"))

# 标准化（虽然是 log2 数据，仍推荐进行 ScaleData）
seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)

# PCA 降维
seu <- RunPCA(seu, features = rownames(seu), verbose = FALSE)
ElbowPlot(seu)

# 可视化 PCA（SeuratTheme()# 可视化 PCA（按分组）
DimPlot(seu, reduction = "pca", group.by = "Group") +
  ggtitle("PCA of Proteomics Data")

# 可选：进一步用 UMAP 可视化
seu <- RunUMAP(seu, dims = 1:11, verbose = FALSE)
DimPlot_scCustom(seu, reduction = "umap", group.by = "Group") +
  ggtitle("UMAP of Proteomics Data")


colnames(seu@meta.data)
FeaturePlot_scCustom(seu, features = c('IGSF3'), figure_plot = T)
DimPlot_scCustom(seu, reduction = "umap", group.by = "TCGA.Cancer.Code")
FeaturePlot_scCustom(seu, features = c('IGSF3'), split.by = "Sample.Type")

### 用于差异分析 ###
library(DEP)
impute(data, fun = "MinProb", q = 0.01)  
impute(data, fun = "QRILC")

### 用于聚类或热图 ###
impute.knn()   # KNN
fillna(min_val - 1)   # 最小值减一

table(metadata$Sample.Type)
table(metadata$Group)
table(metadata$Location)

############# 箱线图
library(ggpubr)
data_imputed_t <- as.data.frame(data_imputed_t)
GENE = data_imputed_t["IGSF3",]

z = as.data.frame(t(GENE)) #转置 转化为数据框

z$type = metadata$TCGA.Cancer.Code
z$type <- factor(z$type, levels = c("-", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
                                    "HNSC", "KIRC", "KIRC/KIRP", "KIRP", "LAML", "LIHC", "LIHC/CHOL",
                                    "LUAD", "LUAD/LUSC", "LUSC", "OV", "PAAD", "PRAD", "READ",
                                    "SARC", "SKCM", "STAD", "THCA"))


###
# z$type = metadata$Group
# z$type <- factor(z$type, levels = c("AML", "B-ALL", "Bladder", "Brain", "Breast",
                                      "Cervix", "Colon", "Control", "DLBCL", "Esophagus",
                                      "Gallbladder", "Head-Neck", "Kidney", 
                                      "Leukemia (unassigned)", "Liver", "Lung", 
                                      "Ovary", "Pancreas", "Prostate", "Sarcoma", 
                                      "Skin", "Stomach", "T-ALL", "Thymus", "Thyroid"))

colnames(z)

ggboxplot(z, width = 0.75,
          x = "type", y = "IGSF3", 
          size = 0.5, 
          add = "jitter",
          add.params = list(size = 1, alpha = 0.5)) + #更改点的大小和透明度
  labs(x = "", y = "Relative Expression", color = "") +
  theme_set(theme_classic2()) +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), #调整x轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  theme(legend.position = "none") + 
  ggtitle("IGSF3")







