getwd()
setwd("path/to/your/ukb/data")
library(Seurat)
data1 <- seurat_control
data2 <- seurat_mdd
data3 <- seurat_oc
data4<- seurat_occontrol
### 把数据变成list
sceList <- list(data1,data2,data3,data4)
mat <- c("control","mdd","oc","occontrol")
################################################
### 2.Merge 合并数据
combined = merge(sceList[[1]],y = sceList[-1],add.cell.ids = mat)
orig <- NULL
for(i in 1:length(sceList)){
  orig = append(orig,rep(mat[i],ncol(sceList[[i]])))
}
combined[["orig.ident"]] = orig
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

### 该操作会在metadata数据里面增加一列叫做percent.mt
metadata <- combined@meta.data

### 质控数据可视化，使用VlnPlot函数
### nFeature_RNA, number of Feature, 每个细胞中有多少个基因
### nCount_RNA, number of counts, 每个细胞中有多少个counts
### percent.mt, 我们自己增加的列,  每个细胞中线粒体基因的比例
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4)

### 正式筛选，筛选的是细胞，最终细胞减少
### nFeature_RNA > 200
### nFeature_RNA < 2500
### percent.mt < 5
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, features = rownames(combined))
combined <- RunPCA(combined, features = VariableFeatures(object = combined),reduction.name = "pca")
combined <- FindNeighbors(combined, dims = 1:30, reduction = "pca")
combined <- IntegrateLayers(
  object = combined, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
### 设置多个resolution选择合适的resolution
combined <- FindClusters(combined, resolution = seq(0.2,1.2,0.1))
### 选定分辨率0.5作为分群
library(clustree)
clustree(combined)
combined@meta.data$seurat_clusters <- combined@meta.data$integrated_snn_res.0.5
combined@meta.data$seurat_clusters
Idents(combined) <- "seurat_clusters"
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p3 <- DimPlot(combined,reduction = "umap.harmony",group.by = c("orig.ident", "harmony_clusters"),combine = FALSE,label.size = 2)
wrap_plots(c(p3), ncol = 3, byrow = F)
names(combined@assays$RNA@counts)
library(hdf5r)
library(SeuratDisk)
getwd()
setwd("/home/xing54321/RAW_count_matrix_per_sample/MDD")
combined <- JoinLayers(combined)
library(Seurat)
library(SeuratDisk)
# 假设 combined 是你的 Seurat 对象
# 使用 layer 参数获取 counts 数据
counts_data <- GetAssayData(combined[["RNA"]], layer = "counts")
# 创建新的 Assay 对象
combined[["RNA"]] <- CreateAssayObject(counts = counts_data)
# 保存为 h5Seurat 格式
SaveH5Seurat(combined, filename = "combined1.h5Seurat", overwrite = TRUE)
combined <- LoadH5Seurat("combined1.h5Seurat")
### 多个分辨率的分群信息会保存在metadata中
metadata <- combined@meta.data
library(clustree)
clustree(combined)
### 选定分辨率0.5作为分群
combined@meta.data$seurat_clusters <- combined@meta.data$RNA_snn_res.0.5
table(combined@meta.data$seurat_clusters)
Idents(combined) <- "seurat_clusters"

### 作图展示
DimPlot(combined, reduction = "umap.harmony", label = T)
### 学习Dimplot的两个重要参数 group.by 和split.by
### group.by, 在一张图中展示多个信息
DimPlot(combined, reduction = "umap.harmony", group.by = "group1")
### split.by, 分成多个信息来展示
DimPlot(combined, reduction = "umap.harmony", split.by = "group1")

### 保存数据
## 细胞比例图
library(reshape2)
library(ggplot2)
library(dplyr)
scRNAsub<-combined
prop_df <- table(scRNAsub@meta.data$seurat_clusters,scRNAsub@meta.data$group1) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

# 作图
# 生成足够的颜色
num_clusters <- length(unique(prop_df$Cluster))
col_vector <- hcl(seq(0, 360, length.out = num_clusters), 100, 65)

# 作图
prop <- ggplot(data = prop_df, aes(x = Number, y = Sample, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = col_vector) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Ratio") +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))

print(prop)
table(combined$group.patient)
# T细胞亚群注释参考
# https://blog.csdn.net/weixin_52505487/article/details/126687526
Idents(scRNAsub)=scRNAsub$seurat_clusters
T_marker <- c("CCR7","LEF1", "TCF7",'SELL','KLF2', #CD4_Naive
              "ANXA1", "CXCR4", "IL2", #CD4_EM
              "BCL6", "CXCR5","ICA1", #CD4_FH
              'IL23R',"CCR6",'CAPG','RORC','IL17A', #TH17
              'FOXP3','IL2RA','IL1R2',#CD4_REG
              'CD8A','GZMK',  # CD8_EM
              'GZMA','CCL5',  #CD8_CM
              'HAVCR2','PDCD1','LAG3', # CD8_exhau
              'EPCAM','CD19','CD3E') 

genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'FCGR3A',
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2', ## fibo 
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                   'MKI67' , 'TOP2A', 
                   'PECAM1', 'VWF',  ## endo 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )

marker_genes <- c("MS4A1", "CD79A","CD19") ### B
marker_genes <- c("GNLY", "NKG7") ### NK
marker_genes <- c("CD3E","CD8A","CD4","IL7R") ### T
marker_genes <- c("CD14", "FCGR3A", "LYZ") ### Monocyte
marker_genes <- c("FCER1A", "PPBP") ### DC
marker_genes <- c("PPBP") ###血小板
marker_genes <- c("CLEC9A","ITGAM","ITGAE","FCER1A") ## cDC
marker_genes <- c("IL3RA","HLA-DRA") ## pDC
marker_genes <- c('FOXP3','TNFRSF9','CTLA4') ## treg
marker_genes <- c('CD11c','CD123','CD303') ## 树突状


## Plasma Cells
DotPlot(scRNAsub, 
        features = marker_genes,
        group.by = "seurat_clusters") + coord_flip()
scRNAsub.markers <- FindAllMarkers(object = scRNAsub, only.pos = FALSE)
write.table(scRNAsub.markers, file = "cluster markers.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
library(dplyr)
top3 <- scRNAsub.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E',   # T cells
                   'KLRD1', # natural killer (NK) cells
                   'FCGR3B', #neutrophils
                   'MS4A1' ,'IGHG4', # b cells 
                   'RPE65','RCVRN',
                   'CD68', 'CD163', 'CD14',   'MKI67' ,'TOP2A',
                   'MLANA', 'MITF',  'DCT','PRAME' , 'GEP' )
library(stringr)  
p_paper_markers <- DotPlot(scRNAsub, features = genes_to_check,
                           assay='RNA'  )  + coord_flip()

p_paper_markers
ggsave(plot=p_paper_markers,
       filename="check_paper_marker_by_seurat_cluster.pdf",width = 12)



# T Cells (CD3D, CD3E, CD8A), 
# B cells (CD19, CD79A, MS4A1 [CD20]), 
# Plasma cells (IGHG1, MZB1, SDC1, CD79A), 
# Monocytes and macrophages (CD68, CD163, CD14),
# NK Cells (FGFBP2, FCG3RA, CX3CR1),  
# Photoreceptor cells (RCVRN), 
# Fibroblasts (FGF7, MME), 
# Endothelial cells (PECAM1, VWF). 
# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
#   immune (CD45+,PTPRC), epithelial/cancer (EpCAM+,EPCAM), 
# stromal (CD10+,MME,fibo or CD31+,PECAM1,endo) 

library(ggplot2) 
genes_to_check = c('SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'FGF7','MME', 'ACTA2',
                   'PECAM1', 'VWF', 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )
library(stringr)  
p_all_markers <- DotPlot(scRNAsub, features = genes_to_check,
                         assay='RNA'  )  + coord_flip()

p_all_markers
ggsave(plot=p_all_markers, filename="check_all_marker_by_seurat_cluster.pdf")



genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CCR7', 'SELL' , 'TCF7','CXCR6' , 'ITGA1',
                   'FOXP3', 'IL2RA',  'CTLA4','GZMB', 'GZMK','CCL5',
                   'IFNG', 'CCL4', 'CCL3' ,
                   'PRF1' , 'NKG7') 
library(stringr)  
p <- DotPlot(scRNAsub, features = genes_to_check,
             assay='RNA'  )  + coord_flip()

p
ggsave(plot=p, filename="check_Tcells_marker_by_seurat_cluster.pdf")

# mast cells, TPSAB1 and TPSB2 
# B cell,  CD79A  and MS4A1 (CD20) 
# naive B cells, such as MS4A1 (CD20), CD19, CD22, TCL1A, and CD83, 
# plasma B cells, such as CD38, TNFRSF17 (BCMA), and IGHG1/IGHG4
genes_to_check = c('CD3D','MS4A1','CD79A',
                   'CD19', 'CD22', 'TCL1A',  'CD83', #  naive B cells
                   'CD38','TNFRSF17','IGHG1','IGHG4', # plasma B cells,
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'PTPRC' ) 
p <- DotPlot(scRNAsub, features = genes_to_check,
             assay='RNA'  )  + coord_flip()

p
ggsave(plot=p, filename="check_Bcells_marker_by_seurat_cluster.pdf")


genes_to_check = c('CD68', 'CD163', 'CD14',  'CD86', 'LAMP3', ## DC 
                   'CD68',  'CD163','MRC1','MSR1','ITGAE','ITGAM','ITGAX','SIGLEC7', 
                   'MAF','APOE','FOLR2','RELB','BST2','BATF3')
p <- DotPlot(scRNAsub, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()

p
ggsave(plot=p, filename="check_myeloids_marker_by_seurat_cluster.pdf")


# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
# - alveolar type I cell (AT1; AGER+)
# - alveolar type II cell (AT2; SFTPA1)
# - secretory club cell (Club; SCGB1A1+)
# - basal airway epithelial cells (Basal; KRT17+)
# - ciliated airway epithelial cells (Ciliated; TPPP3+) 

genes_to_check = c(  'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' ,
                     'AGER','SFTPA1','SCGB1A1','KRT17','TPPP3',
                     'KRT4','KRT14','KRT8','KRT18',
                     'CD3D','PTPRC' ) 
p <- DotPlot(scRNAsub, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()

p
ggsave(plot=p, filename="check_epi_marker_by_seurat_cluster.pdf")


genes_to_check = c('TEK',"PTPRC","EPCAM","PDPN","PECAM1",'PDGFRB',
                   'CSPG4','GJB2', 'RGS5','ITGA7',
                   'ACTA2','RBP1','CD36', 'ADGRE5','COL11A1','FGF7', 'MME')
p <- DotPlot(scRNAsub, features = unique(genes_to_check),
             assay='RNA'  )  + coord_flip()

p
ggsave(plot=p, filename="check_stromal_marker_by_seurat_cluster.pdf")

scRNAsub <- RenameIdents(scRNAsub,
                         "0"="CD8_EM",
                         "1"="CD4 T naive", 
                         "2"="CD8_CM", 
                         "3"= "CD4_EM", 
                         "4"= "macro and mono", 
                         "5"= "CD8_EM",
                         "6"= "B cell", 
                         "7"= "macro and mono", 
                         "8"= "platelet",
                         "9"= "CD8_CM", 
                         "10"= "DC", 
                         "11"= "CD4_EM",
                         "12"= "CD4 T naive", 
                         "13"= "CD8_EM",
                         "14"= "CD8_EM",
                         "15"= "Th17",
                         "16"= "B cell",
                         "17"= "macro and mono",
                         "18"= "B cell" ,
                         "19"= "NK", 
                         "20"= "CD4 T naive",
                         "21"= "mast cell",
                         "22"= "CD8_EM",
                         "23"= "macro and mono")

scRNAsub@meta.data$celltype = Idents(scRNAsub)

Idents(scRNAsub) <- scRNAsub@meta.data$seurat_clusters
scRNAsub@meta.data$seurat_clusters <- Idents(scRNAsub)


#设置idents主要识别标签
Idents(scRNAsub)=scRNAsub@meta.data$celltype
#999999（灰色）
#b38f8f（浅棕色）
#a67b5b（深棕色）
#c2a2a2（淡粉色）
#a8a8a8（深灰色）
#b39d9d（浅灰色）
#a6a6a6（中灰色）
#c2b3b3（淡灰色）
#a69d9d（深灰色）
#b3a6a6（浅灰色）

# 加载必要的库
library(Seurat)
library(RColorBrewer)

# 假设 scRNAsub 是你的 Seurat 对象
scRNAsub@meta.data$seurat_clusters <- Idents(scRNAsub)

# 检查 celltype 的唯一值数量
unique_celltypes <- unique(scRNAsub@meta.data$celltype)
num_celltypes <- length(unique_celltypes)
print(num_celltypes)

# 生成足够的颜色
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[1:num_celltypes]  # 确保颜色数量与细胞类型数量一致

# 绘制聚类图
DimPlot(scRNAsub, group.by = "celltype", label = TRUE, label.size = 3, cols = col_vector, pt.size = 0.5)
DimPlot(scRNAsub, group.by = "celltype",reduction = "umap.harmony", label = T)

## 细胞比例图再次

library(reshape2)
library(ggplot2)
prop_df <- table(scRNAsub@meta.data$celltype,scRNAsub@meta.data$group1) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

sample_color <- col_vector[1:10] 

prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
prop

library(Seurat)
library(SeuratDisk)
# 假设 scRNAsub 是你的 Seurat 对象
# 使用 layer 参数获取 counts 数据
counts_data <- GetAssayData(scRNAsub[["RNA"]], layer = "counts")
# 创建新的 Assay 对象
scRNAsub[["RNA"]] <- CreateAssayObject(counts = counts_data)
# 保存为 h5Seurat 格式
SaveH5Seurat(scRNAsub, filename = "scRNAsub1.h5Seurat", overwrite = TRUE)
scRNAsub1 <- LoadH5Seurat("scRNAsub1.h5Seurat")


