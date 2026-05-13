### 下游分析（轨迹、细胞通讯、代谢、bulk）
gene<-unique(table1$Gene)
## 回到单细胞看表达-----------------
scRNA=readRDS('./scRNA_anno.RDS')
load('GCST90435611 table1.Rdata')
scRNA<-scRNA_T
table(scRNA$group1)
Idents(scRNA)<-scRNA$group1
## 关键基因
expr_matrix <- GetAssayData(
  object = scRNA,
  assay = "RNA",  # 保持不变，指定assay
  layer = "data"  # 用layer替代原来的slot，取值与之前一致："counts"、"data"、"scale.data"
)
rownames(expr_matrix)
# 2. 提取关键基因的表达量（确保gene与矩阵行名匹配）
# 先过滤掉不存在于表达矩阵中的基因（避免报错）
gene <- intersect(gene, rownames(expr_matrix))
expr_subset <- as.matrix(expr_matrix[gene, ])  # 此时可转为矩阵
rownames(expr_subset)
annotation_col <- data.frame(Group = Idents(scRNA))
rownames(annotation_col) <- colnames(scRNA)
# 3. （可选）标准化（若用的是data槽，可进一步做z-score；若用scale.data则跳过）
expr_df <- t(scale(t(expr_subset)))  # 按基因标准化
expr_df <- as.data.frame(t(expr_df))
expr_df$Group <- annotation_col$Group
# 后续按分组计算均值（与之前相同）
library(dplyr)
expr_group_mean <- expr_df %>%
  group_by(Group) %>%
  summarise(across(everything(), mean))  # 计算每组基因的平均表达量

# 转换为热图矩阵（行=基因，列=分组）
expr_mean_mat <- as.matrix(t(expr_group_mean[, -1]))
colnames(expr_mean_mat) <- expr_group_mean$Group
# 4. 准备分组注释（此时列是分组，注释信息与列名一致）
annotation_col_mean <- data.frame(
  Group = colnames(expr_mean_mat),  # 分组名
  row.names = colnames(expr_mean_mat)  # 行名=分组名（与热图列名匹配）
)

# 5. 绘制热图
pheatmap(
  mat = expr_mean_mat,  # 分组均值矩阵（核心修改）
  annotation_col = annotation_col_mean,  # 分组注释
  show_rownames = TRUE,  # 显示基因名
  show_colnames = TRUE,  # 显示分组名（此时列是分组，建议显示）
  scale = "none",  # 已标准化，无需再处理
  treeheight_row = 15,  # 基因聚类树高度
  treeheight_col = 15,  # 分组聚类树高度（可选，可关闭）
  main = "关键基因在各组中的平均表达热图",  # 标题修改
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),  # 颜色渐变
  annotation_colors = list(  # 分组颜色（与之前一致）
    Group = c(
      "control" = "#E41A1C",
      "mdd" = "#377EB8",
      "oc" = "#4DAF4A",
      "occontrol" = "#F4A8A8"  # 6位颜色代码
    )
  ),
  fontsize_row = 8,  # 基因名字体大小
  border_color = NA  # 去掉单元格边框
  # 可选：若不想对分组聚类，添加 cluster_cols = FALSE
)

library(Seurat)
library(viridis)
# 创建一个热图，展示特定基因在不同细胞群中的表达情况
# 标准化数据
DotPlot(scRNA,features = gene,cols = c('#dadada','#bc3c29'))

scRNA_T<-scRNAsub
scRNA_T=readRDS('./scRNA_T.RDS')
gene<-"CLSTN3"
# 加载ggplot2包（提供theme()等绘图函数）
library(ggplot2)

# 绘制DotPlot并调整基因名字号
DotPlot(
  object = scRNA,
  features = as.character(gene),
  assay = "RNA"
) + 
  theme(
    axis.text.x = element_text(
      size = 4,  # 建议字号不小于5，否则可能看不清（2太小）
      face = "italic"
    ),
    axis.title.x = element_blank()
  ) +
  RotatedAxis()  # 旋转x轴标签（避免重叠）
dev.off()
library(viridis)
FeaturePlot(scRNA,features = 'CLSTN3',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'))
FeaturePlot(scRNA_T,features = 'CLSTN3',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'))


##  轨迹的基因表达分析----------------------------------

library(mixtools)
#devtools::install_github("SGDDNB/GeneSwitches")
library(GeneSwitches)
#install.packages('fastglm')
library(SingleCellExperiment)


table(scRNA_T$celltype)
## time and expression in CD8_CM
scRNA$celltype
scRNA_em=subset(scRNA,celltype=='CD8_EM')

cellinfo <- scRNA_em@meta.data

## 构建SingleCellExperiment对象
sce <- as.SingleCellExperiment(scRNA_em)

## run
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)

sce_slingshot <- slingshot(sce , clusterLabels = 'celltype', reducedDim = 'UMAP.HARMONY', 
                           start.clus = c(3,5), shrink = 0.2)






dev.off()
## 可视化
cl1 <- cellinfo$celltype
SlingshotDataSet(sce_slingshot)
plot(reducedDims(sce_slingshot)$UMAP.HARMONY,col = brewer.pal(12,"Paired")[cl1],pch=16,asp=1)
## 下面这行关键，否则容易报错！！（特别是Matrx>1.5-0的同学）
igraph::igraph.options(sparsematrices = FALSE)

## 曲线折线仍选一个即可
lines(SlingshotDataSet(sce_slingshot), lwd=2,col = 'black',start.clus = c(3,5))#,type = 'lineages'

legend("right",legend = unique(sce$celltype),
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(3,2,4), pch = 16)

library(slingshot)
library(Seurat)
library(RColorBrewer)

# 假设sce_slingshot是已经完成Slingshot分析的Seurat对


### 开关基因（驱动基因分析）
allexpdata <- as.matrix(scRNA_em@assays$RNA@data);dim(allexpdata)
allcells<-colData(sce_slingshot);dim(allcells)

allcells$slingPseudotime_1

cells <- allcells[!is.na(allcells$slingPseudotime_1),];dim(cells)
expdata <- allexpdata[,rownames(cells)];dim(expdata)

#filter genes expressed in less than 5 cells
#过滤少于五个细胞表达的基因
expdata <- expdata[apply(expdata > 0,1,sum) >= 5,];dim(expdata)

rd_UMAP <- Embeddings(object = scRNA_em, reduction = "umap.harmony");dim(rd_UMAP)#原 object = seu3obj.integrated
rd_UMAP <- rd_UMAP[rownames(cells), ];dim(rd_UMAP)
all(rownames(rd_UMAP) == colnames(expdata))

library(mixtools)
library(GeneSwitches)
library(SingleCellExperiment)


## create SingleCellExperiment object with log-normalized single cell data
## 使用对数规范化的单细胞数据创建SingleCellExperiment对象
library(SingleCellExperiment)

# 假设expdata是表达矩阵
sce <- SingleCellExperiment(assays = List(expdata = expdata))

# 定义不同的阈值
thresholds <- c(0.1, 0.2, 0.3, 0.4)

# 遍历每个阈值，计算二值化数据并绘制直方图
for (threshold in thresholds) {
  binary_data <- ifelse(assays(sce)[["expdata"]] > threshold, 1, 0)
  h <- hist(binary_data, breaks = c(-0.5, 0.5, 1.5), plot = FALSE)
  plot(h, freq = TRUE, xlim = c(-0.5, 1.5), ylim = c(0, max(h$density)), 
       main = paste("Histogram of gene expression (Threshold =", threshold, ")", sep=""),
       xlab = "Gene expression", 
       col = "orange", 
       border = "grey")
  abline(v = threshold, col = "blue")
}



library(SingleCellExperiment)

# 假设expdata是表达矩阵
sce <- SingleCellExperiment(assays = List(expdata = expdata))

# 定义二值化阈值
threshold <- 0.4

# 二值化基因表达数据
binary_data <- ifelse(assays(sce)[["expdata"]] > threshold, 1, 0)
str(binary_data)
# 计算二值化数据的直方图数据
# 将矩阵转换为向量
binary_vector <- as.vector(binary_data)

# 计算二值化数据的直方图数据
# 计算二值化数据的直方图数据
h <- hist(binary_vector, breaks = c(-0.5, 0.5, 1.5), plot = FALSE)
str(h)
# 绘制直方图
# 绘制直方图
barwidth <- diff(h$breaks)[1]  # 计算条形宽度
counts <- h$counts  # 获取计数
mids <- h$mids  # 获取区间中点

# 使用条形图绘制直方图
barplot(counts, names.arg = mids, col = "orange", border = "grey",
        xlim = c(-0.5, 1.5), ylim = c(0, max(counts) * 1.1),  # 调整y轴范围
        main = paste("Histogram of gene expression (Threshold =", threshold, ")", sep=""),
        xlab = "Gene expression", ylab = "Frequency")

# 添加垂直线
abline(v = threshold, col = "blue")
## 二值化分析
sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = 0.2)
# sce_p1 <- binarize_exp(sce_p1, ncores = 3)
sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = TRUE, show_warning = FALSE, zero_ratio = 0.65, ds_cutoff = 0.65)

table(rowData(sce_p1)$prd_quality)

# 过滤出开关基因
sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92);dim(sg_allgenes)


sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92,
                                genelists = gs_genelists);dim(sg_gtypes)#, genetype = c("Surface proteins", "TFs"))

sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),]);dim(sg_vis)

## 自己关注的基因
gl <- unique(table1$Gene)
intersect(sg_vis$geneID, gl)
sg_my <- rowData(sce_p1)[gl,];head(sg_my)
sg_my@rownames
rowData(sce_p1)
sg_my$feature_type <- "Mendelian genes"
sg_vis <- rbind(sg_vis, sg_my)
# 筛选出特定的基因
sg_vis$geneID
sg_vis_filtered <- subset(sg_vis, geneID %in% c("CLSTN3"))
# 绘制时间线图
str(sce_p1)
p <- plot_timeline_ggplot(sg_vis_filtered, timedata = sce_p1$Pseudotime, txtsize = 3.5)
p <- plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3.5)
p
dev.off()
# 首先，找到CLSTN3基因在sg_vis中的位置
library(ggplot2)

# 确保sg_vis数据框中包含pseudotime列
str(sg_vis)
sg_vis

# 首先，找到CLSTN3基因在sg_vis中的位置
clstn3_index <- which(sg_vis$geneID == "CLSTN3")

# 绘制轨迹图
p <- plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3.5)

# 显示图形
print(p)

## R2大于0，上调型开关基因，R2小于0 ，下调型开关基因
a=sce_p1@assays@data$expdata['CLSTN3',] ## !!!
b=sce_p1$Pseudotime

df=data.frame(gene=a,time=b)

ggstatsplot::ggscatterstats(data=df,x='time',y='gene')



####细胞通讯--------阴阳性群
gc()

scRNA=readRDS('./scRNA_anno.RDS')
table(scRNA$T_celltype)
T_celltype <- character(length = length(scRNA$celltype))

# 根据条件分配分组
T_celltype[grepl("CD8_EM", scRNA$celltype)] <- "T cell"
T_celltype[grepl("CD4 T naive", scRNA$celltype)] <- "T cell"
T_celltype[grepl("CD8_CM", scRNA$celltype)] <- "T cell"
T_celltype[grepl("CD4_EM", scRNA$celltype)] <- "T cell"
T_celltype[grepl("Th17", scRNA$celltype)] <- "T cell"
T_celltype[grepl("macro and mono", scRNA$celltype)] <- "macro and mono"
T_celltype[grepl("B cell", scRNA$celltype)] <- "B cell"
T_celltype[grepl("platelet", scRNA$celltype)] <- "platelet"
T_celltype[grepl("DC", scRNA$celltype)] <- "DC"
T_celltype[grepl("NK", scRNA$celltype)] <- "NK"
T_celltype[grepl("mast cell", scRNA$celltype)] <- "mast cell"

# 将分组信息添加到 Seurat 对象的元数据中
scRNA <- AddMetaData(scRNA, metadata = T_celltype, col.name = "T_celltype")
table(scRNA@meta.data$T_celltype)

# 验证结果
table(scRNA@meta.data$T_celltype)

scRNA_T <- subset(scRNA,subset= T_celltype == 'T cell')
cellinfo <- scRNA_T@meta.data

scRNA <- AddMetaData(scRNA, metadata = T_celltype, col.name = "T_celltype")
table(scRNA@meta.data$T_celltype)
scRNA_other=subset(scRNA,T_celltype !='T cell')
scRNA_em=subset(scRNA,celltype == 'CD8_EM')
saveRDS(scRNA_other,"scRNA_other.RDS")
saveRDS(scRNA_EM,"scRNA_EM.RDS")
rm(scRNA)
gc()
scRNA_T=readRDS('./scRNA_T.RDS')


gc()

##!!!
scRNA_em$gene_group=ifelse(scRNA_em@assays$RNA@counts['CLSTN3',]>0,'CLSTN3+CD8_EM','CLSTN3-CD8_EM')
table(scRNA_T$celltype)
scRNA_otherT=subset(scRNA_T,celltype != 'CD8_EM')
table(scRNA_otherT$celltype)
# 加列
scRNA_other$gene_group =scRNA_other$celltype
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_chat=merge(scRNA_em,c(scRNA_other,scRNA_otherT))

rm(scRNA_em,scRNA_otherT)
rm(scRNA_T)
rm(scRNA_other)
gc()
table(scRNA_em$gene_group)
##
table(scRNA_chat$group1)
scRNA_chat_ov=subset(scRNA_chat,group1=='oc')

set.seed(123)
a=sample(1:ncol(scRNA_chat_ov),2000)
scRNA_chat_ov=scRNA_chat_ov[,a]

meta =scRNA_chat_ov@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_ov@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "gene_group")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

##时常deff.off!!!!
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, sources.use = c('CLSTN3+CD8_EM','CLSTN3-CD8_EM'),
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                           sources.use = c('CLSTN3+CD8_EM','CLSTN3-CD8_EM'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble


### 代谢相关------------------

scRNA_T=readRDS('./scRNA_T.RDS')

gc()
scRNA_em=subset(scRNA_T,T_celltype=='CD8_EM')
scRNA_em$gene_group=ifelse(scRNA_em@assays$RNA@counts['CLSTN3',]>0,'CLSTN3+CD8_EM','CLSTN3-CD8_EM')

scRNA_otherT=subset(scRNA_T,T_celltype != 'CD8_EM')

# 加列
table(scRNA_otherT$T_celltype)
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_metab=merge(scRNA_em,c(scRNA_otherT))

gc()

rm(scRNA_chat)
gc()
table(scRNA_metab$group1)
scRNA_metab_OC=subset(scRNA_metab,group1=='oc')

set.seed(123)
a=sample(1:ncol(scRNA_metab_ov),2000)
scRNA_metab_ov=scRNA_metab_ov[,a]




#### scMetabolism评估巨噬细胞代谢活性
#devtools::install_github("YosefLab/VISION")
#devtools::install_github("wu-yc/scMetabolism")
table(scRNA_metab_ov$gene_group)
library(scMetabolism)
library(ggplot2)
library(rsvd)
scRNA_metab_ov<-sc.metabolism.Seurat(obj = scRNA_metab_OC, method = 'AUCell', imputation = F, ncores = 2, metabolism.type = "KEGG")
input.pathway <- rownames(scRNA_metab_ov@assays[["METABOLISM"]][["score"]])[61:90]
DotPlot.metabolism(obj =scRNA_metab_ov,
                   pathway = input.pathway, phenotype = "gene_group", norm = "y")


gc()
##差异基因------------------------------------
library(Seurat)
## Idents()
Idents(scRNA_em)=scRNA_em$gene_group
df=FindAllMarkers(scRNA_em,only.pos = T,logfc.threshold =0.25)
write.csv(df,'CLSTN3_marker.csv',quote = F)

## 富集分析怎么做,负值csv中的基因列到网站，参见下面教程
#https://mp.weixin.qq.com/s/ClHOFvw3GSM9wvmIPip4VA


## bulk---------------------------------------------
gc()

library(data.table)
rt=fread('GSE112087_counts-matrix-EnsembIDs-GRCh37.p10.txt',data.table = F)
rownames(rt)=rt$V1
rt$V1=NULL


#IOBR有一系列依赖包
#https://mp.weixin.qq.com/s/nVziQeInS-4QxVNPQCilVQ
#   devtools::install_github("IOBR/IOBR")
library(IOBR)
#
#gc()
#rm(scRNA_chat_ov)
#rm(scRNA_em)
#rm(scRNA_T)
gc()

rt=count2tpm(rt,idType = 'Ensembl',org = 'hsa')

rt=as.data.frame(rt)

max(rt)
rt=log2(rt+1)

a1=grep('ov',colnames(rt))

exp1=rt[,a1]
exp2=rt[,-a1]

rt=cbind(exp2,exp1)


load('table1.Rdata')

data=rt[unique(table1$Gene),]

### 注释文件
anno=data.frame(row.names =colnames(rt),group=c(rep('Healthy',58),
                                                rep('ov',62)))
pheatmap::pheatmap(data,cluster_cols = F,
                   scale = 'row',show_colnames = F,annotation_col = anno)

df=data.frame(gene=as.numeric(rt['CLSTN3',]),group=anno$group)
ggpubr::ggboxplot(data = df, x = 'group',y='gene',color = 'group',palette = 'jco',notch = T,size = 1)+
  stat_compare_means()

