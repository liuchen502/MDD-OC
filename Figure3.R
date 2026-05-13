# 关键细胞亚群的深入分析


## 读取数据
scRNA_T <- LoadH5Seurat("scRNA_T.h5Seurat")
scRNAsub<-LoadH5Seurat("scRNAsub1.h5Seurat")
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)
rm(scRNAsub)
rm(scRNAsub_T)
## 细胞轨迹分析----------------------------

## 载入示例数据
scRNAsub<-scRNA
scRNAsub@reductions
DimPlot(scRNAsub,reduction = 'umap.harmony')
table(scRNAsub$T_celltype)
t_cell_types <- c("CD8_EM", "CD4 T naive", "CD8_CM", "CD4_EM", "Th17","macro and mono","platelet"
                  ,"B cell","DC","NK","mast cell")
# 创建一个新的分组向量
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
table(scRNA_T@meta.data$celltype)

scRNA_T <- subset(scRNAsub,subset= T_celltype == 'T cell')
cellinfo <- scRNA_T@meta.data


## 构建SingleCellExperiment对象
sce <- as.SingleCellExperiment(scRNA_T)

## run
# 查看更新后的降维结果
reducedDimNames(sce)
sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', reducedDim = 'UMAP.HARMONY', 
                           start.clus = c(3,5), shrink = 0.2)



lin1 <- getLineages(sce_slingshot, 
                    clusterLabels = "seurat_clusters", 
                 start.clus ="0",#可指定起始细胞簇，用处不大
                    end.clus="5",#可指定终点细胞簇,用处不大
                    reducedDim = "UMAP.HARMONY")


lin1
dev.off()
## 可视化
cl1 <- cellinfo$T_celltype
plot(reducedDims(sce_slingshot)$UMAP,col = brewer.pal(12,"Paired")[cl1],pch=16,asp=1)

## 下面这行关键，否则容易报错！！（特别是Matrx>1.5-0的同学）

igraph::igraph.options(sparsematrices = FALSE)

## 曲线折线仍选一个即可
#lines(SlingshotDataSet(sce_slingshot), lwd=2,col = 'black')#,type = 'lineages'
lines(SlingshotDataSet(sce_slingshot), lwd=2, type = 'lineages', col = 'black')


legend("right",legend = unique(sce$T_celltype),
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(3,2,4), pch = 16)





# 细胞通讯----------------------
# 在SLE和衰老中分别做，看看功能是不是一样的

## cell-cell chat----
table(scRNAsub$T_celltype)
table(scRNA_T$T_celltype)
scRNA_other <- subset(scRNAsub, subset = T_celltype != 'T cell')
table(scRNA_other$T_celltype)
# 加一列！！！
scRNA_other$T_celltype =scRNA_other$celltype
scRNA_chat=merge(scRNA_other,scRNA_T)
rm(scRNA_T)
rm(scRNA_other)
gc()
table(scRNA_chat$group1)
scRNA_chat_mdd=subset(scRNA_chat,group1=='mdd')


## 抽2000细胞做
set.seed(123)
a=sample(1:ncol(scRNA_chat_mdd),2000)
scRNA_chat_mdd=scRNA_chat_mdd[,a]

meta =scRNA_chat_mdd@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_mdd@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")

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
data("PPI.human")
cellchat <- smoothData(cellchat, adj = PPI.human)
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
                 weight.scale = T, label.edge= T, sources.use = 'macro and mono',
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                           sources.use = 'macro and mono',
                           remove.isolate = FALSE)+coord_flip()
p_bubble

saveRDS(cellchat,"mdd cellchat.rds")

## oc
scRNA_chat_oc=subset(scRNA_chat,group1=='oc')
table(scRNA_chat_oc$celltype)

## 抽2000细胞做
set.seed(123)
a=sample(1:ncol(scRNA_chat_oc),2000)
scRNA_chat_oc=scRNA_chat_oc[,a]

meta =scRNA_chat_oc@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_oc@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")

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
data("PPI.human")
cellchat <- smoothData(cellchat, adj = PPI.human)
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
table(cellchat@meta$celltype)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, sources.use = 'macro and mono',
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                           sources.use = 'macro and mono',
                           remove.isolate = FALSE)+coord_flip()
p_bubble

saveRDS(cellchat,"oc cellchat.rds")

