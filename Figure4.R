## 关键基因的孟德尔随机化分析
library(dplyr)
# 读取数据
scRNA=LoadH5Seurat("scRNAsub1.h5Seurat")
scRNAsub=LoadH5Seurat("scRNA_T.h5Seurat")
scRNAsub<-scRNA_T
library(Seurat)
table(scRNAsub$T_celltype)
scRNA_other=subset(scRNA,T_celltype !='T cell')
scRNA_EM=subset(scRNAsub,celltype == 'CD8_EM')

##！！ 加一列
gc()
scRNA_other$celltype=scRNA_other$T_celltype
scRNA_compare=merge(scRNA_other,scRNA_EM)
table(scRNA_compare$celltype)

gc()
rm()


## 1. EM相对其他T亚群
df_EM=FindMarkers(scRNAsub,ident.1 = 'CD8_EM',only.pos = T,logfc.threshold = 1)

## 2. EM相对其他细胞
df_T=FindMarkers(scRNA_compare,ident.1 = 'CD8_EM',only.pos = T,logfc.threshold = 1)

ss=intersect(rownames(df_EM),rownames(df_T))

ss
saveRDS(scRNA,"scRNA.rds")
saveRDS(scRNAsub,"scRNA_T.rds")
saveRDS(scRNA_other,"scRNA_other.rds")
## 保存，后面要用
save(ss,file ='key_marker_gene1.Rdata')

load('key_marker_gene1.Rdata')
write.csv(ss,file ='key_marker_gene.csv',quote = F)

## 自助教程：
## marker基因的富集分析怎么做,负值csv中的基因列到网站，参见下面教程
#https://mp.weixin.qq.com/s/ClHOFvw3GSM9wvmIPip4VA

gene=ss
gene=gene[!duplicated(gene)]

## 转ENSEMBL
library(clusterProfiler)
library(org.Hs.eg.db)
df=bitr(geneID = gene,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
id=df$ENSEMBL
id=id[!duplicated(id)]
id

## 载入包
devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#https://gwas.mrcieu.ac.uk/datasets/
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyMjM2NjI5OTQyQHFxLmNvbSIsImlhdCI6MTc2MTI5NTE1MCwiZXhwIjoxNzYyNTA0NzUwfQ.uhtezOVDFu-yvlxgPdQSdAbikB_u2cOvNXifcbE9khuBw7UuZZ6InR_d9fB91ZwoOWgWRGAKhORuv8iRnElgp-YZ718gsDoAOU-4z6ZljDoTYRpKEekoJX3k8_hW8xpyY8K7RIzVJO-hrzZqFK8a2i6oAVWkC3ehNEOy_LHFsouOMiKxs7vCsBaM4rNQ5tjrVExxEsiVUJefH-PL_tbrxfZa_NXtJMBJC-skbD1CosRXIDGX7163l5ZN_Hlvz3Rtt8bUe2iZslLElR06Zx-D46T4EuAK6fuCOKBCU74Hj6oOXv1J8aM2h8kTK2L_A5Pl0TIUCB3LMoJhxx89-JTIKQ") #这个token就是保存下来一长串的值，复制进去就可以了
OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyMjM2NjI5OTQyQHFxLmNvbSIsImlhdCI6MTc2MTI5NTE1MCwiZXhwIjoxNzYyNTA0NzUwfQ.uhtezOVDFu-yvlxgPdQSdAbikB_u2cOvNXifcbE9khuBw7UuZZ6InR_d9fB91ZwoOWgWRGAKhORuv8iRnElgp-YZ718gsDoAOU-4z6ZljDoTYRpKEekoJX3k8_hW8xpyY8K7RIzVJO-hrzZqFK8a2i6oAVWkC3ehNEOy_LHFsouOMiKxs7vCsBaM4rNQ5tjrVExxEsiVUJefH-PL_tbrxfZa_NXtJMBJC-skbD1CosRXIDGX7163l5ZN_Hlvz3Rtt8bUe2iZslLElR06Zx-D46T4EuAK6fuCOKBCU74Hj6oOXv1J8aM2h8kTK2L_A5Pl0TIUCB3LMoJhxx89-JTIKQ"
# 暴露数据的SNP，工具变量
exposure_id=paste0('eqtl-a-',id)
exposure_dat <- extract_instruments(exposure_id, p1=5e-08, clump=TRUE)

R2a=2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2b=2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2=R2a/(R2a+R2b)

exposure_dat$F_statistics=R2*(exposure_dat$samplesize.exposure-2)/(1-R2)

exposure_dat=exposure_dat[exposure_dat$F_statistics>10,]

exposure_dat$exposure=stringr::str_sub(exposure_dat$exposure,1,15)

exposure_dat$ENSEMBL =exposure_dat$exposure

library(tidyverse)
exposure_dat= left_join(exposure_dat,df, by = "ENSEMBL")
exposure_dat$SYMBOL
exposure_dat$exposure=NULL
exposure_dat$exposure=exposure_dat$SYMBOL
exposure_dat$SYMBOL=NULL
save(exposure_dat,file ='exposure_dat1.Rdata')

load('exposure_dat.Rdata')

rm(scRNA_EM)
rm(scRNA_compare)
rm(scRNA)
gc()
library(data.table)

#提取结局数据
#outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST90011866")
outcome_dat<-fread("CopyOfGCST90435611.csv")
head(exposure_dat)
setnames(outcome_dat,
         old = c("variant_id","beta", "standard_error", "effect_allele", "other_allele", "effect_allele_frequency", "p_value", "n"),
         new = c("SNP","beta.outcome", "se.outcome", "effect_allele.outcome",
                 "other_allele.outcome", "eaf.outcome", "pval.outcome", "samplesize.outcome"))
# 取交集
common_snp <- merge(exposure_dat,outcome_dat,by="SNP")$SNP
outcome_dat <- outcome_dat[outcome_dat$SNP %in% common_snp, ]
exposure_dat <- exposure_dat[exposure_dat$SNP %in% common_snp,]
outcome_dat$id.outcome <- as.character("gcst90435611")
outcome_dat$outcome <- as.character("ov")

harmonised_dat <- harmonise_data(exposure_dat, outcome_dat,action=2)
## MR分析

mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
{
  mr_res <- mr(dat)
  
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  
  return(mr_res)
}




mr_res <- mr_modified(harmonised_dat, prop_var_explained = T)

save(exposure_dat,outcome_dat,mr_res,harmonised_dat,file ='mr_input_res.Rdata')

load('gcst90435611 mr_input_res.Rdata')
load("mr_input_res.Rdata")
#  异质性 ，用处不大
# 当存在显著的异质性（p<0.05）时，结果不太可靠
heter_dat=mr_heterogeneity(harmonised_dat)
write.csv(heter_dat, file="table.heterogeneity.csv", row.names=F)
heter_dat1<-heter_dat %>% 
  filter(Q_pval < 0.05)
# 水平多效性检验，用处不大
#如果该截距项与0差异很大（pval < 0.05），说明存在水平多效性
pleio_dat=mr_pleiotropy_test(harmonised_dat)
write.csv(pleio_dat, file="table.pleiotropy.csv", row.names=F)
pleio_dat1<-pleio_dat %>% 
  filter(pval < 0.05)


### table
table1 <- mr_res %>% 
  filter(pval < 0.05,
         method %in% c("Wald ratio","Inverse variance weighted")) %>% 
  left_join(exposure_dat, by = "exposure")

table1 <- table1 %>% 
  generate_odds_ratios()%>% 
  mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95),
         `P value` = scales::scientific(pval),
         `PVE` = paste0(sprintf("%.2f",100*pve),"%"),
         `F statistics` = sprintf("%.2f",F_statistics)) %>% 
  dplyr::select(Gene = exposure, `ENSEMBL ID` = ENSEMBL,
                SNP, `Effect allele` = effect_allele.exposure, 
                `OR (95% CI)`, `P value`, 
                PVE, `F statistics`)

## 保存，后面常用
save(table1,file ='table1.Rdata')
load("GCST90435611 table1.Rdata")

library(ggrepel)
### 可视化，如果有:=报错可以跳过这段用下面井号注释掉的

volcano_plot <- function(.data, 
                         number_comparasion = 1,
                         title = "eQTL",
                         col_beta = "b",
                         col_size = "pve",
                         col_label = "exposure",
                         legend.position = "none")
{
  p_thershold <- 0.05/number_comparasion
  
  p <- .data %>% 
    rename(beta := !!col_beta,
           size := !!col_size,
           label := !!col_label) %>% 
    mutate(x = beta,
           y = -log10(pval),
           label = ifelse(pval < p_thershold, label, NA)) %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(size = size), alpha = 0.5, color = "#0072b5") +
    geom_vline(xintercept = 0, linetype = 2)+
    geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6.5),
          legend.position = legend.position)+
    labs(x = "ln(OR)", 
         y = parse(text = "-log[10]*(italic(P)-value)"),
         title = title) +
    scale_size(name = "PVE",
               breaks = c(0.2*1:3)) +
    ggrepel::geom_label_repel(aes(label = label),size = 3)
  plot(p)
}
volcano_plot(mr_res, number_comparasion = 1, title = "eQTL", col_beta = "b", col_size = "pve", col_label = "exposure")
colnames(mr_res)

mr_res %>% 
  filter(method %in% c("Wald ratio","Inverse variance weighted")) %>% 
  volcano_plot(number_comparasion = 1)

#### 如果有:=报错可以用下面 ，去掉井号
#  library(dplyr)
#  library(ggplot2)
#  library(ggrepel)

#  volcano_plot <- function(.data, 
#                           number_comparasion = 1,
#                           title = "eQTL",
#                          legend.position = "none") {
#   
#   p_thershold <- 0.05/number_comparasion
#   
#  p <- .data %>% 
#    mutate(y = -log10(pval),  
#           label = ifelse(pval < p_thershold, exposure, NA)) %>% 
#    ggplot(aes(x = b, y = y)) +
#    geom_point(aes(size = pve), alpha = 0.5, color = "#0072b5") +
#    geom_vline(xintercept = 0, linetype = 2) +
#    geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
#    theme_classic() +
#    theme(panel.grid = element_blank(),
#          legend.title = element_text(size = 6.5),
#         legend.text = element_text(size = 6.5),
#         legend.position = legend.position) +
#    labs(x = "ln(OR)", 
#         y = parse(text = "-log[10]*(italic(P)-value)"),
#        title = title) +
#    scale_size(name = "PVE",
#                breaks = c(0.2*1:3)) +
#    ggrepel::geom_label_repel(aes(label = label), size = 3)
#  plot(p)
# }

# 画图
# mr_res %>% 
#   filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>% 
#  volcano_plot(number_comparasion = 1)

## 其他可视化（一般不太好用，因为一个基因只有很有限的eqtl！！
# 取一个有多个工具变量的基因
table(mr_res$exposure)
mr_res1=mr_res[mr_res$exposure=='CLSTN3',]
harmonised_dat1=harmonised_dat[harmonised_dat$exposure=='CLSTN3',]

#绘制散点图
#每一个点其实代表的就是一个IV，每个点上的线实际反映的是95%置信区间，横坐标是SNP对暴露的效应，纵坐标是SNP对结局的效应
mr_scatter_plot(mr_res1, harmonised_dat1)

#森林图 每一条水平实线反映的是单个SNP利用Wald ratio方法估计出来的结果
res_single=mr_singlesnp(harmonised_dat1)
mr_forest_plot(res_single)

#漏斗图  这和我们之前发现这些SNP间的异质性很大的结果相符。
mr_funnel_plot(singlesnp_results = res_single)

#留一法敏感性分析
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))


mr_res2=mr_res[mr_res$exposure %in% table1$Gene,]
result_or=generate_odds_ratios(mr_res2)
write.table(result_or[,4:ncol(result_or)],"gcst90435611OR.txt",row.names = F,sep = "\t",quote = F)

#将统计结果绘制森林图


library(grid)
library(forestploter)

mydata
mydata=read.table("gcst90435611OR.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
str(mydata)
mydata <- mydata %>%
  filter(method == "Inverse variance weighted" & pval < 0.05)

mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                                                                                    mydata$or_uci95))
forest(mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)

common_snp <- merge(mydata,mr_res,by="exposure")$exposure
mr_res <- mr_res[mr_res$exposure %in% common_snp, ]
exposure_dat <- exposure_dat[exposure_dat$SNP %in% common_snp,]


library(ggplot2)
library(dplyr)
library(ggrepel)

volcano_plot <- function(.data, 
                         number_comparasion = 1,
                         title = "eQTL",
                         col_beta = "b",
                         col_size = "pve",
                         col_label = "exposure",
                         legend.position = "none",
                         max_labeled_genes = 10) {  
  p_threshold <- 0.05 / number_comparasion
  
  # 重命名列
  data_renamed <- .data %>% 
    rename(beta = !!col_beta,
           size = !!col_size,
           label = !!col_label)
  
  # 确保 CLSTN3 总是被标记
  clstn3_data <- data_renamed %>% filter(label == "CLSTN3")
  if (nrow(clstn3_data) > 1) {
    clstn3_data <- clstn3_data %>% 
      slice_max(n = 1, order_by = -log10(pval))  # 选择 p 值最小的条目
  }
  
  # 选择其他显著的基因
  significant_genes <- data_renamed %>% 
    filter(pval < p_threshold & label != "CLSTN3") %>% 
    slice_max(n = max_labeled_genes - 1, order_by = -log10(pval))  # 选择其他显著的基因
  
  highlight_genes <- unique(c(clstn3_data$label, significant_genes$label))
  
  p <- data_renamed %>% 
    mutate(x = beta,
           y = -log10(pval),
           significant = label %in% highlight_genes,
           below_threshold = y < -log10(p_threshold)) %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(size = size, color = significant), alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "#0072b5", "FALSE" = "#CCCCCC")) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = -log10(p_threshold), linetype = 2) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6.5),
          legend.position = legend.position) +
    labs(x = "ln(OR)", 
         y = parse(text = "-log[10]*(italic(P)-value)"),
         title = title,
         color = "Significance") +
    scale_size(name = "PVE",
               breaks = c(0.2 * 1:3)) +
    ggrepel::geom_label_repel(data = . %>% filter(significant & !below_threshold), 
                              aes(label = label), 
                              size = 3, 
                              max.overlaps = 100,  # 增加 max.overlaps 参数
                              point.padding = 0.5,  # 增加点的填充
                              box.padding = 0.35)  # 增加标签的填充
  
  print(p)
}

# 调用函数
volcano_plot(mr_res, number_comparasion = 1, title = "eQTL", col_beta = "b", col_size = "pve", col_label = "exposure")





