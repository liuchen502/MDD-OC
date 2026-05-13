## 孟德尔后：验证集、双向孟德尔、共定位分析
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyMjM2NjI5OTQyQHFxLmNvbSIsImlhdCI6MTc2MTM3Mjc1OCwiZXhwIjoxNzYyNTgyMzU4fQ.FoVW8xinXScCOIh0Aji60Yn9Z5BFj946TI2JEucWmPf9h0S8d0hIkpOwTTJu-yx-n3nNblsdIBGb7deihCZxVBQKj44WtXbcdVtwYwXV8dDb7cCBLU6Ij0vuflOV3H9ZXCihQNrStyfhMUp7dlJ8sMRZG5oIlv99VyVaj3E6P9xGSOFxDe354bDugt6AkJCF3uDtuDBu0bVoxHtofe6ttllkRw35vAfk1LKDgNXJzf0b84RhxMvZG_-6y5wqYjrTExzPY9XnPGzsnQ9AV6bPiSxcZMfb-AOHp3oXaVcHfnNiPLTOtZX8giTgeETyPiNXoFWTXqjHtOUpbaOkRckE9g") #这个token就是保存下来一长串的值，复制进去就可以了
OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyMjM2NjI5OTQyQHFxLmNvbSIsImlhdCI6MTc2MTM3Mjc1OCwiZXhwIjoxNzYyNTgyMzU4fQ.FoVW8xinXScCOIh0Aji60Yn9Z5BFj946TI2JEucWmPf9h0S8d0hIkpOwTTJu-yx-n3nNblsdIBGb7deihCZxVBQKj44WtXbcdVtwYwXV8dDb7cCBLU6Ij0vuflOV3H9ZXCihQNrStyfhMUp7dlJ8sMRZG5oIlv99VyVaj3E6P9xGSOFxDe354bDugt6AkJCF3uDtuDBu0bVoxHtofe6ttllkRw35vAfk1LKDgNXJzf0b84RhxMvZG_-6y5wqYjrTExzPY9XnPGzsnQ9AV6bPiSxcZMfb-AOHp3oXaVcHfnNiPLTOtZX8giTgeETyPiNXoFWTXqjHtOUpbaOkRckE9g"
# 记得先改镜像
library(tidyverse)

load('exposure_dat.Rdata')

gene=read.table('./OR.txt',sep = '\t',header = T)
gene<-mydata
gene=gene$exposure

exposure_dat=exposure_dat[exposure_dat$exposure %in% gene,]

### 验证集----

#提取结局数据,换一个队列
library(TwoSampleMR)
outcome_dat<-fread("GCST90011821_buildGRCh37.csv")
head(outcome_dat)
setnames(outcome_dat,
         old = c("snp1","BETA", "SE", "effect_allele", "other_allele", "EAF", "P", "n1"),
         new = c("SNP","beta.outcome", "se.outcome", "effect_allele.outcome",
                 "other_allele.outcome", "eaf.outcome", "pval.outcome", "samplesize.outcome"))
# 取交集
common_snp <- merge(exposure_dat,outcome_dat,by="SNP")$SNP
outcome_dat <- outcome_dat[outcome_dat$SNP %in% common_snp, ]
exposure_dat <- exposure_dat[exposure_dat$SNP %in% common_snp,]
outcome_dat$id.outcome <- as.character("gcst90011821")
outcome_dat$outcome <- as.character("ov")
# 取交集
harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)



## MR

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




mr_res_vali <- mr_modified(harmonised_dat, prop_var_explained = T)

save(mr_res_vali,harmonised_dat,file ='gcst90011821 mr_input_res_vali.Rdata')
load("gcst90011821 mr_input_res_vali.Rdata")

result_or=generate_odds_ratios(mr_res_vali)
write.table(result_or[,4:ncol(result_or)],"gcst90011821 OR_vali.txt",row.names = F,sep = "\t",quote = F)

#将统计结果绘制森林图

library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
library(TwoSampleMR)
library(forestplot)
library(MVMR)
library(grid)
library(forestploter)


mydata=read.table("gcst90011821 OR_vali.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
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


## 其他可视化（一般不太好用，因为一个基因只有很有限的eqtl
mr_res1=mr_res_vali[mr_res_vali$exposure=='CLSTN3',]
harmonised_dat1=harmonised_dat[harmonised_dat$exposure=='CLSTN3',]
# 遍历每个基因


### 反向孟德尔
# 暴露(ieu-a-1120)
getwd()
setwd("/home/xing54321/MDD")
library(data.table)
expo <- fread("CopyOfGCST90435611.csv")
colnames(expo)
## 4. 选工具变量 --------------------------------------------------------
# 4.1 全基因组显著
sig <- expo %>% filter(p_value < 5e-4)
head(sig)
setnames(sig,
         old = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele","p_value","effect_allele_frequency","beta","standard_error","n"),
         new = c("SNP","CHR.exposure","base_pair_location.exposure","effect_allele.exposure", "other_allele.exposure",
                 "pval.exposure","beta.exposure","eaf.exposure","se.exposure","samplesize.exposure"))

colnames(sig)
# 4.2 clumping（用 1000G EUR 参考，默认 r2=0.001, kb=10000）
instruments <- clump_data(sig,
                          clump_r2 = 0.001,
                          clump_kb = 10000,
                          pop = "EUR")

## 6. 保存备用 -----------------------------------------------------------
saveRDS(instruments, "GCST90435611_EUR_exposure.rds")
write_tsv(instruments, "GCST90435611_EUR_exposure.tsv")
bimr_oc<-sig

#提取结局数据
outcome_gene<- extract_outcome_data(snps=bimr_oc$SNP, outcomes="eqtl-a-ENSG00000139182")

# 取交集
bimr_oc =bimr_oc [bimr_oc $SNP %in% outcome_gene$SNP,]

bimr_oc$id.exposure<-"gcst90435611"
bimr_oc$exposure<-"exposure"
harmonised_oc_gene <- harmonise_data(bimr_oc, outcome_gene)

bimr_mr_oc_gene <- mr(harmonised_oc_gene)

result_oc=generate_odds_ratios(bimr_mr_oc_gene)
write.table(result_oc[,4:ncol(result_oc)],"clstn3 bi_OR-反向孟德尔.txt",row.names = F,sep = "\t",quote = F)

#将统计结果绘制森林图


library(grid)
library(forestploter)


mydata=read.table("clstn3 bi_OR-反向孟德尔.txt",header = T,sep = "\t")
## !!
mydata$outcome='clstn3'
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,12,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =6 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)



### 共定位分析----

load('mr_input_res.Rdata')

#如果表型是二分类变量（case和control），输入文件二选一：
#1）rs编号`rs_id`、P值`pval_nominal`、SNP的效应值`beta`、效应值方差`varbeta`；（推荐）
#2）rs编号`rs_id`、P值`pval_nominal`、case在所有样本中的比例`s`，MAF也要，写在list最后

#如果表型是连续型变量，输入文件三选一：
#1）rs编号`rs_id`、P值`pval_nominal`、表型的标准差`sdY`；
#2）rs编号`rs_id`、P值`pval_nominal`、效应值`beta`,效应值方差 `varbeta`, 样本量`N`,次等位基因频率 `MAF`；
#3）rs编号`rs_id`、P值`pval_nominal`、次等位基因频率 `MAF`；(推荐)


## 下载
getwd()
data <- vcfR::read.vcfR("CLSTN3 eqtl-a-ENSG00000139182.vcf.gz")

# 1.SNP ID rs开头
# 2.effect allele，此处相当于ALT
# 3.other allele，此处相当于REF
# 4.beta
# 5.se
# 6.pval


#整理数据
## 处理gt
gt=data@gt
gt=as.data.frame(gt)

colnames(gt)
gt$FORMAT[1]

library(tidyverse)

gt$`eqtl-a-ENSG00000139182`[1]

##!!!
gt=separate(gt,col='eqtl-a-ENSG00000139182',into = c('ES', 'SE',
                                                     'LP','AF','SS',
                                                     'ID'),sep = '\\:')

gc()



gt=na.omit(gt)
colnames(gt)=c('format','beta','se','logpvalue','eaf','samplesize','snp')
gt$beta=as.numeric(gt$beta)
gt$se=as.numeric(gt$se)
gt$logpvalue=as.numeric(gt$logpvalue)
gt$eaf=as.numeric(gt$eaf)
gt$samplesize=as.numeric(gt$samplesize)

gc()
gt$format=NULL


fix=data@fix
fix=as.data.frame(fix)
colnames(fix)
colnames(fix)=c('chr','pos','snp','ref','alt')
fix=fix[,1:5]


## 合并gt fix
eqtl=left_join(fix,gt,by='snp')
eqtl=na.omit(eqtl)
## 查找染色体和位置
#22
#39605007

#MAFs代表次要等位基因的频率，因此其范围在0到0.5之间。值大于0.5的所有EAF通过从1中减去它们而转化为MAF。
eqtl$maf = ifelse(eqtl$eaf < 0.5, 
                  eqtl$eaf,
                  1 - eqtl$eaf)
eqtl$eaf=NULL


#	
#22
#39605007
#CLSTN3
eqtl=eqtl[eqtl$chr==12,]
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)

eqtl$pos=as.numeric(eqtl$pos)

## 上下1mkb
eqtl=eqtl[eqtl$pos > 7129093-1000000 ,]
eqtl=eqtl[eqtl$pos < 7158945+1000000 ,]







#MIAT
eqtl=eqtl[eqtl$chr==22,]
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)

eqtl$pos=as.numeric(eqtl$pos)

## 上下1mkb
eqtl=eqtl[eqtl$pos > 27053446-1000000 ,]
eqtl=eqtl[eqtl$pos < 27072441+1000000 ,]

#PRSS23
eqtl=eqtl[eqtl$chr==11,]
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)

eqtl$pos=as.numeric(eqtl$pos)

## 上下1mkb
eqtl=eqtl[eqtl$pos > 86791059-1000000 ,]
eqtl=eqtl[eqtl$pos < 86952910+1000000 ,]

#PIM1
eqtl=eqtl[eqtl$chr==6,]
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)

eqtl$pos=as.numeric(eqtl$pos)

## 上下1mkb
eqtl=eqtl[eqtl$pos > 37170152-1000000 ,]
eqtl=eqtl[eqtl$pos < 37175428+1000000 ,]

my_eqtl=eqtl[,c('snp','p_value','maf')]

colnames(my_eqtl)=c('snp','pvalues','MAF')
my_eqtl=na.omit(my_eqtl)

my_eqtl=my_eqtl[my_eqtl$MAF>0 ,]

##-----接下来gwas疾病


library(TwoSampleMR)
## ！！！！！
outcome <- fread("CopyOfGCST90435611.csv")
colnames(outcome)
## 4. 选工具变量 --------------------------------------------------------
# 4.1 全基因组显著
setnames(outcome,
         old = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele","p_value","effect_allele_frequency","beta","standard_error","n"),
         new = c("SNP","CHR.outcome","base_pair_location.outcome","effect_allele.outcome", "other_allele.outcome",
                 "pval.outcome","beta.outcome","eaf.outcome","se.outcome","samplesize.outcome"))
colnames(outcome)[1]<-"snp"
common_snp <- merge(my_eqtl,outcome,by="snp")$snp

coloc_or_dat<-outcome[outcome$snp %in% common_snp, ]
colnames(coloc_or_dat)
coloc_or_dat %>%  mutate(chr.outcome = as.numeric(CHR.outcome),
                         pos.outcome = as.numeric(base_pair_location.outcome),
                         outcome = "ov", 
                         id.outcome = "gcst90435611")

##gwas:疾病
gwas=coloc_or_dat
gwas$beta=as.numeric(gwas$beta.outcome)
gwas$se=as.numeric(gwas$se.outcome)
#X=gwas$beta/gwas$se
#P=2*pnorm(q=abs(X), lower.tail=FALSE) #Z=1.96  P=0.05
#gwas$pvalue=P
## check Lp=log10(P)
gwas$varbeta=(gwas$se)^2

gwas=gwas[,c('snp','pval.outcome',"beta",'varbeta')]

colnames(gwas)=c('snp','pvalues','beta','varbeta')

gwas=na.omit(gwas)

##---开始coloc
library(coloc)

input <- merge(my_eqtl, gwas, by="snp", all=FALSE, suffixes=c("_eqtl","_gwas"))


head(input)
library(coloc)

# 如报错，去除重复值
input=input[!duplicated(input$snp),]
library(coloc)

dataset1 <- list(
  snp   = input$snp,
  pvalues = input$pvalues_gwas,
  type    = "cc",
  s       = 0.05,
  N       = 391798
)

dataset2 <- list(
  snp   = input$snp,
  pvalues = input$pvalues_eqtl,
  type    = "quant",
  N       = 31684
)

res <- coloc.abf(dataset1, dataset2, MAF = input$MAF)
need_result=res$results %>% filter(SNP.PP.H4 > 0.95)
need_result
print(res)
saveRDS(res,"CLSTN3 共定位.rds")
a<-res$results
write.table(a, file = "CLSTN3 共定位.tsv", sep = "\t", row.names = TRUE,col.names = TRUE, quote = FALSE)
## 或者用下面的方法

#my_eqtl=as.list(my_eqtl)
#my_eqtl[['type']]='quant'
#my_eqtl[['N']]=31430

#gwas=as.list(gwas)
#gwas[['type']]='cc'
#gwas[['N']]=12653
#gwas$pvalues=NULL
#coloc.abf(dataset1 = gwas,dataset2 = my_eqtl)


## 解读------
#H0：该区域的两个性状都没有遗传关联
#H1 / H2：只有表型1或表型2在该区域具有遗传关联
#H3：两个特征都相关，但因果变量不同
#H4：两个特征都相关并且共享一个因果变量
# 安装并加载 ggplot2 包
if (!require(ggplot2)) {
  install.packages("ggplot2")
}
library(ggplot2)
a$snp <- as.factor(a$snp)
# 绘制两个性状的 p 值的曼哈顿图
# 绘制曼哈顿图，突出显示特定的 SNP
p1_plot <- ggplot(a, aes(x = as.numeric(snp), y = -log10(pvalues.df1))) +
  geom_point(aes(color = as.factor(snp)), alpha = 0.6) +  # 绘制所有点，颜色按 SNP 分类
  geom_point(data = subset(a, snp == "rs3759416"), aes(x = as.numeric(snp), y = -log10(pvalues.df1)), 
             color = "red", size = 4, shape = 18) +  # 突出显示特定的 SNP
  theme_minimal() +
  labs(title = "Manhattan Plot for Trait 1",
       x = "SNP Index",
       y = "-log10(p-value)",
       color = "SNP") +
  theme(legend.position = "none")

# 显示图形
print(p1_plot)


p2_plot <- ggplot(a, aes(x = as.numeric(snp), y = -log10(pvalues.df2))) +
  geom_point(aes(color = as.factor(snp)), alpha = 0.6) +  # 绘制所有点，颜色按 SNP 分类
  geom_point(data = subset(a, snp == "rs3759416"), aes(x = as.numeric(snp), y = -log10(pvalues.df2)), 
             color = "red", size = 4, shape = 18) +  # 突出显示特定的 SNP
  theme_minimal() +
  labs(title = "Manhattan Plot for Trait 2",
       x = "SNP Index",
       y = "-log10(p-value)",
       color = "SNP") +
  theme(legend.position = "none")

# 显示图形
print(p2_plot)


# 绘制两个性状的贝叶斯因子的比较图
bf_plot <- ggplot(a, aes(x = as.numeric(snp))) +
  geom_point(aes(y = lABF.df1, color = "oc lABF"), alpha = 0.6) +
  geom_point(aes(y = lABF.df2 + 5, color = "eqtl clstn3 lABF"), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Comparison of oc and eqtl clstn3 lABF",
       x = "SNP Index",
       y = "lABF",
       color = "Legend") +
  theme(legend.position = "right")

# 显示图形
print(bf_plot)


# 绘制后验概率的图
pp_plot <- ggplot(a, aes(x = as.numeric(snp), y = SNP.PP.H4, color = as.factor(snp))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Posterior Probability Plot",
       x = "SNP Index",
       y = "Posterior Probability",
       color = "SNP") +
  theme(legend.position = "none")

# 显示图形
print(pp_plot)


A<-`CLSTN3 共定位`
b<-A$results
str(b)
library(dplyr)

# 筛选SNP.PP.H4较高的SNP（如前100名，避免点过于拥挤）
top_snps <- b %>%
  arrange(desc(SNP.PP.H4)) %>%
  slice_head(n = 100)

# 绘图
ggplot(top_snps, aes(x = -log10(pvalues.df1), y = -log10(pvalues.df2))) +
  geom_point(aes(size = SNP.PP.H4, color = SNP.PP.H4)) +  # 大小和颜色均表示共定位概率
  scale_color_viridis_c(option = "plasma", name = "SNP.PP.H4") +  # 渐变色
  geom_text(aes(label = ifelse(SNP.PP.H4 > 0.8, snp, "")),  # 标记高概率SNP
            hjust = 1.1, vjust = 0, size = 3) +
  labs(
    x = expression(-log[10]("p-value (df1)")),
    y = expression(-log[10]("p-value (df2)")),
    title = "Colocalization Probability vs. Association Signals",
    size = "SNP.PP.H4"
  ) +
  theme_minimal() +
  theme(legend.position = "right")



