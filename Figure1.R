


#查看工作空间
getwd()
setwd(""path/to/your/ukb/data"")
#R包下载、检验、安装
if (!require("haven")) install.packages("haven")
if (!require("tableone")) install.packages("tableone");library("tableone")
if (!require("survival")) install.packages("survival");library("survival")
if (!require("dplyr")) install.packages("dplyr");library("dplyr")
if (!require("VIM")) install.packages("VIM");library("VIM")
if (!require("mice")) install.packages("mice");library("mice")
if (!require("mediation")) install.packages("mediation");library("mediation")
if (!require("ggplot2")) install.packages("ggplot2");library("ggplot2")
if (!require("rms")) install.packages("rms");library("rms")
if (!require("jstable")) install.packages("jstable");library("jstable")
if (!require("forestploter")) install.packages("forestploter");library("forestploter")



## —————————————————— 5.缺失数据填补与回归分析 ————————————————————####

#读取需多重插补的数据
tt3<-fread("occontrol.csv")       ## 
names(tt3)
str(datana)
datana<-tt3
#检查缺失值情况
var <- names(datana)
tab <- list()
total <- nrow(datana)
for (i in 1:length(var)){
  tab[i]=data.frame(Variables=var[i],
                    Total=total,
                    Freq=total-sum(is.na(datana[,var[i]])),
                    Missing=sum(is.na(datana[,var[i]])),
                    miss_p=sprintf("%0.2f",sum(is.na(datana[,var[i]]))/total*100)) %>% list()}
na<-do.call(rbind, tab)

#利用VIM包检查缺失值情况可视化
library(VIM)
mice_plot <- aggr(datana, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(datana), cex.axis=.7,
                  gap=5, ylab=c("Missing data","Pattern"),)


#对数据中的分类变量进行因子化
datana<-datana[,c(-15,-16)]
colnames(datana)
v<- c("group","Age","BMI","Smoking","Alcohol","Oral_Contraceptive",
      "Infertility_History","PCOS","PID","Hysterectomy","BRCA_Mutation","Anxiety"
      ,"Hormone_Use","Family_History","Me0pausal_Status")

v<- c("group","Age","Smoking","Alcohol","Oral_Contraceptive",
      "Infertility_History","PCOS","PID","Hysterectomy","BRCA_Mutation","Anxiety"
      ,"Hormone_Use")
library(data.table)

setDT(datana)
datana[, (v) := lapply(.SD, factor), .SDcols = v]
str(datana)

#根据不同类型数据选择不同插补方法
meth <- make.method(datana)
meth <- c("age" = "default", "dp" = "default", "chronic_num" = "default", "other" = "default")
# 使用索引向量赋值
meth[c("BMI", "Depression_Score_PHQ9", "CRP_mg","CA125_U","HE4_pmol")] <- "pmm"                                   #定量数据
meth[c("group","Age","Oral_Contraceptive",
       "Infertility_History","PCOS","PID","Hysterectomy","BRCA_Mutation","Anxiety"
       ,"Hormone_Use")] <- "logreg"   #二分
meth[c("Smoking","Alcohol")] <- "rf"    #多分类


# m为模型数量，通常为5，maxit为模型迭代次数，通常为50，此处为节省时间进行缩减
library(mice)
imputed<-mice(datana, m=3, maxit = 20, seed = 314)
imputed$method
complete_data <- complete(imputed, 1)

fit<-with(imputed,coxph(Surv(fultime,stat)~dp+
                          age+education+
                          drink+smoke+hypertension+diabetes+chronic_num))
pooled<-pool.syn(fit)
summary(pooled)

## —————————————————— 5.箱线图 ————————————————————####
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggpubr)) install.packages("ggpubr")
library(ggplot2)
library(ggpubr)

t_test_result <- t.test(complete_data$dp ~ complete_data$ec, data = complete_data)
print(t_test_result)
complete_data$ec <- as.factor(complete_data$ec)
complete_data$dp <- as.numeric(complete_data$dp)
# 绘制箱线图并添加p值
ggboxplot(complete_data, x = "ec", y = "dp",
          color = "ec", palette = "jco",
          add = "jitter") +
  stat_compare_means(method = "t.test", label = "p.signif")


### ——— RCS ——— ####

## 利用ggplot2包和rms包，绘制RCS图
library(ggplot2)
library(rms)

## 为后续程序设定数据环境
dd <- datadist(complete_data)
options(datadist='dd')

#拟合回归方程，选择4个节点
# 检查 drink 变量的水平
# 重新设置 datadist
dd <- datadist(complete_data)
options(datadist = 'dd')
table(complete_data$drink)
rcs1<-cph(Surv(fultime,stat) ~ rcs(dp,4)+
            dp+education+drink+smoke+hypertension+hyperlipid+diabetes+chronic_num,
          data = complete_data)
anova(rcs1)               # P Nonlinear<0.05为非线性趋势

#计算HR值
HR<-Predict(rcs1,dp,fun=exp,ref.zero = TRUE)      

#绘制RCS曲线
p1<-ggplot()+
  geom_line(data=HR, aes(dp,yhat),      # 添加线条
            linetype=1,size=1,            # 设置线的类型为实线、线的粗细为1
            alpha = 0.9,colour="red")+    # 线条颜色为红色，透明度为0.9
  geom_ribbon(data=HR, aes(dp,ymin = lower, ymax = upper),  # 添加置信区间
              alpha = 0.3,fill="red")+    # 线条颜色为红色，透明度为0.3
  geom_hline(yintercept=1,                # 添加一条水平线y=1
             linetype=2,size=1)+          # 设置线的类型为虚线、线的粗细为1
  geom_vline(xintercept=54.124 ,          # 添加一条垂直线 X=23.611 (查表HR=1对应的VD值)
             linetype=2,size=1,color="black")+   
  theme_classic()+                        # 设置图形为经典主题
  labs(title = "CVD", x="PM2.5", y="HR (95%CI)")+   # 设置图形标题、X轴和Y轴标签
  geom_text(aes(x = 70, y = 2.5),                   # 添加文本
            label = expression(paste(italic("P"), " for nonlinear < 0.001")),
            color = "black")

p1

### —————————————————— 亚组分析、交互作用（森林图） ————————————————————####

### ———— 亚组分析、交互作用 ————####

## 确定需要进行的亚组分析变量
# select表示从complete_data2数据集中选择需要的变量
# mutate表示对变量进行转换和重新编码
library(dplyr)
df <- complete_data %>%
  dplyr::select(fultime, stat, dp, age, gender, residence, smoke, drink, hypertension, diabetes) %>%
  mutate(
    stat = as.integer(stat == 1),
    age = factor(ifelse(age < 60, "<60", ">=60"), levels = c("<60", ">=60")),
    gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male")),
    residence = factor(residence, levels = c(0, 1), labels = c("Rural", "Urban")),
    smoke = factor(smoke, levels = c(0, 1, 2), labels = c("Never", "Former","Current")),
    drink = factor(drink, levels = c(0, 1, 2), labels = c("Never", "Former","Current")),
    hypertension = factor(hypertension, levels = c(0, 1), labels = c("No", "Yes")),
    diabetes = factor(diabetes, levels = c(0, 1), labels = c("No", "Yes"))
  )


## 使用jstable包进行亚组分析
# 利用TableSubgroupMultiCox函数进行Cox生存分析
# formula = Surv(fultime,stat)~VD指定生存分析的模型
# var_subgroups指定将哪些变量设置为亚组
library(jstable) 
res<-TableSubgroupMultiCox(formula=Surv(fultime,stat)~dp,
                           var_subgroups=c("age","gender","residence","smoke","drink","hypertension","diabetes"),
                           data=df)
res
write.csv(res,"subgroup_result.csv")     ## 输出亚组分析和交互作用结果


### ———— 森林图 ————####

## 整理res数据来绘制森林图 
plot_df <- res          ## 对res数据集重新命名为plot_df，变量名存在空格，因此对变量名重新命名
names(plot_df)<-c("Variable","Count" ,"Percent", "PointEstimate","Lower" ,"Upper" ,"Pvalue" ,"Pforinteraction" )  
plot_df[, c(2, 3, 7, 8)][is.na(plot_df[, c(2, 3, 7, 8)])] <- ""  ## 选取第2、3、7、8列的数据在森林图中显示，并替换缺失值NA为一个空格字符
plot_df$""<- paste(rep("", nrow(plot_df)), collapse = "")        ## 添加空白列，用于存放森林图的图形部分
plot_df[,4:6]<-apply(plot_df[,4:6],2,as.numeric)                 ## 将第4到第6列的数据转换为数值型
plot_df$"HR(95%CI)"<-ifelse(is.na(plot_df$"PointEstimate"),"",
                            sprintf("%.2f(%.2fto%.2f)",
                                    plot_df$"PointEstimate",plot_df$Lower,plot_df$Upper))#计算HR(95%CI)，以便显示在图形中
plot_df                ## 查看数据


## 利用forestploter包绘制森林图
library(forestploter)
p<-forest(
  data=plot_df[,c(1,2,3,9,10,7,8)],  ## 选择需要用于绘图的列
  lower=plot_df$Lower,               ## 置信区间下限
  upper=plot_df$Upper,               ## 置信区间上限
  est=plot_df$`PointEstimate`,       ## 点估计值
  ci_column=4,                       ## 点估计对应的列
  ref_line=1,                        ## 设置参考线位置
  xlim=c(0,5)                        ## x轴的范围
)
plot(p)                              ## 输出森林图