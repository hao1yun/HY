#usethis::edit_r_environ()获取token
#ieugwasr::get_opengwas_jwt()返回token

#hypothyroidism--ukb-b-19732
#polycystic ovarian syndrome--finn-b-E4_POCS

setwd("E:\\杨旭-文章\\论文\\第2篇\\IBS-N")#run/ctrl+enter
getwd()
.libPaths()#R包安装路径

##安装TwoSampleMR包----
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("devtools")
library(devtools)
devtools::install_github("MRCIEU/TwoSampleMR")
#本地安装
##安装MRPRESSO包----
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
#devtools::package_info("ieugwasr",dependencies = F)查看包版本
#devtools::package_info("TwoSampleMR",dependencies = F)
library(TwoSampleMR)
library(MRPRESSO)
#usethis::edit_r_environ()该句跳转输入token(JWT)
# 加载ieugwasr包
library(ieugwasr)
# 获取OpenGWAS的JWT（验证身份）
#token <- ieugwasr::get_opengwas_jwt()
##1.提取暴露[exp]数据------
exp <-extract_instruments(outcomes = 'BBJ-a-119',p1 = 5e-08,
                           clump = TRUE,r2 = 0.001,kb = 10000,
                          access_token = NULL )#提取暴露数据
exp<- extract_instruments(outcomes='ebi-a-GCST90016564')
exp<-extract_instruments(outcomes= 'ebi-a-GCST90016564',
  p1 = 5e-06,
  clump = TRUE,
  r2 = 0.001,
  kb = 5000,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  force_server = FALSE
)
dim(exp)#【行, 列】
colnames(exp)
write.csv(exp,'exp.csv')
save(exp,file="exp.rda")
##2.提取结局[out]数据------
out<-extract_outcome_data(snps = exp$SNP,
                          outcomes = 'BBJ-a-119',
                          proxies=TRUE)#提取结局数据
write.csv(out,'out.csv')
save(out,file="out.rda")
##3.等位基因对齐 + 去除回文序列------
####action = 3，意味着剔除回文序列
dat<-harmonise_data(exp,out,action = 3) 
head(dat)
write.csv(dat,'dat.csv')
#save(dat,file="dat.rda")
##excel里操作
#①phenoscanner【二三假设】
#②如缺失暴露的样本量，需导出补齐samplesize.exposure
dat <- read.csv('dat.csv')
##4.计算F-statistics------
###如缺少samplesize.exposure，需自动导出补齐[上述步骤]
PVEfx <- function(BETA, SE, N){
  pve <- (BETA^2)/((BETA^2) + ((SE^2)*N))
  return(pve) 
}
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$se.exposure, N = dat$samplesize.exposure)
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE))  #F值

##如上数据缺失算不出F值后，使用b^2/se^2
dat$`F`<-((dat$beta.exposure)^2)/((dat$se.exposure)^2)

write.csv(dat,'datF.csv')
save(dat,file="datF.rda")
##5.MR分析------
#查看TwoSample内置的MR方法
mr_method_list()
results<-mr(dat,method_list = c("mr_ivw","mr_egger_regression",
                                "mr_weighted_median",
                                "mr_simple_mode","mr_weighted_mode"))
results
write.csv(results,'results.csv')
save(results,file="results.rda")
#(OR,95%CI_or,p_value)--二分类变量【1】
#(beta,95%CI_beta,p_value)--连续性变量【0】
##6.二分类变量--OR值------
OR <-generate_odds_ratios(results) 
OR
write.csv(OR,'OR.csv')
save(OR,file="OR.rda")

##7.异质性检验------
heterogeneity <- mr_heterogeneity(dat)  #异质性分析啦
heterogeneity

write.csv(heterogeneity,'heterogeneity.csv')
save(heterogeneity,file="heterogeneity.rda")
##8.水平多效性【MR-Egger截距】------
pleio <- mr_pleiotropy_test(dat)
pleio
write.csv(pleio,'pleio.csv')
save(pleio,file="pleio.rda")
##9.水平多效性【MR-PRESSO】------

#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)   #1水平多效性  2.挑出离群snp   3.计算跳出  /不跳出 离群snp 结果是否差异
mr_presso(BetaOutcome = 'beta.outcome',
          BetaExposure = 'beta.exposure', 
          SdOutcome = 'se.outcome', 
          SdExposure = 'se.exposure', 
          data = dat, OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 1000, seed = NULL)
mr_presso##word/txt导出

library(ggsci)
library(ggplot2)
##10.留一法+作图------
###逐个剔除每个SNP，计算剩余SNP的效应
single <- mr_leaveoneout(dat)
mr_leaveoneout_plot(single)
dev.off()#清空画板
##导出图
##11.散点图------
# 散点图
scatter_plot<-mr_scatter_plot(results,dat)[[1]]+
  scale_color_lancet()+
  scale_fill_lancet()+
  theme_bw()
# 更改散点图的点的颜色
scatter_plot[["layers"]][[3]][["aes_params"]]$colour<-"black"
# 更改散点图的点的透明度
scatter_plot[["layers"]][[3]][["aes_params"]]$alpha<-0.5
# 更改散点图的横坐标名字
scatter_plot[["labels"]][["x"]]<-"SNP effect on irritable bowel syndrome (IBS)"
# 更改散点图的纵坐标名字
scatter_plot[["labels"]][["y"]]<-"SNP effect on Neuroticism"
scatter_plot
dev.off()#清空画板
##12.森林图---------
forest_plot<-mr_forest_plot(mr_singlesnp(dat))[[1]]+
  scale_color_lancet()+
  scale_fill_lancet()+
  theme_bw()+
  theme(legend.position = 'none')
forest_plot
dev.off()#清空画板
##13.漏斗图------
funnel_plot<-mr_funnel_plot(mr_singlesnp(dat,all_method=c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_simple_mode","mr_weighted_mode")))[[1]]+
  theme_bw()+
  scale_color_lancet()+
  scale_fill_lancet()
# 更改漏斗图点的颜色
funnel_plot[["layers"]][[1]][["aes_params"]]$colour <- "black"
# 更改点的透明度
funnel_plot[["layers"]][[1]][["aes_params"]]$alpha<-0.5
funnel_plot
dev.off()#清空画板
##14.密度图----
density_plot<-mr_density_plot(mr_singlesnp(dat),results)[[1]]+
  theme_bw()+
  scale_color_lancet()+
  scale_fill_lancet()
# 更改密度图点的颜色
density_plot[["layers"]][[3]][["aes_params"]][["colour"]]<-"black"
# 更改密度图点的透明度
density_plot[["layers"]][[3]][["aes_params"]]$alpha<-0.5
# 更改密度图点的图例的名字
density_plot[["labels"]][["colour"]]<-"Method"
# 更改密度图点的横轴标的名字
density_plot[["labels"]][["y"]][[1]]<-"Density"
density_plot
dev.off()#清空画板
##15. Contamination mixture法----------
if (!requireNamespace("MendelianRandomization"))
  install.packages("MendelianRandomization")
library(MendelianRandomization)
library(TwoSampleMR)

# 正常两样本流程
# 将mr_keep为TRUE的snp行留下来
dat<-dat[dat$mr_keep,]

# 将TwosampleMR包的数据转换成MendelianRandomization的数据
MRinput_dat<-MendelianRandomization::mr_input(bx = dat$beta.exposure, 
                                              bxse = dat$se.exposure,
                                              by = dat$beta.outcome,
                                              byse = dat$se.outcome, 
                                              exposure = dat$exposure[1],
                                              outcome = dat$outcome[1], 
                                              snps = dat$SNP,
                                              effect_allele = dat$effect_allele.exposure, 
                                              other_allele = dat$other_allele.exposure,
                                              eaf = dat$eaf.exposure)


# 开跑
res_MR<-mr_conmix(MRinput_dat,psi = 0, CIMin = NA, CIMax = NA, CIStep = 0.01, alpha = 0.05)

res_TSMR<-mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))

# 结果转化函数
resConMix_to_TSMR<-function(res_mr,res_tsmr){
  res_tsmr[nrow(res_tsmr)+1,]<-NA
  res_tsmr[nrow(res_tsmr),c("id.exposure","id.outcome","outcome","exposure")]<-res_tsmr[1,c("id.exposure","id.outcome","outcome","exposure")]
  res_tsmr[nrow(res_tsmr),c("method")]<-class(res_mr)
  res_tsmr[nrow(res_tsmr),c("nsnp")]<-res_mr@SNPs
  res_tsmr[nrow(res_tsmr),c("b")]<-res_mr@Estimate
  res_tsmr[nrow(res_tsmr),c("se")]<-res_mr@Psi
  res_tsmr[nrow(res_tsmr),c("pval")]<-res_mr@Pvalue
  return(res_tsmr)
}

# 结果整合
res_total<-resConMix_to_TSMR(res_MR,res_TSMR)

# 散点图
scatter_plot<-mr_scatter_plot(res_total,dat)[[1]]+
  scale_color_lancet()+
  scale_fill_lancet()+
  theme_bw()
# 更改散点图的点的颜色
scatter_plot[["layers"]][[3]][["aes_params"]]$colour<-"black"
# 更改散点图的点的透明度
scatter_plot[["layers"]][[3]][["aes_params"]]$alpha<-0.5
# 更改散点图的横坐标名字
scatter_plot[["labels"]][["x"]]<-"SNP effect on hypothyroidism"
# 更改散点图的横坐标名字
scatter_plot[["labels"]][["y"]]<-"SNP effect on POCS"
scatter_plot
dev.off()#清空画板