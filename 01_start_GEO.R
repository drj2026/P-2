rm(list = ls())
#打破下载时间的限制,改前60秒，改后10w秒
options(timeout = 100000) 
options(scipen = 20)#不要以科学计数法表示

#传统下载方式
library(GEOquery)#下载读取和提取GEO数据信息的一个包
eSet = getGEO("GSE13911", destdir = '.', getGPL = F)#可以实现下载并读取文件
#网速太慢，下不下来怎么办
#1.从网页上下载/发链接让别人帮忙下，放在工作目录里
#2.试试geoChina,只能下载2019年前的表达芯片数据
#library(AnnoProbe)#使用的是技能树服务器
#eSet = geoChina("GSE7305") #选择性代替第8行
eSet <- getGEO(filename = "GSE13911_series_matrix.txt.gz", getGPL = FALSE)
#研究一下这个eSet
class(eSet)
dim(eSet)
exprSet <- exprs(eSet)
dim(exprSet)          # 行 = 探针/基因，列 = 样本
exprSet[1:5, 1:5]     # 查看前5行5列
length(eSet)

eSet = eSet[[1]] 
class(eSet)

#(1)提取表达矩阵exp
exp <- exprs(eSet)
#⭐第一个要检查的地方👇，表达矩阵行列数，正常是几万行，列数=样本数，
#如果0行说明不是表达芯片或者是遇到特殊情况，不能用此流程分析
dim(exp)
#⭐二个要检查的地方👇
range(exp)#看数据范围决定是否需要log，是否有负值，异常值，如有负值，结合箱线图进一步判断
#⭐可能要修改的地方👇查看geo网站gsm，看是否去了log几
exp = log2(exp+1) #需要log才log，不需要log要注释掉这一句。+1是避免有0和<1值带来的影响
#⭐第三个要检查的地方👇
boxplot(exp,las = 2) #看是否有异常样本
#第一个办法：删掉异常样本exp[,-第几列]
exp = limma::normalizeBetweenArrays(exp)
#样本数较大可以随机选择几个样本
#sample(1:100,10)表示1到100样本随机取10个
boxplot(exp[,sample(1:130,20)],las = 2)
#(2)提取临床信息
pd <- pData(eSet)
#⭐多分组中提取两分组的代码示例，二分组不需要
#if(F){
  #因为现在这个例子不是多分组，所以编造一列做示例。
  pd$fake = paste0(rep(c("a","b","c","d"),each = 5),1:5)
  k1 = str_detect(pd$fake,"b");table(k1)
  k2 = str_detect(pd$fake,"c");table(k2)
  pd = pd[k1|k2,]
}

library(stringr)
k1 = str_detect(pd$title,"Normal");table(k1)
k2 = str_detect(pd$title,"Ovarian Cancer");table(k2)
pd = pd[k1|k2,]

#(3)让exp列名与pd的行名顺序完全一致
p = identical(rownames(pd),colnames(exp));p#identical表明内容顺序完全一致
if(!p) {
  s = intersect(rownames(pd),colnames(exp))
  exp = exp[,s]
  pd = pd[s,]
}

#(4)提取芯片平台编号，后面要根据它来找探针注释
gpl_number <- eSet@annotation;gpl_number#@是对象提取子集的符号
save(pd,exp,gpl_number,file = "step1output.Rdata")

# 原始数据处理的代码，按需学习
# https://mp.weixin.qq.com/s/0g8XkhXM3PndtPd-BUiVgw
