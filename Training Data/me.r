options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
GDSC2_Expr = readRDS(file=file.path('GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path("GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 
load("luda-limma-8-22-75mut.rda")
exp<-  log2(df + 1) 
p7 <- read.table("p7.txt") 
testExpr <- exp[rownames(exp)%in%p7$V1,]
testExpr = as.matrix(testExpr)
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

res <- fread("calcPhenotype_Output\\DrugPredictions.csv")
res <- as.data.frame(res)
rownames(res) <- res$V1
res <- res[,-1]
group <- c(rep("Wild", 429), rep("Mut",75)) %>% factor(., levels = c("Wild", "Mut"), ordered = F)

# 创建一个空的列表来存储所有的箱线图
boxplots_list <- list()

# 循环遍历198列
for (i in 1:198) {
  # 提取当前列的数据并将其存储在 Cam 中
  Cam <- data.frame(senstivity = res[, i], Type = group)
  
  # 执行统计检验，例如 t.test
  result <- t.test(senstivity ~ Type, data = Cam)
  
  # 创建箱线图
  boxplot <- ggboxplot(Cam, x = "Type", y = "senstivity", fill = "Type",
                       xlab = "Type",
                       ylab = paste0(colnames(res)[i], " senstivity (IC50)"),
                       legend.title = "Type",
                       palette = c("green", "red")) +
    # 根据显著性水平为每个箱线图添加星号标签
    annotate("text", x = 1.5, y = max(Cam$senstivity), 
             label = ifelse(result$p.value < 0.001, "***", 
                            ifelse(result$p.value < 0.01, "**", 
                                   ifelse(result$p.value < 0.05, "*", "ns"))))
  
  # 将箱线图添加到列表中
  boxplots_list[[i]] <- boxplot
}

# 将所有的箱线图保存为 PDF 文件
pdf("all_boxplots17.pdf", width = 8, height = 6)
for (i in 1:198) {
  print(boxplots_list[[i]])
}
dev.off()
