

#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

gene="ARL15"               #目标基因的名称
corFilter=0.3            #相关系数的过滤条件
pFilter=0.05             #相关性检验pvalue的过滤条件
expFile="combined_RNAseq_counts.txt"      #表达数据文件
drugFile="DrugPredictions.csv"            #药物敏感文件

#读取输入文件，并对输入文件进行整理
rt0=read.csv(drugFile,header=T,check.names=F,row.names = 1)
rt0=as.matrix(t(rt0))
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=log2(data+1)

#提取目标基因表达量
x=as.numeric(data[gene,])
outTab=data.frame()
#对基因进行循环，进行相关性检验
for(j in rownames(rt0)){
	if(gene==j){next}
    y=as.numeric(rt0[j,])
	corT=cor.test(x, y, method = 'pearson')
	cor=corT$estimate
	pvalue=corT$p.value
	outTab=rbind(outTab, cbind(Query=gene, Gene=j, cor, pvalue))
	#保存满足条件的基因
	if((abs(cor)>corFilter) & (pvalue<pFilter)){
		#可视化
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
			xlab(paste0(gene, " expression"))+ ylab(paste0(j, "   drug sensitivity"))+
			geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
			stat_cor(method = 'pearson', aes(x =x, y =y))
		pdf(file=paste0(j, ".pdf"), width=5, height=4.6)
		print(p1)
		dev.off()
	}
}



