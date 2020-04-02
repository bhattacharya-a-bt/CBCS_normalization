source('nanostring_RUV_functions.R')

## dir is the directory with raw RCC files from Sabry et al 2019 (GSE130286)
setwd(dir)


require(NanoStringQCPro)
require(ggplot2)
require(EnvStats)
files.RCC = list.files()
files.RCC = files.RCC[grepl('RCC',files.RCC)]
ng = nrow(readRcc(files.RCC[1])$Code_Summary)
ncol = length(files.RCC)

raw_expression = as.data.frame(matrix(nrow = ng,ncol = length(files.RCC)+2))
colnames(raw_expression)[1:2] = c('Gene','Class')
pData = as.data.frame(matrix(nrow = length(files.RCC),ncol = 11))
colnames(pData) = c('BCAC_ID','SampleID','Owner','Comments','Date','GeneRLF','SystemAPF','imagingQC',
                    'bindingDensityQC','limitOfDetectionQC','positiveLinearityQC')
raw_expression[,1:2] = readRcc(files.RCC[1])$Code_Summary[,c(2,1)]

for (i in 1:length(files.RCC)){
  
  print(i)
  rcc = readRcc(files.RCC[i])
  raw = rcc$Code_Summary
  
  raw_expression[,i+2] = as.numeric(raw$Count)
  colnames(raw_expression)[i+2] = strsplit(files.RCC[i],'_')[[1]][1]
  pData[i,2:7] = as.vector(rcc$Sample_Attributes)
  pData$imagingQC[i] = imagingQC(rcc)
  pData$bindingDensityQC[i] = bindingDensityQC(rcc,.05,2.25)
  pData$limitOfDetectionQC[i] = limitOfDetectionQC(rcc)
  pData$positiveLinearityQC[i] = positiveLinQC(rcc)
  
}

pData$Group = as.factor(rep(c('NK_alone','CTV_6hrs',
                              'CTV_16hrs','IL2'),each = 3))

raw = raw_expression[,-c(1:2)]
fData = raw_expression[,c(1:2)]
rownames(raw) = fData$Gene
cIdx <- fData$Gene[fData$Class == "Housekeeping"]
pData$HK_Gene_Miss = colSums(raw[cIdx,] == 0)
rownames(fData) = fData$Gene
rownames(raw) = fData$Gene
rownames(pData) = colnames(raw)

#### CHECK IF HK ARE ASSOCIATED WITH PRIMARY PHENO
hk_raw = raw[cIdx,]
pval = vector(length = nrow(hk_raw))

require(MASS)

for (i in 1:nrow(hk_raw)){
  
  reg = glm.nb(as.numeric(hk_raw[i,]) ~ as.factor(pData$Group))
  pval[i] = coef(summary(reg))[2,4]
  
}

sum(pval <= .05)


### grid of datasets
library(RUVSeq)
library(DESeq2)
library(limma)
library(matrixStats)

k = 1
vsd = RUV_total(raw,pData,fData,k = k)$vsd
set = RUV_total(raw,pData,fData,k = k)$set
save(vsd,file = paste0("deangelo_ruv_vsd_k_",k,".rda"))
save(set,file = paste0("deangelo_ruv_set_k_",k,".rda"))


i = 1

results = paste0("deangelo_deg_k",i,".csv")
load(paste0("deangelo_ruv_set_k_",i,".rda"))
dds <- DESeqDataSetFromMatrix(countData = counts(set)[1:652,],
                              colData = pData(set),
                              design = ~ W_1 + Group)
dds <- DESeq(dds)
res_il2 <- as.data.frame(results(dds,contrast = c('Group','IL2','NK_alone')))
res_il2$Treatment = 'IL-2'
res6 <- as.data.frame(results(dds,contrast = c('Group','CTV_6hrs','NK_alone')))
res6$Treatment = 'CTV-1'
res16 <- as.data.frame(results(dds,contrast = c('Group','CTV_16hrs','NK_alone')))
res16$Treatment = 'CTV-1'




#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)
library(RUVSeq)
library(DESeq2)
library(limma)
library(matrixStats)
require(NanoStringNorm)
# load series and platform data from GEO

library(Biobase)
library(GEOquery)

dds.nsolver <- DESeqDataSetFromMatrix(countData = counts(set)[1:652,],
                              colData = pData(set),
                              design = ~Group)
pos = raw[grepl('POS',rownames(raw)),]
hk = raw[grepl('Housekeeping',raw_expression$Class),]
pos.factors = mean(as.numeric(apply(pos,2,geoMean)))/as.numeric(apply(pos,2,geoMean))
hk.factors = mean(as.numeric(apply(hk,2,geoMean)))/as.numeric(apply(hk,2,geoMean))
sizeFactors(dds.nsolver) <- pos.factors * hk.factors

dds.nsolver <- DESeq(dds.nsolver)
res_il2.nsolver <- as.data.frame(results(dds.nsolver,contrast = c('Group','IL2','NK_alone')))
res_il2.nsolver$Treatment = 'IL-2'
res6.nsolver <- as.data.frame(results(dds.nsolver,contrast = c('Group','CTV_6hrs','NK_alone')))
res6.nsolver$Treatment = 'CTV-1'
res16.nsolver <- as.data.frame(results(dds.nsolver,contrast = c('Group','CTV_16hrs','NK_alone')))
res16.nsolver$Treatment = 'CTV-1'

pval.plot.il2 = data.frame(Method = rep(c('nSolver','RUVSeq'),each = nrow(res_il2)),
                           P = c(res_il2.nsolver$pvalue,res_il2$pvalue))
pval.plot.il2$Contrast = 'IL-2'

pval.plot.ctv6 = data.frame(Method = rep(c('nSolver','RUVSeq'),each = nrow(res_il2)),
                            P = c(res6.nsolver$pvalue,res6$pvalue))
pval.plot.ctv6$Contrast = 'CTV-1'

pval.plot = rbind(pval.plot.il2,pval.plot.ctv6)
pval.plot$Contrast = as.factor(pval.plot$Contrast)
pvals = ggplot(data = pval.plot,
               aes(x = P,color = Method,
                   fill = Method)) + geom_histogram(alpha = .2) +
  facet_wrap(~Contrast) + theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=24),
        plot.title = element_text(size = 30),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        strip.text = element_text(size=24),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = "grey", 
                                    fill = NA, size = .1),
        legend.position = 'bottom') +
  scale_color_manual(values = c("#377EB8","#E41A1C","#4DAF4A")) + 
  scale_fill_manual(values = c("#377EB8","#E41A1C","#4DAF4A")) +
  ylab('Count')

#### IL2
il2.ruv = subset(res_il2,padj < 0.05)
il2.nsolver = subset(res_il2.nsolver,padj < 0.05)
nrow(il2.ruv)
nrow(il2.nsolver)
length(intersect(rownames(il2.ruv),rownames(il2.nsolver)))

il2.both = intersect(rownames(il2.ruv),rownames(il2.nsolver))
only.il2.ruv = rownames(il2.ruv)[!(rownames(il2.ruv) %in% il2.both)]
only.il2.nsolver = rownames(il2.nsolver)[!(row.names(il2.nsolver) %in% il2.both)]

boxplot.df.il2 = data.frame(Method = c(rep('Both',length(il2.both)),
                                       rep('RUVSeq',length(only.il2.ruv)),
                                       rep('nSolver',length(only.il2.nsolver))),
                            Median = c(rowMedians(as.matrix(log2(raw[rownames(raw) %in% il2.both,]+1))),
                                       rowMedians(as.matrix(log2(raw[rownames(raw) %in% only.il2.ruv,]+1))),
                                       rowMedians(as.matrix(log2(raw[rownames(raw) %in% only.il2.nsolver,]+1)))))
boxplot.df.il2$Contrast = 'IL2'

#### CTV6
ctv6.ruv = subset(res6,padj < 0.05)
ctv6.nsolver = subset(res6.nsolver,padj < 0.05)
nrow(ctv6.ruv)
nrow(ctv6.nsolver)
length(intersect(rownames(ctv6.ruv),rownames(ctv6.nsolver)))


ctv.both = intersect(rownames(ctv6.ruv),rownames(ctv6.nsolver))
only.ctv.ruv = rownames(ctv6.ruv)[!(rownames(ctv6.ruv) %in% il2.both)]
only.ctv.nsolver = rownames(ctv6.nsolver)[!(row.names(ctv6.nsolver) %in% il2.both)]

colnames(raw_expression)[1:2] = c('Name','Code.Class')
NanoString.mRNA.norm <- NanoStringNorm(
  x = raw_expression,
  anno = NA,
  CodeCount = 'sum',
  Background = 'none',
  SampleContent = 'housekeeping.sum',
  round.values = FALSE,
  take.log = TRUE,
  return.matrix.of.endogenous.probes = TRUE
)

boxplot.df.ctv = data.frame(Method = c(rep('Both',length(ctv.both)),
                                       rep('RUVSeq',length(only.ctv.ruv)),
                                       rep('nSolver',length(only.ctv.nsolver))),
                            Median = c(rowMedians(as.matrix(log2(raw[rownames(raw) %in% ctv.both,]+1))),
                                       rowMedians(as.matrix(log2(raw[rownames(raw) %in% only.ctv.ruv,]+1))),
                                       rowMedians(as.matrix(log2(raw[rownames(raw) %in% only.ctv.nsolver,]+1)))))
boxplot.df.ctv$Contrast = 'CTV-1'

boxplot.df = rbind(boxplot.df.ctv,boxplot.df.il2)
boxplot.df$Contrast = as.factor(boxplot.df$Contrast)
boxplot = ggplot(data = boxplot.df,
                 aes(x = Method,
                     y = Median)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=24),
        plot.title = element_text(size = 30),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        strip.text = element_text(size=24),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = "grey", 
                                    fill = NA, size = .1),
        legend.position = 'bottom') + ylab('Median raw log-counts') + facet_wrap(~Contrast)




res_il2.plot = subset(res_il2,rownames(res_il2) %in% c(il2.both,rownames(il2.nsolver),
                                                       rownames(il2.ruv)))
res_il2.nsolver.plot = subset(res_il2.nsolver,
                              rownames(res_il2) %in% 
                                c(il2.both,rownames(il2.nsolver),rownames(il2.ruv)))
res_il2.plot.tot = data.frame(Gene = rownames(res_il2.plot),
                              nSolverFC = res_il2.nsolver.plot$log2FoldChange,
                              nSolverSE = res_il2.nsolver.plot$lfcSE,
                              RUVFC = res_il2.plot$log2FoldChange,
                              RUVSE = res_il2.plot$lfcSE)
res_il2.plot.tot$Method = ifelse(res_il2.plot.tot$Gene %in% il2.both,'Both',
                                 ifelse(res_il2.plot.tot$Gene %in% rownames(il2.nsolver),
                                        'nSolver only','RUVSeq only'))
res_il2.plot.tot$Contrast = 'IL-2'

res_ctv.plot = subset(res6,rownames(res6) %in% c(ctv.both,only.ctv.nsolver,
                                                       only.ctv.ruv))
res_ctv.nsolver.plot = subset(res6.nsolver,
                              rownames(res6) %in% 
                                c(ctv.both,only.ctv.nsolver,
                                  only.ctv.ruv))
res_ctv.plot.tot = data.frame(Gene = rownames(res_ctv.nsolver.plot),
                              nSolverFC = res_ctv.nsolver.plot$log2FoldChange,
                              nSolverSE = res_ctv.nsolver.plot$lfcSE,
                              RUVFC = res_ctv.plot$log2FoldChange,
                              RUVSE = res_ctv.plot$lfcSE)
res_ctv.plot.tot$Method = ifelse(res_ctv.plot.tot$Gene %in% ctv.both,'Both',
                                 ifelse(res_ctv.plot.tot$Gene %in% only.ctv.nsolver,
                                        'nSolver only','RUVSeq only'))
res_ctv.plot.tot$Contrast = 'CTV-1'
res.plot.tot = rbind(res_il2.plot.tot,
                     res_ctv.plot.tot)


res.plot.tot$Contrast = as.factor(res.plot.tot$Contrast)
res.plot.tot$Method = factor(res.plot.tot$Method,
                             levels = c('Both','nSolver only','RUVSeq only'))
res.plot.tot = res.plot.tot[order(res.plot.tot$Method,decreasing = T),]
fc.plot = ggplot(data = res.plot.tot,
       aes(x = RUVFC, y = nSolverFC, color = Method)) + geom_point(aes(size = RUVSE)) +
  scale_color_manual(values = c("#4DAF4A","#377EB8","#E41A1C")) + 
  geom_hline(yintercept = 0,linetype = 2) +
  geom_vline(xintercept = 0,linetype = 2) +  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=24),
        plot.title = element_text(size = 30),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        strip.text = element_text(size=24),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = "grey", 
                                    fill = NA, size = .1),
        legend.position = 'bottom') + facet_wrap(~Contrast) +
  xlab('Log-fold change (RUVSeq)') + ylab('Log-fold change (nSolver)') +
  labs(size = 'SE (RUVSeq)') + geom_abline(slope = 1,intercept = 0,linetype = 2, color='grey') +
  guides(colour = guide_legend(override.aes = list(size=4)))
