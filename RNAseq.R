library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(RColorBrewer)
library(gplots)
library(tidyverse)
## loading data 
seqdata <- read.csv("Counts.csv", stringsAsFactors = FALSE, row.names = 1)
sampleInfo <- read.csv("groups.csv", stringsAsFactors = FALSE)

## design matrix
condition <- factor(sampleInfo$Condition)
design <- model.matrix(~ 0 + condition)
design

## create DGElist object
countdata <- seqdata
group = factor(sampleInfo$Condition)

dge <- DGEList(counts = countdata, group = group, remove.zeros = FALSE)

##add annotaion
ENTREZID <- mapIds(org.Mm.eg.db, row.names(dge), keytype = "SYMBOL", column = "ENTREZID")
row.names(dge$counts) <- ENTREZID
ann <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(dge$counts), 
                             columns = c("ENTREZID", "SYMBOL", "GENENAME"))    
head(ann)
dge$genes <- ann

i <- is.na(dge$genes$ENTREZID)    
dge <- dge[!i, ]
dge <- dge[!is.na(rownames(dge$counts)), ]

## Filtering Low-Expressed Genes 
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.size=FALSE]
dge <- calcNormFactors(dge)

#############################################################

### Exploring the data

#library sizes
barplot(dge$samples$lib.size/1e06, names = colnames(dge), las = 2, 
        ann = FALSE, cex.names = 0.75, col = "lightskyblue", space = 0.5)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

#MD plot
par(mfrow=c(2,3))
for (i in 1:6) {
    plotMD(dge, column = i,
           xlab = "Average log CPM (All samples)",
           ylab = "log-ratio (this vs others")
    abline(h=0, col="red", lty=2, lwd=2)
}

#boxplot
par(mfrow=c(1,2))
logcounts <- cpm(dge, log = TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las = 2)
abline(h=median(logcounts), col="blue")
title("Unnormalized logCPMs")
logCountsNorm <- cpm(dge, log = TRUE)
boxplot(logCountsNorm, xlab="", ylab="Log2 counts per million", las=2)
abline(h=median(logCountsNorm), col="blue")
title("Normalized logCPMs")
#dev.off()

# multidimentional scaling (MDS) plot
pseudoCounts <- log2(dge$counts + 1)  
colConditions <- brewer.pal(3, "Set2")  
colConditions <- colConditions[match(sampleInfo$Condition,  levels(factor(sampleInfo$Condition)))]  
patients <- c(8, 15, 16)[match(sampleInfo$Patient,  levels(factor(sampleInfo$Patient)))]  
plotMDS(pseudoCounts, pch = patients, col = colConditions, xlim =  c(-2,2))  
legend("topright", lwd = 2, col = brewer.pal(3, "Set2")[1:2],  legend = levels(factor(sampleInfo$Condition)))  
legend("bottomright", pch = c(8, 15, 16),  legend = levels(factor(sampleInfo$Patient)))  
#dev.off()

# heatmap
logcountsNorm <- cpm(dge,log=TRUE)  
var_genes <- apply(logcountsNorm, 1, var)  
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]  
highly_variable_lcpm <- logcountsNorm[select_var,]  
mypalette <- brewer.pal(11,"RdYlBu")  
morecols <- colorRampPalette(mypalette)  
col.con <- c(rep("purple",3),  rep("orange",3))[factor(sampleInfo$Condition)]  
heatmap.2(highly_variable_lcpm,  col=rev(morecols(50)),trace="none",  main="Top 10 most variable genes",  ColSideColors=col.con,scale="row",  margins=c(12,8),srtCol=45)  
#dev.off() 
#############################################################

logCPM <- cpm(dge, log = TRUE, prior.count = 3)

# voom transformation
v <- voomWithQualityWeights(dge, design = design, plot = TRUE)

my.contrasts <- makeContrasts(conditionKO-conditionWT, levels = design)

# Fit linear model
vfit <- lmFit(v, design)

# apply contrast
vfit <- contrasts.fit(vfit, contrasts = my.contrasts)

vfit <- eBayes(vfit, robust = TRUE)
options(digits = 3)
topTable(vfit, sort.by = "P")
DEG <- topTable(vfit,n="all", sort.by = "P")
sum(DEG$adj.P.Val<0.05)
DEG <- dplyr::filter(DEG, adj.P.Val < 0.05)
write.csv(DEG, file = "DEG.csv")


dt <- decideTests(vfit)
summary(dt)

### GLM

fitq <- glmQLFit(dge, design)
res1 <- glmQLFTest(fitq, contrast = my.contrasts)
topTags(res1)
res100 <- topTags(res1, n=100)
res100 <- res100$table

res_all <- topTags(res1, n="all", adjust.method = "BH")
res_all <- res_all$table
res_all <- dplyr::filter(res_all, FDR < 0.05)


write.fit(vfit, adjust = "BH", file="results.txt")
res <- read.delim("results.txt",stringsAsFactors = FALSE, header = T)
res <- dplyr::filter(res, P.value.adj  < 0.05)
####################################################################
# visualizations

# heatmap
fitq <- glmQLFit(dge, design)
qlfq <- glmQLFTest(fitq, contrast = my.contrasts)
logCPM <- cpm(dge, prior.count = 2, log = TRUE)
rownames(logCPM) <- dge$genes$SYMBOL
colnames(logCPM) <- paste(dge$samples$group, 1:3, sep = "-")
o <- order(qlfq$table$PValue)
logCPM <- logCPM[o[1:100], ]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM, col = col.pan, Rowv = TRUE, scale = "none", 
          trace = "none", dendrogram = "both", cexRow = 0.6, cexCol = 1.4,
          margins = c(5,5), lhei = c(2,10), lwid = c(2,8),
          key = T,keysize = 1.5, density.info = "none")

# MD plot
plotMD(vfit, column = 1, status = dt[,1], main = colnames(vfit)[1], xlim=c(-8, 13))
#abline(h=0,col="darkgrey")

# venn diagram
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"), 
            names = c(conditionKO="KO", conditionWT="WT"),  include = c("up", "down"),
            counts.col=c("red", "blue"))

write.fit(vfit, dt, adjust = "BH", file="results.csv")

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)


#########################################################
## Ontology and pathways

## KEGG
fitq <- glmQLFit(dge, design)
qlfq <- glmTreat(fitq, contrast = my.contrasts, lfc=2)
keg <- kegga(qlfq, species = "Mm")
kegg20 <- topKEGG(keg, sort = "up", n = 20)
write.csv(kegg20, file = "kegg20.csv")

## GO
go <- goana(qlfq, species = "Mm")
topGO20 <- topGO(go, sort = "up", n=20)
write.csv(topGO20, file = "topGO20.csv")

topGO100 <- topGO(go, sort = "up", n = 100)
write.csv(topGO100, file = "topGO100.csv")
