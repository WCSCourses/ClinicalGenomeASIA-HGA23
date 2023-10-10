#Bioconductor(https://www.bioconductor.org/) is an open source repository of R packages for analysis of high-throughput genomic data. 
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

#Installing an R package from Bioconductor
BiocManager::install("GEOquery") #Retrieving gene expression data (microarray or RNA-seq read counts) from GEO database

BiocManager::install("limma") 
BiocManager::install("edgeR") 
BiocManager::install("DESeq2") # Three popularly used R packages for DE analysis of RNA-seq data

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 3) #Increasing the capacity of your computer.
#The number 131072 depends on your system

library(GEOquery) #How to attach an installed package.

#bipolar disease RNA-seq data
sfiles = getGEOSuppFiles('GSE80336') #Downloading RNA-seq read count data from a supplementary file
fnames = rownames(sfiles)
BD = read.delim(fnames, header=TRUE, sep="\t")  
head(BD)
colnames(BD)

#Downloading meta data
gse = getGEO("GSE80336",GSEMatrix=TRUE,AnnotGPL=TRUE) 
gse = gse[[1]] #Try 'gse' in the command line to see what are in it.

#Extracting feature (gene) description data. No feature data for GSE80336.
#fset = fData(gse) 
#colnames(fset) #Find "Symbol", "Gene Symbol" or "Gene symbol"
#ID = fset$'Gene Symbol' #Double quotation also works
#rownames(eset) = ID #Each row of eset is named with 'gene symbols'

#Extracting phenotype data
pset = pData(gse) 
colnames(pset) #See what are included
pset$title
age=as.numeric(pset$`age (years):ch1`)
sex=pset$`Sex:ch1`


## DE analysis of RNA-seq data using limma

library(edgeR)
library(limma)
#Take only expression values and gene symbols from the matrix
BD1=BD[,5:40]
rownames(BD1)=BD$GeneSymbol # Only the numbers are data
head(BD1)

#Taking the index of control and biplolar cases
control=grep('C',colnames(BD1))
bipolar=grep('BD',colnames(BD1))

#Defining design matrix
condition=c() #Newly define condition vector
condition[control]='C'
condition[bipolar]='BD'
condition=factor(condition)
design = model.matrix(~0+condition)
#design = model.matrix(~0+condition+sex). #Also trying covariate. Compare the difference in "voom.result.table"

norm.factor = edgeR::calcNormFactors(BD1, method = 'TMM') # TMM normalization using edgeR
voom.data = limma::voom(BD1, design = design, lib.size = colSums(BD1) * norm.factor)  #log-transformed data & precision weight factor
head(voom.data)
dim(voom.data$E)

#Selecting the pair of conditions to compare. The first term is the "test" condition.
contrast.matrix = makeContrasts(conditionBD-conditionC, levels=design)

fit = lmFit(voom.data, design = design) # Performing differential expression analysis
fit = contrasts.fit(fit,contrast.matrix) 
fit = eBayes(fit)

#Taking results
voom.pvalues = fit$p.value
voom.adjpvalues = p.adjust(voom.pvalues, method = 'BH')
voom.logFC = fit$coefficients

voom.result.table = data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC)
head(voom.result.table)
colnames(voom.result.table)=c("pvalue", "adjpvalue", "logFC")
head(voom.result.table)

lfc = voom.result.table$logFC
names(lfc)=rownames(voom.result.table)
write.table(lfc,"lfc.txt", sep="\t") # Input data for preranked GSEA

DE_index=which(voom.result.table$adjpvalue<0.25) #
length(DE_index)

#Comparing the results with cpm normalized data
BDcpm=edgeR::cpm(BD1)
head(BDcpm)
aa=data.frame("cont" = apply(BDcpm[,1:18],1,mean), "tst" = apply(BDcpm[,19:36],1,mean))
cpmLogFC=log2(aa[,2]/aa[,1]) #a rough calculation of log2FC



### DESeq2 analysis
library(DESeq2)

control=grep('C',colnames(BD1)) #sample indexes for each condition
bipolar=grep('BD',colnames(BD1))

condition=c() #Defining group labels
condition[control]="con"
condition[bipolar]="tst" #"tst">"con", higher order is recognized as test condition
condition=factor(condition)
## Important!!! DESeq2 automatically assumes letter of larger order as "test" condition. You can also use numbers.
## If condition[control]="C"; condition[bipolar]="BD" as in the limma example, then BD is considered as control condition and sign of analysis results will be reversed. 

DESeq2.ds = DESeqDataSetFromMatrix(countData = BD1, colData = data.frame(condition), design = ~condition)
#DESeq2.ds = DESeqDataSetFromMatrix(countData = BD1, colData = data.frame("condition" = condition, "sex"= factor(sex), "age" = factor(age)), design = ~age+condition)
#Also try design = ~age+condition
#Important!! Last variable must be "condition" to be compared for DESeq2
#However, the first variable must be "condition" to be compared for limma

DESeq2.ds = DESeq(DESeq2.ds, fitType = "parametric", test = "Wald", betaPrior = TRUE )
dim(DESeq2.ds)

#Truncating outlier counts and recalculating
DESeq2.ds.clean = DESeq2::replaceOutliersWithTrimmedMean(DESeq2.ds)
DESeq2.ds.clean = DESeq2::DESeq(DESeq2.ds.clean, fitType = "parametric", test = "Wald", betaPrior = TRUE)
DESeq2.ds = DESeq2.ds.clean
dim(DESeq2.ds)

DESeq2.results = DESeq2::results(DESeq2.ds, independentFiltering = TRUE, cooksCutoff = TRUE) #results (p.adj, p_value, log2foldchange) from genes with all zero counts become NA. 
head(DESeq2.results) #Shows which conditions were compared

#filter out genes with NA p values
DESeq2.results=DESeq2.results[-which(is.na(DESeq2.results$padj)==TRUE),]
head(DESeq2.results)

DE_index2=which(DESeq2.results$padj<0.1)
length(DE_index2)
DEG=rownames(DESeq2.results)[DE_index2]

write.table(DEG,"DEG.txt", row.names = F, col.names = F, quote = F) #Writing only gene names. Input for enrichR

lfc2 = DESeq2.results$log2FoldChange
names(lfc2)=rownames(DESeq2.results)
write.table(lfc2,"lfc2.txt", sep="\t") # Input data for preranked GSEA





