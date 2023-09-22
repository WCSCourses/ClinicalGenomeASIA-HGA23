#Bioconductor(https://www.bioconductor.org/) provides tools for the analysis and comprehension of high-throughput genomic data. 
#Bioconductor uses the R statistical programming language, and is open source and open development.
#BiocManager package is used to install and manage packages from the Bioconductor project.

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GEOquery") #How to install an R package from Bioconductor using BiocManager.
BiocManager::install("limma")
BiocManager::install("edgeR")

library(GEOquery) #How to activate an installed package.
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 3) #Increases the capacity of your computer.
#The number 131072 may change depending on your system

gse = getGEO("GSE80336",GSEMatrix=TRUE,AnnotGPL=TRUE) #Use getGEO function to download a microarray gene expression data from GEO database.
# GSE45827: In a cohort study of primary invasive breast cancer 
# (41 Basal, 30 HER2, 29 Luminal A and 30 Luminal B) as well as 11 normal tissues samples and 14 cell lines.

gse = gse[[1]] #Try 'gse' in the command line to see what are in it.
eset = exprs(gse) # Extract the matrix of gene expression data  
dim(eset)
print(eset[1:5,1:5]) #Try to see how 'eset' looks like
length(which(eset<0)) # Checking log normalization

#Extracting feature (gene) description data
fset = fData(gse) 
colnames(fset) #Find "Symbol", "Gene Symbol" or "Gene symbol"
ID = fset$'Gene Symbol' #Double quotation also works
rownames(eset) = ID #Each row of eset is named with 'gene symbols'

#Extracting phenotype data
pset = pData(gse) 
colnames(pset) #See what are included
pset$title
age=as.numeric(pset$`age (years):ch1`)
sex=pset$`Sex:ch1`

BiocManager::install("edgeR") #You may have to remove a file or folder, C:\Users\User\Documents\R\win-library\4.0/00LOCK

library(GEOquery)
library(edgeR)

sfiles = getGEOSuppFiles('GSE80336') #bipolar disease RNA-seq data
fnames = rownames(sfiles)
# there is only one supplemental file
BD = read.delim(fnames[1],header=TRUE)
head(BD)

#Take only expression values and gene symbols from the matrix
BD1=BD[,5:40]
dim(BD1)
rownames(BD1)=BD$GeneSymbol


#Take the index of Control and Bioplolar disorder cases
control=grep('C',colnames(BD1))
bd=grep('BD',colnames(BD1))
fac=factor(c(rep('C',length(control)), rep('BD',length(BioPolar))))

#Limma-Voom analysis: 
#Read count data are converted to continuous data using voom transform which is input to limma for DE analysis

#calculate normalization factor
norm.factor = edgeR::calcNormFactors(BD1, method = 'TMM') # calculating normalization factor
design = model.matrix(~fac+age)
voom.data = limma::voom(BD1, design = design, lib.size = colSums(BD1) * norm.factor) #voom transformation of read count data
voom.data$genes = BD$GeneSymbol
dim(voom.data$E) # This transformed expression data can be treated like microarray data

voom.fitlimma = limma::lmFit(voom.data, design = design) # Perform differential expression analysis
voom.fitbayes = limma::eBayes(voom.fitlimma)

#Take results
voom.pvalues = voom.fitbayes$p.value[, 2]
voom.adjpvalues = p.adjust(voom.pvalues, method = 'BH')
voom.logFC = voom.fitbayes$coefficients[, 2]
#voom.score = 1 - voom.pvalues
voom.result.table = data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC)
rownames(voom.result.table) = BD$GeneSymbol
dim(voom.result.table)
head(voom.result.table)


### DESeq2 analysis
#Download BD1 data as presented above
BiocManager::install("DESeq2")
library(DESeq2)

DESeq2.ds =DESeq2::DESeqDataSetFromMatrix(countData = BD1, colData = data.frame(condition = fac,age=age), design = ~ condition+age)
DESeq2.ds = DESeq2::DESeq(DESeq2.ds, fitType = "parametric", test = "Wald", betaPrior = TRUE )

impute.outliers=TRUE
if(impute.outliers==TRUE){
  DESeq2.ds.clean = DESeq2::replaceOutliersWithTrimmedMean(DESeq2.ds)
  DESeq2.ds.clean = DESeq2::DESeq(DESeq2.ds.clean, fitType = "parametric", test = "Wald", betaPrior = TRUE)
  DESeq2.ds = DESeq2.ds.clean
}

DESeq2.results = DESeq2::results(DESeq2.ds, independentFiltering = TRUE, cooksCutoff = TRUE) #results (p.adj, p_value, log2foldchange) from genes with all zero counts become NA. 

#Assign gene names to result
rownames(DESeq2.results)=BD$GeneSymbol

#filter out genes with NA p values
DESeq2.results=DESeq2.results[-which(is.na(DESeq2.results$padj)==TRUE),]


#There are 6 columns. 'baseMean', 'log2FoldChange', 'lfcSe', 'stat', 'pvalue' and 'padj'.
DESeq2.pvalues = DESeq2.results$pvalue
DESeq2.adjpvalues = DESeq2.results$padj
DESeq2.logFC = DESeq2.results$log2FoldChange
#DESeq2.score = 1 - DESeq2.pvalues
DESeq2.result.table = data.frame('pvalue' = DESeq2.pvalues, 'adjpvalue' = DESeq2.adjpvalues, 'logFC' = DESeq2.logFC)
rownames(DESeq2.result.table) = rownames(DESeq2.results)
head(DESeq2.results)
dim(DESeq2.results)
sig.genes=rownames(BD1)[DESeq2.DE_index] #This gene list can be used for enrichr analysis

#compare DE genes from voom and DESeq2.
DESeq2.DE_index=which(DESeq2.result.table$adjpvalue<0.1)
length(DESeq2.DE_index)
Voom.DE_index=which(voom.result.table$adjpvalue<0.1)
length(Voom.DE_index) #limma is more conservative
#Shared DE genes
Shared_DEgenes=intersect(rownames(DESeq2.result.table)[DESeq2.DE_index],rownames(voom.result.table)[Voom.DE_index])
length(Shared_DEgenes)

  