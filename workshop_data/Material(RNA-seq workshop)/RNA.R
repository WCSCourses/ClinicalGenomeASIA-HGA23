#Bioconductor(https://www.bioconductor.org/) provides tools for the analysis and comprehension of high-throughput genomic data. 
#Bioconductor uses the R statistical programming language, and is open source and open development.
#BiocManager package is used to install and manage packages from the Bioconductor project.

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GEOquery") #How to install an R package from Bioconductor using BiocManager.
BiocManager::install("limma")

library(GEOquery) #How to activate an installed package.
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 3) #Increases the capacity of your computer.
#The number 131072 may change depending on your system

gse = getGEO("GSE229705",GSEMatrix=TRUE,AnnotGPL=TRUE) 
#gse = getGEO("GSE45827",GSEMatrix=TRUE,AnnotGPL=TRUE) #Use getGEO function to download a microarray gene expression data from GEO database.
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
col=cbind(pset$`tissue:ch1`,pset$`recurrence:ch1`)
unique(pset$`tissue:ch1`)
unique(pset$`recurrence:ch1`)
library(magrittr)
p_subset<-pset%>%dplyr::filter(`tissue:ch1`=='Tumor')%>%
  dplyr::filter(`recurrence:ch1`%in%c("Recurrence","No Progression (>5y)"))

column_select<-p_subset$title%>%gsub(pattern=' [(]RNA-seq[)]',replacement = '')%>%gsub(pattern='-',replacement='.')

BD_subset<-BD[,column_select]
head(BD[,column_select])
condition<-factor(p_subset$`recurrence:ch1`)


column_select[!(column_select%in%colnames(BD))]
# 'No progresssion'%in%c("Recurrence","No Progression (>5y)")
group = pset$source_name_ch1 #Select one to use for sample names
colnames(eset)=group

################################################################################ 
#Differential expression analysis using limma package
################################################################################

library(limma) #Load limma package for differential expression analysis
design = model.matrix(~0 + group) #Input group information of samples
colnames(design) = c("Basal","CellLine","Her2","LuminalA", "LuminalB", "Normal")
#Rename the six subtypes of breast cancer

head(design)

fit = lmFit(eset,design) #Linear model fitting
contrast.matrix = makeContrasts(Basal-LuminalB,LuminalB-LuminalA,LuminalA-Normal, levels=design) #Select which pairs of conditions to compare
fit = contrasts.fit(fit,contrast.matrix) 
fit = eBayes(fit) #Bayesian differential expression analysis

#Pair-wise analysis results
results1 = topTable(fit, coef=1, adjust="BH", number=Inf) #Basal vs. LuminalB
results2 = topTable(fit, coef=2, adjust="BH", number=Inf) #LuminalB vs. LuminalA
results3 = topTable(fit, coef=3, adjust="BH", number=Inf) #LuminalA vs. Normal, default number = 10

head(results2)
ls(results2)

#Function/pathway analysis of the second results
DEG=which((results2$adj.P.Val<0.01)&(results2$logFC>2)) #Four-fold upregulated genes
length(DEG)
results2[DEG,1]

#Extracting gene list for function analysis
DEG=which((results2$adj.P.Val<0.01)&(results2$logFC>2))
slash=grep("///",results2[,1])
DEG1=setdiff(DEG,slash)
length(DEG1)
DEgene=unique(results2[DEG1,1])
length(DEgene)

#Save a result as CSV file
write.csv(results2[DEG,], file="lumB_lumA.csv")

results = decideTests(fit) #Save analysis results 
vennDiagram(results) #Visualize all the pair-wise anaysis results using venn diagram




############################################################################## 
#Differential expression analysis of RNA-seq data
##############################################################################

BiocManager::install("edgeR") #You may have to remove a file or folder, C:\Users\User\Documents\R\win-library\4.0/00LOCK

library(GEOquery)
library(edgeR)

sfiles = getGEOSuppFiles('GSE229705') #lung cancer
#sfiles = getGEOSuppFiles('GSE216345') 
#sfiles = getGEOSuppFiles('GSE80336') #bipolar disease RNA-seq data
fnames = rownames(sfiles)[2]
# there is only one supplemental file
BD = read.delim(fnames[1],header=TRUE, sep=",")
head(BD)
colnames(BD)
pset$title%>%gsub(pattern=' [(]RNA-seq[)]',replacement = '')%>%
  gsub(pattern='-',replacement='.')

library(tidyr)
library(magrittr)


#Take only expression values and gene symbols from the matrix
BD1=BD[,5:40]
rownames(BD1)=BD$GeneSymbol


#Take the index of Control and Bioplolar disorder cases
control=grep('C',colnames(BD1))
BioPolar=grep('BD',colnames(BD1))


#Limma-Voom analysis: 
#Read count data are converted to continuous data using voom transform which is input to limma for DE analysis

#calculate normalization factor
norm.factor = edgeR::calcNormFactors(BD1, method = 'TMM') # calculating normalization factor
design = model.matrix(~factor(c(rep('C',length(control)), rep('BD',length(BioPolar)))))
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

head(voom.result.table)



### Optional

### DESeq2 analysis
#Download BD1 data as presented above
BiocManager::install("DESeq2")
library(DESeq2)

DESeq2.ds =DESeq2::DESeqDataSetFromMatrix(countData = BD1, colData = data.frame(condition = factor(c(rep('C',length(control)),rep('BD',length(BioPolar))))), design = ~ condition)
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
DESeq2.results

sig.genes=rownames(BD1)[DESeq2.DE_index] #This gene list can be used for enrichr analysis

#compare DE genes from voom and DESeq2.
DESeq2.DE_index=which(DESeq2.result.table$adjpvalue<0.1)
length(DESeq2.DE_index)
Voom.DE_index=which(voom.result.table$adjpvalue<0.1)
length(Voom.DE_index)
#Shared DE genes
Shared_DEgenes=intersect(rownames(DESeq2.result.table)[DESeq2.DE_index],rownames(voom.result.table)[Voom.DE_index])
length(Shared_DEgenes)

  