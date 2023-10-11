### Enrichment analysis for a given list of genes
BiocManager::install("enrichR")
library(enrichR)
setEnrichrSite("Enrichr")
websiteLive = TRUE
listEnrichrDbs() # Available DBs within EnrichR

#GO Biological Processes
db_GO = c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
input = DEG
#input = c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr")
input=as.vector(as.list(read.delim("DEG.txt", header=F, sep = "\t"))) #check!!

res = enrichr(input, db_GO)
res[[1]]
res[[3]]
write.table(res[[3]],"erGO_BP.txt",sep="\t")


# Plot results for Gene Ontology
if (websiteLive) plotEnrich(res[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_Molecular_Function")
plotEnrich(res[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_Cellular_Component")
plotEnrich(res[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "GO_Biological_Process")

# Pathway DBs
db_pathway = c("KEGG_2019_Human", "WikiPathways_2019_Human", "Reactome_2016")
#input = c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr")
res = enrichr(input, db_pathway)
plotEnrich(res[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "KEGG 2019")
plotEnrich(res[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "WikiPathway 2019")
plotEnrich(res[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "Reactome 2016")

# miRNAs
db_mir = c("miRTarBase_2017", "TargetScan_microRNA_2017")
#input = c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr")
res = enrichr(input, db_mir)
plotEnrich(res[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "miRTarBase")
plotEnrich(res[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "TargetScan")


### Enrichment analysis for all genes without cutoff
# "result.table" consist of rownames(genes) and three columns(pvalue, adjpvalue and logFC)

BiocManager::install("fgsea") # Function or pathway enrichment analysis for all genes without using a cutoff
BiocManager::install("qusage") 
BiocManager::install("tidyverse") 
library(fgsea)
library(qusage)
library(magrittr)
library(tidyverse)

dbl = qusage::read.gmt('c2.cp.kegg.v2023.1.Hs.symbols.gmt') # list of C2 pathways(Curated genesets) from MSigDB
prerank = lfc2

fgseaRes = fgsea(dbl, prerank, nPermSimple=100000, scoreType='std', minSize=5, maxSize=500, nproc=8,eps=1e-16)
head(fgseaRes)
#leadingEdge: gene subset that most contributed the pathway pattern

#Editing leadingEdge
fgseaRes$leadingEdge = sapply(seq_len(nrow(fgseaRes)), FUN=function(x)paste0(fgseaRes$leadingEdge[x][[1]], collapse = ', '))
fgseaRes %<>% as.data.frame()
head(fgseaRes)

# Save gsea result table
#write.table(fgseaRes,file=paste0(custom_dir,'gsea_table.txt'))
write.table(fgseaRes,"gsea_table.txt", sep="\t")

# Plot Random walk(enrichment plot)
pathway_name = fgseaRes[order(fgseaRes$pval)[12],]$pathway #select most significant pathway
fgsea::plotEnrichment(dbl[[pathway_name]],
                      prerank) + labs(title=pathway_name) # plot enrichment
head(fgseaRes[order(fgseaRes$pval)[1],1:6]) #check the enrichment score 


