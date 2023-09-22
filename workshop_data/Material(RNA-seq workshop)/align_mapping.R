library(Rsubread)
library(magrittr)
library(tidyverse)
setwd('/hdd2/workshop/mouse_align_mapping/')
# fastq files can be downloaded via SRA toolkit
# fasterq-dump SRR1552444

# We use processed fastq files containing only 1000 reads.


# Build index for alignment. Use chromosome 1 of latest version mm39.
# The whole reference can be downloaded at
# human : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.fna.gz 
# or at
# mouse : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz

Rsubread::buildindex(basename="mm39_chr1",
                      reference = "chr1.fna.gz") 


example.fastq=list.files(pattern='fastq.gz$') #Take all fastq files in the working directory.
align(index="mm39_chr1",readfile1=example.fastq,phredOffset = 33) #Align selected fastq files with reference index. Output is in bam file format.

bam.files <- list.files( pattern = ".BAM$", full.names = TRUE)
# qs <- qualityScores(all.fastq,offset=33)

props <- propmapped(files=bam.files) # check mapping ratio. In general, mapped ratio is over 80%.
fc <- featureCounts(bam.files, annot.inbuilt="mm39") #Generate gene level count matrix from bam files. Default Gene_ids are entrez_ids.
count.matrix<-fc$counts
# write.table(count.matrix,file=paste0('/custom_path/','count.txt'))

library(biomaRt)
gene_convert<-rownames(count.matrix)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset("mmusculus_gene_ensembl", mart)
# listAttributes(mart)
mrna_attributes <- getBM(mart = mart,
                         # attributes = c(platform,
                         attributes = c('entrezgene_id',
                                        'gene_biotype',
                                        'external_gene_name'),
                         filter = 'entrezgene_id',
                         values = rownames(count.matrix),
                         uniqueRows = TRUE)
mrna_attributes%<>%dplyr::filter(gene_biotype=='protein_coding')
mrna_attributes%<>%dplyr::filter(external_gene_name!='')
count.matrix<-count.matrix[rownames(count.matrix)%in%mrna_attributes[['entrezgene_id']],]
count.matrix<-count.matrix[as.character(mrna_attributes$entrezgene_id),]
GeneName=mrna_attributes$external_gene_name
count.matrix<-cbind(GeneName,count.matrix)
rownames(count.matrix)=NULL
write.table(count.matrix,file=paste0('/custom_path/','count.matrix.txt'))