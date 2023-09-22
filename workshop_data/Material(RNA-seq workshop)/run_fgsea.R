# GSE229705 results
# "DESeq2.result.table" consist of rownames(genes) and three columns(pvalue, adjpvalue and logFC)
library(qusage)

dbl<-qusage::read.gmt('/hdd2/workshop/c2.all.v2023.1.Hs.symbols.gmt') # list of C2 pathways(curated genesets) from msigdb
# dbl<-qusage::read.gmt('/hdd2/workshop/c5.all.v2023.1.Hs.symbols.gmt') # list of C5 pathways(Ontology genesets) from msigdb
dbl_name=names(dbl)
dbl<-lapply(dbl_name,FUN=function(x){toupper(dbl[[x]])}) #Convert lowercase alphabet in genenames to uppercase.
names(dbl)=dbl_name


#If you want to analyze specific pathways(e.g. wikipathways or reactome)
{
  dbl_name_sub=dbl_name[str_detect(pattern='^WP_',string = dbl_name)] #select pathways start with 'WP' (wikipathways)
  dbl_name_sub=dbl_name[str_detect(pattern='^REACTOME_',string = dbl_name)] #select pathways start with 'REACTOME' (reactome)
  
  dbl_sub<-lapply(dbl_name_sub,FUN=function(x){toupper(dbl[[x]])}) #Convert lowercase alphabet in genenames to uppercase.
  names(dbl_sub)=dbl_name
}

rank_stat<-DESeq2.result.table$logFC #Take logFC to assign ranks to all genes for running fgsea
names(rank_stat)=rownames(DESeq2.result.table)%>%toupper() #fgsea takes named(gene-name) vector of rank statistics.
#If you want to analyze specific pathways, use dbl_sub, instead of dbl
fgseaRes <- fgsea(dbl, rank_stat, nPermSimple=100000,scoreType='std', minSize=5, maxSize=500, nproc=8,eps=1e-16)
fgseaRes$leadingEdge=sapply(seq_len(nrow(fgseaRes)),FUN=function(x)paste0(fgseaRes$leadingEdge[x][[1]],collapse = ', '))
fgseaRes%<>%as.data.frame()

# Save gsea result table
write.table(fgseaRes,file=paste0('/custom/paths/','gsea_table.txt'))


# Top enriched pathways
#Plot dot graph of top and bottom 10 enriched pathways
fgseaRes.sub<-fgseaRes%>%filter(pathway%in%c(fgseaRes[fgseaRes$NES<0,]$pathway[1:10],fgseaRes[fgseaRes$NES>0,]$pathway[1:10])) #select top 10 and bottom 10 enriched pathways and save these subset of fgsea result table.
deseq2_deg<-rownames(DESeq2.result.table)[DESeq2.result.table$adjpvalue<0.1]
deseq2_deg<-toupper(deseq2_deg)
fgseaRes.sub$ratio=sapply(fgseaRes.sub$pathway,FUN=function(x){sum(deg%in%dbl[[x]])/fgseaRes.sub$size[fgseaRes.sub$pathway==x]})%>%as.vector()
fgseaRes.sub$count=sapply(fgseaRes.sub$pathway,FUN=function(x){sum(deg%in%dbl[[x]])})%>%as.vector()

b<-ggplot(fgseaRes.sub,aes(x=reorder(pathway,NES),y=NES, color=log10(padj), size=count))+
  geom_point(stat = 'identity') +
  geom_hline(yintercept = 0,linetype='dashed', color='black',linewidth=0.4)+
  coord_flip() +
  scale_color_gradient(low='#F54029',high='#2e65a8',
                       limits=c(-13,0.05),
                       breaks=c(0,-3,-6,-9,-12))+
 
  labs(x="Pathway", y="NES",
       title="Top and bottom 10 enriched pathways from GSEA",
  ) +
  # ylim(-3,3)+
  theme(text = element_text(size = 15, color = 'black'),
        axis.text = element_text(size=12, color = 'black'),
        panel.background = element_rect(fill='white',colour = 'black',linewidth=2),
        panel.grid.major =element_line(color = "lightgray"),
        panel.grid.minor =element_line(color = "lightgray"))
print(b)

ggsave(paste0('/custom/paths/','20 enriched gsea_pathways.pdf'),width = 18,height = 10)





# HEATMAP
#select the most enriched pathway and generate expression heatmap to check enrichment.
path=fgseaRes.sub$pathway[1]
dgeObj <- DGEList(BD_subset)
## Perform TMM normalisation
dgeObj <- calcNormFactors(dgeObj)
logCPM <- cpm(dgeObj, log=TRUE, prior.count=3)

logCPM_z<-t(scale(t(logCPM))) #scale logCPM
rownames(logCPM_z)%<>%toupper()
logCPM_z<-logCPM_z[rownames(logCPM_z)%in%dbl[[path]],] #leave the pathway genes

gene_ranks<-rank(rank_stat)
names(gene_ranks)=names(rank_stat)
logCPM_z<-logCPM_z[granks[rownames(logCPM_z.2)]%>%sort()%>%names(),] #Reorder the genes with logFC values.

library(pheatmap)

height=3+nrow(logCPM_z.2)*0.55
width=8
pdf(paste0('/custom/paths/',path,'heatmap.pdf'),width = width,height=height)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
# pheatmap(t(logCPM_z.2[,c(1:3,5:8)]), main = paths)
main_tt=paste0('Expression heatmap of ',path, ' genes')

#control 컬럼과 test 컬럼을 정해주어야합니다.
pheatmap(logCPM_z[,c(control_columns,test_columns)],scale ='row' , main = main_tt,cluster_rows = F,fontsize = 20, fontsize_row = 18, fontsize_col = 18,
         color=rev(cols), breaks = c(min(logCPM_z.2)*(20:1)/20,0,max(logCPM_z.2)*(1:20)/20),border_color = NA)
setHook("grid.newpage", NULL, "replace")
grid.text("Sample",x=0.35, y=-0.03, gp=gpar(fontsize=24))
grid.text("Gene", x=0.9, rot=270, gp=gpar(fontsize=24))
dev.off()



