library(Seurat)
require(GEOquery)
require(DESeq2)
require(ggraph)
#setwd('~/sciebo/Deconv_RNA/code/')
#source("deconv_vis.R")

devtools::load_all("./rnamagnet/")

## BulkRNASeq Data
data1 <- read.csv("~/sciebo/Helene/Deconv/HeleneMouse/count_matrix_megaK_vsn.csv") 
## BulkRNASeq Metadata
meta1 <- read.csv("~/sciebo/Helene/Deconv/HeleneMouse/mega_meta.csv")
## Need to normalize here and map the genes mouse to human
log2.cpm.filtered.norm.df<- log2.cpm.filtered.norm.df_1
log2.cpm.filtered.norm.df<- log2.cpm.filtered.norm.df[1:16]
colnames(log2.cpm.filtered.norm.df) <- sampleLabels
#OR DO USE THE UNFILTERED UNLOG TRANSFORMED COUNT MATRIX

cpm <- cpm(myDGEList)

cpm<- as.data.frame(cpm)
cpm$gene<- rownames(cpm)

dgList_1 <- calcNormFactors(myDGEList.filtered, method="TMM")
library(edgeR)



# Basic function to convert human to mouse gene names

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(biomaRt)
library(tidyverse)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 <- getLDS(attributes = c("hgnc_symbol"), 
                  filters = "hgnc_symbol", 
                  values = log2.cpm.filtered.norm.df_1$gene, 
                  mart = human, 
                  attributesL = c("mgi_symbol"), 
                  martL = mouse, uniqueRows = T)

CO2 <- log2.cpm.filtered.norm.df_1 %>% 
  left_join(genesV2, by = c( "gene" ='HGNC.symbol')) 

CO3<- na.omit(CO2)
CO4<- distinct(CO3,MGI.symbol, .keep_all= TRUE )

row.names(CO4)<- CO4$MGI.symbol

CO5<- CO4[,-c(17,18)]

mean_by_cluster<- read.csv(file = "mean_by_cluster.csv", sep = ",",
                           row.names = 1, header = T)

write.table(mean_by_cluster, file ="ref1.txt",sep="\t", row.names = T, col.names = NA )

# This might help: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

##reference prep.

load('../data/RNAMagnetDataBundle/NicheData10x.rda')
?FindAllMarkers
data <- NicheData10x
tmp <- FindAllMarkers(seurat_object_1, test.use = "roc", return.thresh = 0)
usegenes <- unique(tmp$gene[(tmp$myAUC > 0.6 |tmp$myAUC < 0.2) ])

rawdata<- seurat_object_1@assays$RNA@counts
norm_matrix <- RelativeCounts(rawdata, scale.factor = 1e6, verbose = TRUE)

cpm_seurat<- CreateSeuratObject(norm_matrix, project = "CPM")
cpm_seurat<- AddMetaData(cpm_seurat, metadata = seurat_object_1@meta.data)
cpm_seurat<- SetIdent(cpm_seurat, value = cpm_seurat@meta.data$status)
mean_by_cluster_cpm <- do.call(cbind, lapply(unique(Idents(cpm_seurat)), function(x) {
  apply(GetAssayData(cpm_seurat, slot = "data")[usegenes,colnames(cpm_seurat)][,Idents(cpm_seurat) == x], 1, mean)
}))

colnames(mean_by_cluster_1) <- unique(Idents(data))
design <-as.factor(meta1$stage)
design<- as.factor(targets_2$`MF-Grade`)
rnaseq_data <- rnaseq_data[,1:dim(rnaseq_data)[2]-1]
rownames(mean_by_cluster)<- mean_by_cluster$...1
mean_by_cluster<- mean_by_cluster[,1:dim(mean_by_cluster)[2]-1]


CIBER_barywane<- runCIBERSORT(cpm,mean_by_cluster_cpm,design, mc.cores=1)
annos <- CIBER$CellType
experiment<-CIBER$SampleClass

##apply Cibersortx on average samples (0,1,2,3)

data_1 <- CIBER
rdata <- data_1 %>%
  select(CellType,SampleClass,Fraction) %>%
  reshape2::dcast(CellType~SampleClass,fun.aggregate = sum,value.var = "Fraction")
normrdata <- rdata[,2:5]
normrdata <- t(apply(normrdata, 1, function(x)(x-min(x))/(max(x)-min(x))))
rownames(normrdata) <- rdata$CellType
normrdata<-normrdata[complete.cases(normrdata),]

png("mega_heatmap_cor.png",width = 7,height = 7,units = 'in',res = 300)

ComplexHeatmap::Heatmap(t(normrdata),column_title = "EVs-Decon.", 
                            column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                            col=brewer.pal(9, "Reds"))
dev.off()


png("mega_heatmap_cor.png",width = 7,height = 7,units = 'in',res = 300)
draw(p3c %v% p3)
dev.off()


##apply Cibersortx on all the 16 samples

data_1 <- CIBER
rdata_1 <- data_1 %>%
  select(CellType,SampleID,Fraction) %>%
  reshape2::dcast(CellType~SampleID,fun.aggregate = sum,value.var = "Fraction")
normrdata_1 <- rdata_1[,2:17]
normrdata_1 <- t(apply(normrdata_1, 1, function(x)(x-min(x))/(max(x)-min(x))))
rownames(normrdata_1) <- rdata_1$CellType
normrdata_1<-normrdata_1[complete.cases(normrdata_1),]
sampleLabels<- targets_2$`MF-Grade`
colnames(normrdata_1)<- sampleLabels
ComplexHeatmap::Heatmap(t(normrdata_1),column_title = "EVs-Decon.", 
                        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                        col=brewer.pal(9, "Reds"))


##boxplot fraction

ggplot(CIBER,aes(x=SampleClass,y=Fraction,fill=CellType))+
  geom_boxplot()+
  ggtitle('All clusters per region')+
  facet_wrap(.~CellType, scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


##barplot


paletteLength = 32
myColor = colorRampPalette(c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC"))(paletteLength)


ggplot(CIBER,aes(x=SampleID,y=Fraction,fill=CellType))+
  geom_bar(position="fill", stat="identity")+
  facet_wrap(.~SampleClass,scales = 'free_x')+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))+scale_fill_manual(values = myColor)+
scale_color_manual(values = myColor)



