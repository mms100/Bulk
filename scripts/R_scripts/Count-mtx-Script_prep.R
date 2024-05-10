
#load the need library ----

library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Mmusculus.v79) #replace with your organism-specific database package
library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot)
library(tximeta)
library(SummarizedExperiment)
#import files ----
targets<- read.csv(file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/target.csv")
path<- file.path(paste0("/Volumes/groups/Group_RS/Csnk1a1 omics/Csnk1a1_Haploinsufficiency_+_P53_KO"), targets$File.name,"quant.sf")
setwd("/Volumes/groups/Group_RS/Csnk1a1 omics/Csnk1a1_Haploinsufficiency_+_P53_KO/00_1_Csnk_p53_primary_leukemic_specimen/")
all(file.exists(path)) 

#build matrix with tpm, cpm and raw count and store it as csv

Txi_gene<- tximeta(path, type = "salmon",txOut=TRUE )

gse <- summarizeToGene(Txi_gene)
counts<- assays(gse)[["counts"]]
rownames(counts)<- gse@rowRanges@elementMetadata@listData[["gene_name"]]
counts<- as.data.frame(counts)
sample_label <- targets$Sample.name
colnames(counts)<- sample_label
tpm<- assays(gse)[["abundance"]]
rownames(tpm)<- gse@rowRanges@elementMetadata@listData[["gene_name"]]
tpm<- as.data.frame(tpm)
colnames(tpm)<- sample_label


write.table(counts, file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/Count_normalized_matrix/counts.txt", quote = T, row.names = T, sep = "\t", col.names = NA)
write.table(tpm, file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/Count_normalized_matrix/tpms.txt", quote = T, row.names = T, sep = "\t", col.names = NA)

#then make cpm ----

# Make a DGElist from your counts, and plot ----
myDGEList <- DGEList(counts)
# take a look at the DGEList object 
myDGEList
#DEGList objects are a good R data file to consider saving to you working directory
save(myDGEList, file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/Count_normalized_matrix/myDGEList")


# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(myDGEList) 
colSums(cpm)
log2.cpm <- cpm(myDGEList, log=TRUE)


write.table(cpm, file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/Count_normalized_matrix/cpms.txt", quote = T, row.names = T, sep = "\t", col.names = NA)
write.table(log2.cpm, file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/Count_normalized_matrix/log2_cpm.txt", quote = T, row.names = T, sep = "\t", col.names = NA)



#plotting some genes specific for Mk_Ery_myeloid_pro


#load the count 

targets <- targets[, c(2,3)]

rownames(targets) <- targets$Sample.name

counts <- read_delim(file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/Count_normalized_matrix/counts.txt")

counts <- column_to_rownames(counts, var = "...1")

fake_cpm_seurat <- CreateSeuratObject(counts, meta.data = targets)
fake_cpm_seurat$Condition<- c("Blast", "Blast", "Blast", "Blast", "WT", "WT", "WT", "WT",
                              "WT", "WT")

#apply the dittoseq

Vector <- 
  list(Ery_MK_Prog1=c("Car2", "Eif5a"),
       Ery_MK_Prog2=c("Npm1","Ddx5"),#Velten CLOUD seq
       ImmatureMyeloProg_5=c("TK1","MYBL2"),
       ImmatureMyeloProg_6=c("Eno1","TCF19")#Velten CLOUD seq
  )

#convert from human to mouse markers if you have as list

Vector_mouse <- list()
for(i in 1:length(Vector)){
  Vector_mouse[[i]] <- convert_human_to_mouse_symbols(Vector[[i]])
  names(Vector_mouse)[[i]] <- names(Vector)[[i]]
}


#store all the markeres in one vector

#convert the list to dataframe

#store all the markeres in one vector

HSPCs=c("CD34", "AVP", "CRHBP"),
MEP_CMP=c("TAL1","GATA1","SPI1"),
MEP_1=c("TK1","RRM2","ITGA2B","MYBL2","KLF1"),
MEP_2=c("CSF1","TFR2","CNRIP1"),
MPP=c("KIT" ,"MYC"),
GMP_MDP=c("CSF3R","CSF1R")
EryProg_1=c("CD71","HBB","LMNA","NFIA")

vector <- c("CD34", "AVP", "CRHBP", "TAL1","GATA1","SPI1", "TK1","RRM2","ITGA2B","MYBL2","KLF1", "CSF1","TFR2","CNRIP1",
            "KIT" ,"MYC", "CSF3R","CSF1R", "CD71","HBB","LMNA","NFIA")

Vector_mouse <- convert_human_to_mouse_symbols(vector)


Vector_mouse <- c("Cd34", "Crhbp", "Tal1", "Gata1", "Spi1", "Tk1", "Rrm2", 
  "Itga2b", "Mybl2", "Klf1", "Csf1", "Tfr2", "Cnrip1", "Kit", "Myc", 
  "Csf3r", "Csf1r", "Tfrc", "Lmna", "Nfia")


dittoDotPlot(fake_cpm_seurat, common_genes,group.by  = "Sample.name" , min.color = "#1984c5", max.color = "#c23728", size=3, y.reorder = c(1:10))+
  coord_flip()+theme(axis.text.y = element_text(face= "bold", size= 10))+theme(axis.text.x = element_text(face= "bold", size= 7))

fake_cpm_seurat <- NormalizeData(fake_cpm_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(fake_cpm_seurat)

fake_cpm_seurat <- ScaleData(fake_cpm_seurat, features = all.genes)


dittoBoxPlot(fake_cpm_seurat, c("Tk1", "Hbb.bs", "Mybl2", "Batf3", "Tcf19"), group.by = "Condition", boxplot.lineweight = 0.1, main = "myeloid_progenitors_2")+
  theme(axis.text.x = element_text(face= "bold", size= 7))


dittoBoxPlot(fake_cpm_seurat, c("Cd34", "Crhbp"), group.by = "Condition", boxplot.lineweight = 0.1,, main = "HSPCs")+
  theme(axis.text.x = element_text(face= "bold", size= 7))

dittoBoxPlot(fake_cpm_seurat, c("Tal1", "Gata1", "Spi1"), group.by = "Condition", boxplot.lineweight = 0.1,main= "MEP_CMP")
+theme(axis.text.x = element_text(face= "bold", size= 7))

dittoBoxPlot(fake_cpm_seurat, c("Tk1", "Rrm2","Itga2b", "Mybl2", "Klf1"), group.by = "Condition", boxplot.lineweight = 0.1,main= "MEP_1")+
  theme(axis.text.x = element_text(face= "bold", size= 7))
dittoBoxPlot(fake_cpm_seurat, c("Csf1", "Tfr2", "Cnrip1"), group.by = "Condition", boxplot.lineweight = 0.1, main = "MEP_2")+
  theme(axis.text.x = element_text(face= "bold", size= 7))

dittoBoxPlot(fake_cpm_seurat, c("Kit", "Myc"), group.by = "Condition", boxplot.lineweight = 0.1, main= "MPP")+
  theme(axis.text.x = element_text(face= "bold", size= 7))

dittoBoxPlot(fake_cpm_seurat, c("Csf3r", "Csf1r"), group.by = "Condition", boxplot.lineweight = 0.1, main = "GMP_MDP")+
  theme(axis.text.x = element_text(face= "bold", size= 7))
                                  > dittoBoxPlot(fake_cpm_seurat, c("Csf3r", "Csf1r"), group.by = "Condition", boxplot.lineweight = 0.1,)+theme(axis.text.x = element_text(face= "bold", size= 7))

dittoBoxPlot(fake_cpm_seurat, c("Tfrc", "Lmna", "Nfia"), group.by = "Condition", boxplot.lineweight = 0.1, main = "EryProg_1")+
  theme(axis.text.x = element_text(face= "bold", size= 7))

dittoBoxPlot(fake_cpm_seurat, c("Cd34", "Crhbp"), group.by = "Condition", boxplot.lineweight = 0.1,)+theme(axis.text.x = element_text(face= "bold", size= 7))

dittoBoxPlot(fake_cpm_seurat, c("Tk1", "Spink2", "Mybl2", "Batf3", "Tcf19"), group.by = "Condition", boxplot.lineweight = 0.1, main = "ImmatureMyeloProg_2")+
  theme(axis.text.x = element_text(face= "bold", size= 7))

print(P1, P2)
fake_cpm_seurat$Condition<- c("Blast", "Blast", "Blast", "Blast", "WT", "WT", "WT", "WT",
"WT", "WT")
####ridge_plot####

#for the ridge plot it is nicer to use the normalized count matrix
DGElist <- load(file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/Count_normalized_matrix/myDGEList")
cpm <- cpm(myDGEList)
fake_cpm_seurat <- CreateSeuratObject(counts, meta.data = targets)

#import the markergenelist
Hass_marker<- read.csv("/Users/msaad/Desktop/Copy of Haas markers.csv")
Hass_marker <- Hass_marker[-1,]

genes_list <- list()

for(i in 1: 32){
  genes_list[[i]] <- Hass_marker[,i]
  genes_list[[i]]<- genes_list[[i]][nzchar(genes_list[[i]])]
  names(genes_list)[[i]] <- names(Hass_marker)[[i]]
}

scrna <- fake_cpm_seurat
#plot the ridge plot
for(i in 1:32){
  scrna <- AddModuleScore(object = scrna,features = as.list(genes_list[i]),name = names(genes_list)[i])
}


modulenames <- paste0(names(genes_list),"1")

plist <- list()
for(i in 1:length(modulenames)){
  
  if(modulenames[i] %in% colnames(scrna@meta.data)){
    p <- RidgePlot(scrna,features = modulenames[i],group.by = "Condition")+
      geom_point(alpha = 0.5)+
      ggtitle(names(genes_list)[i])+
      coord_cartesian(clip = "off") +
      geom_density_ridges(quantile_lines=TRUE,
                          quantile_fun=function(x,...)mean(x))
  } else {}
  
  ifelse(modulenames[i] %in% colnames(scrna@meta.data),
         yes = plist[[i]] <- p,
         no = print(paste0("Module name ",modulenames[i]," not scored")))
}


names(plist) <- names(genes_list)


pdf(file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Query/Blast_analysis/marker_plot/ridge_plot.pdf")


for(i in 1:length(modulenames)){
  print(plist[[i]])
}

dev.off()


