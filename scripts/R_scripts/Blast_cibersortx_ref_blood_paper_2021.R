#load needed libararies
library(readr)
library(tidyverse)
library(limma)
library(edgeR)
library(DESeq2)
library(ense)
library(org.Mm.eg.db)



#import the dataset

counts <- read.delim(file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Refereces/Blood_paper_2021/GSE147265_RSEM-counts.txt")

counts <- counts[,-1]
#aggreagate the sum of the duplicated genes symbols


result <- counts %>%
  group_by(Symbol) %>%
  summarise(across(names(counts)[1:61], sum))

result <- as.data.frame(result, check.names = F)

result <- column_to_rownames(result, var = "Symbol")

#load metadata

META <- as.data.frame(read_excel("/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Refereces/Blood_paper_2021/META.xlsx"))


rownames(META) <- META$`Sequencing ID`

#remove any entries that is not present in the count matrix


META <- META[(rownames(META) %in% colnames(result) ), ]


#convert count matrix to cpm

myDGEList <- DGEList(result)
# take a look at the DGEList object 
myDGEList


cpm <- cpm(myDGEList) 
colSums(cpm)

#rename the colmnames of the result to the cell subtype

match(META[,"samplenames"], colnames(cpm))

colnames(cpm) <- as.character(META$Diagnosis)

write.table(cpm, file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Refereces/Blood_paper_2021/cpm.txt", quote = T, row.names = T, sep = "\t", col.names = NA)



#plotting results from website

#import dataset
CIBER_MDS <- read.csv(file = "/Volumes/home/Data analysis/MDS-Cibersortx/MDS_R/Refereces/Blood_paper_2021/CIBERSORTx_Job165_Results.csv")

CIBER_MDS$sample_type <- c(rep(c("Blast"), 4), rep("WT", 6))

CIBER_MDS_pivoted<- pivot_longer(CIBER_MDS,
                                   cols = c(2:6),
                                   names_to = "CellType",
                                   values_to = "Fraction"
)


ggplot(CIBER_MDS_pivoted,aes(x=sample_type,y=Fraction,fill=CellType))+
  geom_boxplot()+
  ggtitle('All clusters per Genotype')+
  facet_wrap(.~CellType, scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, size = 10  ,vjust = 0.5, hjust=1), text = element_text(size = 10))

