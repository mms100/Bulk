require(ggplot2)

barplot_deconv <- function(CIBER, name,width=10,height=10){
  res<- ggplot(CIBER,aes(x=SampleID,y=Fraction,fill=CellType))+
        geom_bar(position="fill", stat="identity")+
        facet_wrap(.~SampleClass,scales = 'free_x')+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45))
   res
   ggsave(name,width = width,height = height)
   return(res)
}

boxplot_deconv <- function(CIBER, name,width=10,height=10){
  res<- ggplot(CIBER,aes(x=SampleClass,y=Fraction,fill=CellType))+
      geom_boxplot()+
      ggtitle('All clusters per region')+
      facet_wrap(.~CellType, scales = "free")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
   res
   ggsave(name,width = width,height = height)
   return(res)
}

save_ciber <- function(CIBER,name){
  write.csv(CIBER,name)
}

rna_pca <- function(pca,cols,name,width=10,height=10){
  res<-qplot(x = pca$x[,1], y = pca$x[,2], color = cols)+
    scale_color_discrete(name="Sample type") +
    xlab("PC1") +
    ylab("PC2")
  res
  ggsave(name,width,height)
  return(res)
}
