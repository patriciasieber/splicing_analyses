library(pheatmap)
library(RColorBrewer)

## matrix for each species with values of v1,v2,v3 for heatmap

#setwd("/mnt/dessertlocal/emanuel/hacken_alternative_splicing/differential_splicing/results_mats/")
setwd("/home/trice/Phd/Hacken2017/output_mats/results_mats_fromHolyList//")


species <- c("dre","hsa","mmu","nfu")

holy_as_matrix <- matrix(nrow=4,ncol=12)
rownames(holy_as_matrix) <- c("blood","brain","liver","skin")
repeat_species <- c(rep("dre",3),rep("hsa",3),rep("mmu",3),rep("nfu",3))
repeat_comp <- rep(c("v1","v2","v3"),4)
colnames(holy_as_matrix) <- paste(repeat_species,repeat_comp,sep="_")


for(sp in species){
  ffiles <- list.files(pattern=paste0("mats_results_corrected_",sp))
  ffiles <- ffiles[grep("_v",ffiles)]
  tissues <- unique(unlist(lapply(ffiles,function(x){return(strsplit(x,"_")[[1]][5])})))
  
  ## matrix with tissues in rows, comparisons in columns
  as_matrix <- matrix(nrow=4,ncol=3)
  rownames(as_matrix) <- c("blood","brain","liver","skin")
  colnames(as_matrix) <- c("v1","v2","v3")
  
  for(i in ffiles){
    splitt <- strsplit(i,"_")[[1]]
    row_nr <- which(rownames(as_matrix) == as.character(splitt[5]))
    col_nr <- which(colnames(as_matrix) == as.character(splitt[6]))
    
    current_entry <- readLines(i)
    
    as_matrix[row_nr,col_nr] <- length(current_entry)
  }
  
  
  holy_as_matrix[,grep(sp,colnames(holy_as_matrix))] <- as_matrix
  
  pdf(paste0("heatmap_",sp,".pdf"),onefile=F)
  pheatmap(as_matrix,cellwidth=40,cellheight=40,cluster_rows=F,cluster_cols=F,
           scale="none",display_numbers=T,margins=c(10,20),
           main=paste0("Alternativly spliced genes in ",sp))
  dev.off()
}

## mmu brain has 0 AS in brain
holy_as_matrix[2,c(7:9)] <- c(0,0,0)



pdf(paste0("heatmap_as_all_species.pdf"),width=10,height=10,onefile=F)
pheatmap(holy_as_matrix,cellwidth=40,cellheight=40,cluster_rows=F,cluster_cols=F,
         scale="none",display_numbers=T,margins=c(10,20),
         main="Alternativly spliced genes in all species")
dev.off()




paletteLength <- 20
myColor <- colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(paletteLength)  ## RdYlBu,RdBu,Spectral
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(holy_as_matrix,na.rm=T), 60, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(holy_as_matrix,na.rm=T)/paletteLength, max(holy_as_matrix,na.rm=T), length.out=floor(paletteLength/2)))
pdf("heatmap_as_all_species_newcolors.pdf",width=10,height=10,onefile=F)
pheatmap(holy_as_matrix,cellwidth=40,cellheight=40,
         cluster_rows=F,cluster_cols=F,color=myColor,breaks=myBreaks,
         scale="none",display_numbers=T,margins=c(10,20),
         main="Alternativly spliced genes in all species")
#color=brewer.pal(9, "Blues")
dev.off()
