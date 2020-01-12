library("gdata")
library("VennDiagram")

## venn diagram of overlapping as genes in different tissues

setwd("/home/aging/output_mats/results_mats_fromHolyList/")

species <- c("dre","hsa","mmu","nfu")

v1_all <- list(dre_v1=c("12vs36","24vs36"),hsa_v1=c("29vs65","50vs65"),mmu_v1=c("09vs24","15vs24"),nfu_v1=c("12vs27","20vs27"))
v2_all <- list(dre_v2=c("12vs42","24vs42"),hsa_v2=c("29vs80","50vs80"),mmu_v2=c("09vs30","15vs30"),nfu_v2=c("12vs39","20vs39"))
v3_all <- list(dre_v3=c("36vs42"),hsa_v3=c("65vs80"),mmu_v3=c("24vs30"),nfu_v3=c("27vs39"))

holy_list <- read.csv("/home/aging/output_mats/mats_holy_list.csv")

for(sp in species){
  holy_list_sp <- holy_list[holy_list$species == sp,]
  tissues <- as.character(unique(holy_list_sp$tissue))
  
  as_genes <- lapply(tissues,function(x){as.character(holy_list_sp[holy_list_sp$tissue == x,]$gene_id)})
  
  pdf(paste0("venndiagram_tissue_corrected_overlap_",sp,".pdf"))
  venn.plot <- venn.diagram(as_genes,NULL,
                            sub.fontfamily="serif",
                            #alpha=c(0.5,0.5),
                            cex=2,cat.fontface=4,
                            category.names=tissues)
  grid.draw(venn.plot)
  dev.off()
}
  
  
## venn diagrams not merged for time points but each comparison v1,v2,v3 separately to compare the tissues of each species
for(sp in species){
  holy_list_sp <- holy_list[holy_list$species == sp,]
  tissues <- as.character(unique(holy_list_sp$tissue))
  
  comparisons_v <- paste0(rep(sp,3),"_",c("v1","v2","v3"))
  current_v1 <- as.character(unlist(v1_all[names(v1_all) == comparisons_v[1]]))
  current_v2 <- as.character(unlist(v2_all[names(v2_all) == comparisons_v[2]]))
  current_v3 <- as.character(unlist(v3_all[names(v3_all) == comparisons_v[3]]))
  
  
  holy_list_sp_v1 <- holy_list_sp[holy_list_sp$comparison == current_v1[1] | holy_list_sp$comparison == current_v1[2],] 
  holy_list_sp_v2 <- holy_list_sp[holy_list_sp$comparison == current_v2[1] | holy_list_sp$comparison == current_v2[2],]
  holy_list_sp_v3 <- holy_list_sp[holy_list_sp$comparison == current_v3,]  
  
  as_genes_list_v1 <- NULL
  as_genes_list_v2 <- NULL
  as_genes_list_v3 <- NULL
  
  for(i in 1:length(tissues)){
    ti <- tissues[i]
    
    holy_list_sp_ti_v1 <- holy_list_sp_v1[holy_list_sp_v1$tissue == ti,]
    current_genes <- as.character(unique(holy_list_sp_ti_v1$gene_id))
    as_genes_list_v1[[i]] <- current_genes
    writeLines(current_genes,paste0("mats_results_corrected_",sp,"_",ti,"_v1"))
    
    holy_list_sp_ti_v2 <- holy_list_sp_v2[holy_list_sp_v2$tissue == ti,]
    current_genes <- as.character(unique(holy_list_sp_ti_v2$gene_id))
    as_genes_list_v2[[i]] <- current_genes
    writeLines(current_genes,paste0("mats_results_corrected_",sp,"_",ti,"_v2"))
    
    holy_list_sp_ti_v3 <- holy_list_sp_v3[holy_list_sp_v3$tissue == ti,]
    current_genes <- as.character(unique(holy_list_sp_ti_v3$gene_id))
    as_genes_list_v3[[i]] <- current_genes
    writeLines(current_genes,paste0("mats_results_corrected_",sp,"_",ti,"_v3"))
  }
    
    pdf(paste0("venndiagram_tissue_corrected_v1_overlap_",sp,".pdf"))
    venn.plot <- venn.diagram(as_genes_list_v1,NULL,
                              sub.fontfamily="serif",
                              #alpha=c(0.5,0.5),
                              cex=2,cat.fontface=4,
                              category.names=tissues)
    grid.draw(venn.plot)
    dev.off()  
  
    
    pdf(paste0("venndiagram_tissue_corrected_v2_overlap_",sp,".pdf"))
    venn.plot <- venn.diagram(as_genes_list_v2,NULL,
                              sub.fontfamily="serif",
                              #alpha=c(0.5,0.5),
                              cex=2,cat.fontface=4,
                              category.names=tissues)
    grid.draw(venn.plot)
    dev.off() 
    
    
    pdf(paste0("venndiagram_tissue_corrected_v3_overlap_",sp,".pdf"))
    venn.plot <- venn.diagram(as_genes_list_v3,NULL,
                              sub.fontfamily="serif",
                              #alpha=c(0.5,0.5),
                              cex=2,cat.fontface=4,
                              category.names=tissues)
    grid.draw(venn.plot)
    dev.off() 
}
