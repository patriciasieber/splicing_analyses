library(gdata)

location_as_analysis <- "/home/bioinformatik/Hacken2017/"

location_as_analysis_dexseq <- "/home/bioinformatik/Hacken2017/output_dexseq/"
location_as_analysis_mats <- "/home/bioinformatik/Hacken2017/output_mats/FDR0.05ILD0.1/"

species_mats <- c("Mus musculus","Homo sapiens","Danio rerio","Nothobranchius furzeri")
species_dexseq <- c("M. musculus","H. sapiens","D. rerio","N. furzeri")
species_shortcut <- c("mmu","hsa","dre","nfu")

all_tissues <- c("blood","brain","liver","skin")

setwd(location_as_analysis)
system("mkdir compare_output_AStools")

gene_list_dexseq <- read.csv(paste0(location_as_analysis_dexseq,"gene_list.csv"))
colnames(gene_list_dexseq) <- c("GeneID","featureID","padj","seqname","start","end","strand","Species","Tissue","AgeComparison","WikiGene.name")
gene_list_mats <- read.csv(paste0(location_as_analysis_mats,"gene_list_fdr+ild.csv"))

geneLists <- c("gene_list_dexseq","gene_list_mats")
for(k in geneLists){
  gene_list <- get(k)
  if(k == "gene_list_dexseq") species <- species_dexseq
  if(k == "gene_list_mats") species <- species_mats
  for(j in 1:length(species)){
    x <- species[j]
    x_short <- species_shortcut[j]
    current_species <- gene_list[gene_list$Species == x,]
    
    ## for each species and tissue ####
    current_tissues <- as.data.frame(table(current_species$Tissue))
    current_tissues <- current_tissues[current_tissues$Freq != 0,]
    names_tissue <- as.character(current_tissues$Var1)
    tissue_gene_list <- lapply(names_tissue,function(y){
      tmp <- current_species[current_species$Tissue == y,]
      return(unique(tmp$GeneID))
    })
    for(i in 1:length(names_tissue)){
      ffile <- file(paste0("compare_output_AStools/unique_genes_",x_short,"_",names_tissue[i],"_",k,".txt"),"w")
      writeLines(as.character(tissue_gene_list[[i]]),ffile)
      close(ffile)
    }
    
    ## for each species and time point ####
    current_ages <- as.data.frame(table(current_species$AgeComparison))
    current_ages <- current_ages[current_ages$Freq != 0,]
    names_ages <- as.character(current_ages$Var1)
    age_gene_list <- lapply(names_ages,function(y){
      tmp <- current_species[current_species$AgeComparison == y,]
      return(unique(tmp$GeneID))
    })
    for(i in 1:length(names_ages)){
      ffile <- file(paste0("compare_output_AStools/unique_genes_",x_short,"_",names_ages[i],"_",k,".txt"),"w")
      writeLines(as.character(age_gene_list[[i]]),ffile)
      close(ffile)
    }
    
    ## for each species, tissue and time point ####
    lapply(names_tissue,function(y){
      current_species_tissue <- current_species[current_species$Tissue == y,]
      lapply(names_ages,function(z){
        current_species_tissue_age <- current_species_tissue[current_species_tissue$AgeComparison == z,]
        ffile <- file(paste0("compare_output_AStools/unique_genes_",x_short,"_",y,"_",z,"_",k,".txt"),"w")
        writeLines(as.character(unique(current_species_tissue_age$GeneID)),ffile)
        close(ffile)
      })
    })
  }
}


## compare outputs for every species, tissue and timepoint for overlap

ffile <- file("compare_output_AStools/overlap_mats_dexseq.csv","w")
writeLines(paste("species","tissue","age_comparison","rMATS","DEXSeq","overlap",sep="\t"),ffile)
for(x in species_mats){
  mats_current_species <- gene_list_mats[gene_list_mats$Species == x,]
  x_dexseq <- species_dexseq[which(species_mats == x)]
  dexseq_current_species <- gene_list_dexseq[gene_list_dexseq$Species == x_dexseq,]
  
  current_tissues <- as.data.frame(table(mats_current_species$Tissue))
  current_tissues <- current_tissues[current_tissues$Freq != 0,]
  current_tissues <- as.character(current_tissues$Var1)
  
  for(y in current_tissues){
    mats_current_species_tissue <- mats_current_species[mats_current_species$Tissue == y,]
    dexseq_current_species_tissue <- dexseq_current_species[dexseq_current_species$Tissue == y,]
    
    current_ages <- as.data.frame(table(mats_current_species_tissue$AgeComparison))
    current_ages <- current_ages[current_ages$Freq != 0,]
    current_ages <- as.character(current_ages$Var1)
    
    current_ages_dexseq <- as.data.frame(table(dexseq_current_species_tissue$AgeComparison))
    current_ages_dexseq <- current_ages_dexseq[current_ages_dexseq$Freq != 0,]
    current_ages_dexseq <- as.character(current_ages_dexseq$Var1)
    
    if(length(current_ages) != length(current_ages_dexseq)){
      if(x == "Mus musculus"){ current_ages_dexseq <- current_ages_dexseq[5:10]
      }else if(x == "Homo sapiens"){ current_ages_dexseq <- current_ages_dexseq[4:6]
      }else print(paste0("different length of comparisons for ",x))
    }
    
    for(z in current_ages){
      z_dexseq <- current_ages_dexseq[which(current_ages == z)]
      
      mats_current_species_tissue_age <- mats_current_species_tissue[mats_current_species_tissue$AgeComparison == z,]
      mats_ids <- unique(as.character(mats_current_species_tissue_age$GeneID))
      dexseq_current_species_tissue_age <- dexseq_current_species_tissue[dexseq_current_species_tissue$AgeComparison == z_dexseq,]
      dexseq_ids <- unique(as.character(dexseq_current_species_tissue_age$GeneID))
      
      overlap <- unlist(lapply(mats_ids,function(w){
        return(dexseq_ids[dexseq_ids == w])
      }))
      overlap_ids <- unique(overlap)
      
      writeLines(paste(x,y,z,length(mats_ids),length(dexseq_ids),length(overlap_ids),sep="\t"),ffile)
    }
  }
}
close(ffile)



## for gene depending on homology ####

## read homology table
# homology <- read.csv("/home/bioinformatik/Hacken2017/homology/homologyTable_mmu_hsa_dre_ensemblIDs_nfu_geneIDs.csv",
#                      stringsAsFactors=F)
# 
# ## build subsets based on homology: hs_mm_dr_nf, hs_mm_dr, hs_mm, only mm
# hom_hs_mm_dr_nf <- homology[(nchar(homology$hsa)>0 & nchar(homology$mmu)>0 & nchar(homology$dre)>0 & nchar(homology$nfu)>0),]
# hom_hs_mm_dr <- homology[(nchar(homology$hsa)>0 & nchar(homology$mmu)>0 & nchar(homology$dre)>0 & nchar(homology$nfu)==0),]
# hom_hs_mm <- homology[(nchar(homology$hsa)>0 & nchar(homology$mmu)>0 & nchar(homology$dre)==0 & nchar(homology$nfu)==0),]
# hom_mm <- homology[(nchar(homology$hsa)==0 & nchar(homology$mmu)>0 & nchar(homology$dre)==0 & nchar(homology$nfu)==0),]
# 
# 
# ## homology for each tissue ####
# for(j in 1:length(all_tissues)){
#   print(j)
#   current_tissue <- all_tissues[j]
#   
#   current_gene_list <- gene_list[gene_list$Tissue == current_tissue,]
#   
#   current_species <- as.data.frame(table(current_gene_list$Species))
#   current_species <- current_species[current_species$Freq != 0,]
#   
#   names_species <- as.character(current_species$Var1)
#   
#   species_gene_list <- lapply(names_species,function(y){
#     tmp <- current_gene_list[current_gene_list$Species == y,]
#     return(unique(tmp$GeneID))
#   })
#   
#   #if(any(names_species == "Mus musculus")){
#   ## divide list into sublists of homology
#   #sub_hs_mm_dr_nf <- hom_hs_mm_dr_nf[hom_hs_mm_dr_nf$mmu.Ensembl.ID %in% as.character(species_gene_list[[which(names_species == "Mus musculus")]]),]
#   #sub_hs_mm_dr <- hom_hs_mm_dr[hom_hs_mm_dr$mmu.Ensembl.ID %in% as.character(species_gene_list[[which(names_species == "Mus musculus")]]),]
#   #sub_hs_mm <- hom_hs_mm[hom_hs_mm$mmu.Ensembl.ID %in% as.character(species_gene_list[[which(names_species == "Mus musculus")]]),]
#   #sub_mm <- hom_mm[hom_mm$mmu.Ensembl.ID %in% as.character(species_gene_list[[which(names_species == "Mus musculus")]]),]
#   
#   
#   ## nach species aufspalten und für jedes subset
#   for(i in 1:length(names_species)){
#     print(i)
#     check_species <- names_species[i]
#     check_species_shortcut <- species_shortcut[which(species == check_species)]
#     current_gene_list <- as.character(species_gene_list[[i]])
#     
#     ## sub_hs_mm_dr_nf
#     current_sub_hs_mm_dr_nf <- hom_hs_mm_dr_nf[hom_hs_mm_dr_nf[[which(names(hom_hs_mm_dr_nf) == check_species_shortcut)]] %in% current_gene_list,]$mmu.Ensembl.ID  
#     if(length(current_sub_hs_mm_dr_nf)>0){
#       ffile <- file(paste0("homology_sub_hs_mm_dr_nf_unique_genes_",current_tissue,"_",check_species_shortcut,".txt"),"w")
#       writeLines(as.character(current_sub_hs_mm_dr_nf),ffile)
#       close(ffile)
#     }
#     
#     ## sub_hs_mm_dr
#     current_sub_hs_mm_dr <- hom_hs_mm_dr[hom_hs_mm_dr[[which(names(hom_hs_mm_dr) == check_species_shortcut)]] %in% current_gene_list,]$mmu.Ensembl.ID  
#     if(length(current_sub_hs_mm_dr)>0){
#       ffile <- file(paste0("homology_sub_hs_mm_dr_unique_genes_",current_tissue,"_",check_species_shortcut,".txt"),"w")
#       writeLines(as.character(current_sub_hs_mm_dr),ffile)
#       close(ffile)
#     }
#     
#     ## sub_hs_mm
#     current_sub_hs_mm <- hom_hs_mm[hom_hs_mm[[which(names(hom_hs_mm) == check_species_shortcut)]] %in% current_gene_list,]$mmu.Ensembl.ID  
#     if(length(current_sub_hs_mm)>0){
#       ffile <- file(paste0("homology_sub_hs_mm_unique_genes_",current_tissue,"_",check_species_shortcut,".txt"),"w")
#       writeLines(as.character(current_sub_hs_mm),ffile)
#       close(ffile)
#     }
#     
#     ## sub_mm
#     current_sub_mm <- hom_mm[hom_mm[[which(names(hom_mm) == check_species_shortcut)]] %in% current_gene_list,]$mmu.Ensembl.ID  
#     if(length(current_sub_mm)>0){
#       ffile <- file(paste0("homology_sub_mm_unique_genes_",current_tissue,"_",check_species_shortcut,".txt"),"w")
#       writeLines(as.character(current_sub_mm),ffile)
#       close(ffile)
#     }
#   }
#   #}
# }
# 
# 

