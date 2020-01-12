location_as_analysis <- "/mnt/alternative_splicing/"
setwd(location_as_analysis)


import.as.events <- function(condition){
  import_mats_output <- function(path_to_mats_output,fdr=0.05,ild=0.1){
    mats_files <- list.files(path_to_mats_output)
    mats_table <- NULL
    for(i in mats_files){
      current_table <- read.table(paste0(path_to_mats_output,i),header=TRUE,sep="\t",fill=T,quote="")
      sub_table <- current_table[c("GeneID","chr","strand","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")]
      mats_table <- rbind(mats_table,sub_table)
    }
    mats_table <- mats_table[mats_table$FDR < fdr,]
    mats_table <- mats_table[abs(as.numeric(mats_table$IncLevelDifference)) > ild,]
    return(mats_table)
  }
  
  mats <- import_mats_output(paste0(condition,"/MATS_output/"),0.05,0.1)
  mats <- unique(mats$GeneID)
  print(paste(condition,", AS events: ",length(mats),sep=""))
  return(as.character(mats))
}


species <- c("dre","hsa","mmu","nfu")
tissues <- c("blood","brain","liver","skin")

for(sp in species){
  print(sp)
  for(ti in tissues){
    print(ti)
    current_dir <- paste0(location_as_analysis,sp,"/",ti,"/")
    if(dir.exists(current_dir)){
      setwd(current_dir)
      comparisons <- list.dirs(path=".",recursive=F,full.names=F)
      
      mats_results <- lapply(comparisons,function(x){
        as_genes <- import.as.events(x)
        writeLines(as_genes,paste0(location_as_analysis,"/results_mats/mats_results_",sp,"_",ti,"_",x))
        return(as_genes)
      })

      ## group results
      ## V1: normal aging (2-4 and 3-4)
      v1 <- unique(unlist(c(mats_results[1],mats_results[3])))
      ## V2: old aging (2-5 and 3-5)
      v2 <- unique(unlist(c(mats_results[2],mats_results[4])))
      ## V3: very old aging (4-5)
      v3 <- unique(unlist(mats_results[5]))
      
      writeLines(v1,paste0(location_as_analysis,"/results_mats/mats_results_",sp,"_",ti,"_v1"))
      writeLines(v2,paste0(location_as_analysis,"/results_mats/mats_results_",sp,"_",ti,"_v2"))
      writeLines(v3,paste0(location_as_analysis,"/results_mats/mats_results_",sp,"_",ti,"_v3"))
      
    }
  }
}



## csv of all species, tissues, and splice patterns ## the holy list ####
location_as_analysis <- "/mnt/alternative_splicing/"
setwd(location_as_analysis)

species <- c("dre","hsa","mmu","nfu")
tissues <- c("blood","brain","liver","skin")

holy_list <- NULL
for(sp in species){
  print(sp)
  for(ti in tissues){
    print(ti)
    current_dir <- paste0(location_as_analysis,sp,"/",ti,"/")
    if(dir.exists(current_dir)){
      setwd(current_dir)
      comparisons <- list.dirs(path=".",recursive=F,full.names=F)
      for(comp in comparisons){
        mats_files <- list.files(path=paste0(current_dir,"/",comp,"/MATS_output/"),pattern="ReadsOnTargetAndJunctionCounts.txt$",full.names=F)
        for(current_file in mats_files){
          current_table <- read.table(paste0(current_dir,"/",comp,"/MATS_output/",current_file),header=TRUE,sep="\t",fill=T,quote="")
          mats_table <- current_table[c("GeneID","chr","strand","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")]
          mats_table <- mats_table[mats_table$FDR < 0.05,]
          mats_table <- mats_table[abs(as.numeric(mats_table$IncLevelDifference)) > 0.1,]
          if(nrow(mats_table) > 0){
            aspattern <- strsplit(current_file,"[.]")[[1]][1]
            len <- nrow(mats_table)
            new_entry <- cbind(rep(sp,len),rep(ti,len),rep(comp,len),rep(aspattern,len),as.character(mats_table$GeneID),as.character(mats_table$chr),as.character(mats_table$strand),mats_table$PValue,mats_table$FDR,
                               mats_table$IncLevel1,mats_table$IncLevel2,mats_table$IncLevelDifference)
            holy_list <- rbind(holy_list,new_entry)
          }
        }
      }
    }
  }
}

colnames(holy_list) <- c("species","tissue","comparison","splicing_pattern","gene_id","chr","strand","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")

write.csv(holy_list,file="/mnt/alternative_splicing/results_mats/mats_holy_list.csv",row.names=F,quote=T)



