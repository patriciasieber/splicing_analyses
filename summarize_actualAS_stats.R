## sum stats for actual AS events

wdir <- "/media/patricia/as_paper/"
setwd(wdir)

species_folders <- list.dirs(recursive=F)


sum_stats <- function(csv_path){
  ## read stress/ infection AS results
  as_input <- read.csv(csv_path,sep="\t",header=T)
  
  ffile <- file(paste0(csv_path,"_stats"),"w")
  writeLines(paste0("Number of AS genes in total: ",nrow(as_input)),ffile)
  
  if(nrow(as_input) > 0){
    ## make unique, count genes with multiple events
    as_input_unique <- unique(as_input)
    unique_genes <- as.character(unique(as_input$gene))
    
    writeLines(paste0("Number of AS genes unique: ",length(unique_genes)),ffile)
    
    if(nrow(as_input_unique) != length(unique_genes)){
      writeLines(paste0((nrow(as_input_unique)-length(unique_genes))," genes appear multiple times with different AS patterns"),ffile)
    }
    
    ## print multiple genes
    multiple <- as.data.frame(table(as_input$gene))
    multi <- as.character(multiple[multiple$Freq > 1,]$Var1)
    if(length(multi) > 0){
      multi <- paste(multi,collapse=",")
      writeLines(paste0("Genes with multiple occurrence: ",multi),ffile)
    }
    
    ## count how many as are supported by stringtie (per condition)
    stringtie_supp <- as.character(as_input_unique$stringtie_supported)
    stringtie_supp <- unlist(lapply(stringtie_supp,function(x){return(strsplit(x,",")[[1]])}))
    count_stringtie_supp <- as.data.frame(table(stringtie_supp))
    count_y <- count_stringtie_supp[count_stringtie_supp$stringtie_supp == "y",]$Freq
    count_n <- count_stringtie_supp[count_stringtie_supp$stringtie_supp == "n",]$Freq
    
    if(is.null(count_y)) count_y <- 0
    if(is.null(count_n)) count_n <- 0
    
    writeLines(paste0("AS is supported by Stringtie prediction in ",count_y," cases"),ffile)
    writeLines(paste0("AS is not supported by Stringtie prediction in ",count_n," cases"),ffile)
    
    ## count AS pattern per condition (control, infection, stress; and only infection and stress)
    as_pattern_control <- as.character(as_input_unique$AS_pattern_cond1)
    as_pattern_treatment <- as.character(as_input_unique$AS_pattern_cond2)
    
    as_pattern_control <- unlist(lapply(as_pattern_control,function(x){return(strsplit(x,",")[[1]])}))
    as_pattern_treatment <- unlist(lapply(as_pattern_treatment,function(x){return(strsplit(x,",")[[1]])}))
    
    table_as_pattern_control <- as.data.frame(table(as_pattern_control))
    table_as_pattern_treatment <- as.data.frame(table(as_pattern_treatment))
    writeLines("-------------------------------",ffile)
    writeLines(paste0("AS pattern control: in total ",sum(table_as_pattern_control$Freq)),ffile)
    if(nrow(table_as_pattern_control) > 0){
      writeLines(paste(table_as_pattern_control$as_pattern_control,table_as_pattern_control$Freq,collapse="\t"),ffile)
    }
    writeLines("-------------------------------",ffile)
    writeLines(paste0("AS pattern treatment: in total ",sum(table_as_pattern_treatment$Freq)),ffile)
    if(nrow(table_as_pattern_treatment) > 0){
      writeLines(paste(table_as_pattern_treatment$as_pattern_treatment,table_as_pattern_treatment$Freq,collapse="\t"),ffile)
    }
    writeLines("-------------------------------",ffile)
    as_pattern_all <- c(as_pattern_control,as_pattern_treatment)
    table_as_pattern_all <- as.data.frame(table(as_pattern_all))
    writeLines(paste0("AS pattern all: in total ",sum(table_as_pattern_all$Freq)),ffile)
    if(nrow(table_as_pattern_all) > 0){
      writeLines(paste(table_as_pattern_all$as_pattern_all,table_as_pattern_all$Freq,collapse="\t"),ffile)
    }
    writeLines("-------------------------------",ffile)
  }
  close(ffile)
}




for(i in species_folders){
  setwd(paste0(wdir,i))
  sum_stats("AS_filter_results_infection.csv")
  sum_stats("AS_filter_results_stress.csv")
}




###
## print AS pattern merged:
astypes_names <- c("AFE","ALE","A5E","A3E","ES","IR","MXE")
astypes_afum <- c(12,6,12,19,4,21,0)
astypes_calb <- c(10,0,1,6,0,24,0)
astypes_cgla <- c(0,0,0,1,0,1,0)
astypes_cpar <- c(2,1,0,0,0,11,0)
astypes_cneo <- c(51,11,0,4,1,14,0)
astypes_lcor <- c(4,0,0,3,0,35,0)
astypes_hcap <- c(13,2,34,49,8,193,0)

astypes_all <- list(astypes_calb,astypes_cpar,astypes_cgla,astypes_afum,astypes_hcap,astypes_cneo,astypes_lcor)

astypes_all_matrix <- matrix(unlist(astypes_all),nrow=7,ncol=7,byrow=F)
colnames(astypes_all_matrix) <- c("Ca","Cp","Cg","Af","Hc","Cn","Lc")

barplot(astypes_all_matrix,main="AS patterns for each species",las=2,
        xlab="Species",col=rainbow(7),beside=FALSE,legend.text=astypes_names)


as_sum <- unlist(lapply(astypes_all,function(x){
  return(sum(x))
}))
names(as_sum) <- c("Ca","Cp","Cg","Af","Hc","Cn","Lc")

barplot(as_sum,main="Number of AS events for each species",
        xlab="Species",col="blue",yaxp=c(0, max(as_sum),9))


## proportion of AS patterns
astypes_relative_all <- lapply(astypes_all,function(x){
  x <- unlist(x)
  summed <- sum(x)
  x_rel <- x/summed
  return(x_rel)
})

astypes_relative_all_matrix <- matrix(unlist(astypes_relative_all),nrow=7,ncol=7,byrow=F)
colnames(astypes_relative_all_matrix) <- c("Ca","Cp","Cg","Af","Hc","Cn","Lc")

pdf("/data/species_comparison/acutalASpattern_relative.pdf")
barplot(astypes_relative_all_matrix,main="Relative occurence of AS pattern for each species",
        xlab="Species",col=rainbow(7),beside=F,legend.text=astypes_names,
        xlim=c(0, ncol(astypes_relative_all_matrix) + 3),
        args.legend=list(x=ncol(astypes_relative_all_matrix)+3,y=max(colSums(astypes_relative_all_matrix)),bty = "n"))
dev.off()


