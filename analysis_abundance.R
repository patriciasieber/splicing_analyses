
## read abundance of each transcript per sample

## for L_corymbifera
species <- "A_fumigatus"
control <- c(1:3,13:14,17) 
stress <- c(4:12,15:16,18:21)
infection <- c()

conditions <- list(control=control,stress=stress,infection=infection)


#for(species in species_all){
  input_path <- paste0("/data/species_comparison/",species,"/")

  if(file.exists(input_path)){
    setwd(input_path)
  }else stop(paste0("directory ",input_path," does not exist"),call.=FALSE)
  
  path_to_plot_output <- input_path

  ## average for every cond
  cond <- list.files(pattern="gene_abundance")
  nr_cond <- length(cond)
  nr_tr_per_cond <- vector(mode="numeric",length=nr_cond)
  nr_genes_per_cond <- vector(mode="numeric",length=nr_cond)
  
  
  for(i in 1:nr_cond){
    gff <- read.csv(paste0(input_path,"/",cond[i]),sep="\t",skip=2,header=F)
    names(gff) <- c("seqnames","source","type","start","end","score","strand","score2","group")
    ## use only transcript entries and their TPM value
    gff <- gff[gff$type == "transcript",]
    
    ## number of expressed transcripts (tpm > 1)
    tpm <- unlist(lapply(as.character(gff$group),function(x){
      tmp <- (strsplit(x,";")[[1]])
      tmp <- tmp[length(tmp)]
      tmp <- as.numeric(strsplit(tmp," ")[[1]][3])
      return(tmp)
    }))
    gff$tpm <- tpm
    gff_expressed <- gff[gff$tpm > 1,]
    
    ## number of expressed transcripts
    nr_tr_per_cond[i] <- nrow(gff_expressed)
    
    ## number of expressed genes
    gene_id <- unlist(lapply(as.character(gff_expressed$group),function(x){strsplit(strsplit(x,";")[[1]][1]," ")[[1]][2]}))
    gff_expressed$gene_id <- gene_id
    nr_genes_per_cond[i] <- length(unique(gene_id))
  }

  ## write file
  ffile <- file(paste0(path_to_plot_output,"/expressed_transcripts_",species,".txt"),"w")
  writeLines("# expression values with TPM > 1 of all samples of the considered cond",ffile)
  writeLines(paste("cond","expressed_transcripts","expressed_genes",sep="\t"),ffile)
  writeLines(paste(cond,nr_tr_per_cond,nr_genes_per_cond,sep="\t"),ffile)
  close(ffile)
  

  ## average number of expressed genes/transcripts per condition
  print(species)
  for(i in 1:length(conditions)){
    current_name <- names(conditions)[i]
    current_list <- conditions[[i]]
    av_genes <- nr_genes_per_cond[current_list]
    av_genes <- sum(av_genes)/length(av_genes)
    av_tr <- nr_tr_per_cond[current_list]
    av_tr <- sum(av_tr)/length(av_tr)
    
    print(paste("Average expressed genes",current_name,av_genes,sep=" "))
    print(paste("Average expressed transcripts",current_name,av_tr,sep=" "))
  }
  
#}

  
  
## haploid for C.albicans
  species <- "C_albicans"
  control <- c(1:2,5:6,9:10,13:14,17:18,21:22,57:59,95,101,107,113:128,153:160) 
  stress <- c(25:56,63:94,129:144,148:152)
  infection <- c(3:4,7:8,11:12,15:16,19:20,23:24,60:62,96:100,102:106,108:112,145:147)
  
  conditions <- list(control=control,stress=stress,infection=infection)
  
  input_path <- paste0("/data/species_comparison/",species,"/")
  
  if(file.exists(input_path)){
    setwd(input_path)
  }else stop(paste0("directory ",input_path," does not exist"),call.=FALSE)
  
  path_to_plot_output <- input_path
  
  ## average for every cond
  cond <- list.files(pattern="gene_abundance")
  nr_cond <- length(cond)
  nr_tr_per_cond <- vector(mode="numeric",length=nr_cond)
  nr_genes_per_cond <- vector(mode="numeric",length=nr_cond)
  
  
  for(i in 1:nr_cond){
    gff <- read.csv(paste0(input_path,"/",cond[i]),sep="\t",skip=2,header=F)
    names(gff) <- c("seqnames","source","type","start","end","score","strand","score2","group")
    ## use only transcript entries and their TPM value
    gff <- gff[gff$type == "transcript",]
    
    ## number of expressed transcripts (tpm > 1)
    tpm <- unlist(lapply(as.character(gff$group),function(x){
      tmp <- (strsplit(x,";")[[1]])
      tmp <- tmp[length(tmp)]
      tmp <- as.numeric(strsplit(tmp," ")[[1]][3])
      return(tmp)
    }))
    #tpm <- unlist(lapply(as.character(gff$group),function(x){as.numeric(strsplit(strsplit(x,";")[[1]][5]," ")[[1]][3])}))
    gff$tpm <- tpm
    gff_expressed <- gff[gff$tpm > 1,]
    
    tr_id <- unlist(lapply(as.character(gff_expressed$group),function(x){strsplit(strsplit(x,";")[[1]][2]," ")[[1]][3]}))
    
    haplotid_tr_id <- unlist(lapply(tr_id,function(x){
      if(grepl("-T",x)){
        return(substr(x,1,nchar(x)-4))
      }else{return(x)}
    }))
    haplotid_tr_id <- unique(haplotid_tr_id)
    
    ## number of expressed transcripts
    nr_tr_per_cond[i] <- length(haplotid_tr_id)
    
    ## number of expressed genes
    haploid_gene_id <- unlist(lapply(as.character(gff_expressed$group),function(x){
      temp <- strsplit(strsplit(x,";")[[1]][1]," ")[[1]][2]
      temp <- substr(temp,1,nchar(temp)-2)
      }))
    nr_genes_per_cond[i] <- length(unique(haploid_gene_id))
  }
  
  ## write file
  ffile <- file(paste0(path_to_plot_output,"/expressed_transcripts_",species,".txt"),"w")
  writeLines("# expression values with TPM > 1 of all samples of the considered cond",ffile)
  writeLines(paste("cond","expressed_transcripts","expressed_genes",sep="\t"),ffile)
  writeLines(paste(cond,nr_tr_per_cond,nr_genes_per_cond,sep="\t"),ffile)
  close(ffile)
  
  
  ## average number of expressed genes/transcripts per condition
  print(species)
  for(i in 1:length(conditions)){
    current_name <- names(conditions)[i]
    current_list <- conditions[[i]]
    av_genes <- nr_genes_per_cond[current_list]
    av_genes <- sum(av_genes)/length(av_genes)
    av_tr <- nr_tr_per_cond[current_list]
    av_tr <- sum(av_tr)/length(av_tr)
    
    print(paste("Average expressed genes",current_name,av_genes,sep=" "))
    print(paste("Average expressed transcripts",current_name,av_tr,sep=" "))
  }
  
  
