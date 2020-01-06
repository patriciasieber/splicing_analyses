## convert bed file into fasta sequnce, translate it into protein sequcence, check change of protein domains (pfam or prosite)
library(readr)
library("Biostrings")

#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library(rtracklayer)

## read gff/gtf of considered isoforms####
## input: list of gffs of isoforms, first line is transcript name
#input_file <- "path/to/file_list"
#path_to_fasta <- "path/to/fasta"

path_to_fasta <- "/data/PhD/tool_comparison/A_fumigatus/Aspergillus_fumigatusa1163.CADRE.31.dna.toplevel.fa"

working_dir <- "/data/PhD/tool_comparison/A_fumigatus/KatrinLapp/isoforms_gtf_correctedWithCuffdiff/"
gene_name <- "AFUB_057840"
system(paste0("mkdir ",working_dir,"/",gene_name))
setwd(paste0(working_dir,"/",gene_name))


input_file <- paste0(working_dir,"isoforms_paths_",gene_name,".txt")
if(file.exists(input_file)){
  gff_list <- readLines(input_file)
  gff_list <- gff_list[nchar(gff_list) != 0]
  
  gff <- lapply(gff_list,function(x){
    tmp <- import.gff(x) 
    tmp <- tmp[tmp$type == "exon",]
    ## sort entries (the exon of genes on the negative strand might be in the "wrong" order within the annotation file)
    tmp <- sort(tmp)
    return(tmp)
  })
  transcript_names <- unlist(lapply(gff,function(x){return(x$transcript_id[1])}))
  
}else{
  input_file <- paste0(working_dir,"all_isoforms_",gene_name,".gtf")
  if(file.exists(input_file)){
    all_isoforms <- import.gff(input_file)
    all_isoforms <- all_isoforms[all_isoforms$type == "exon",]
    transcript_names <- unique(all_isoforms$transcript_id)
    transcript_names <- transcript_names[!is.na(transcript_names)]
    
    gff <- lapply(transcript_names,function(x){
      tmp <- all_isoforms[all_isoforms$transcript_id == as.character(x),]
      tmp <- sort(tmp)
      return(tmp)
    })
    
  }else{print("Could not read any file for the provided gene name.")}
}

log_file <- file(paste0("protein_changes_",gene_name),"w")
genome <- readDNAStringSet(path_to_fasta,format="fasta")



isoform_sequences <- lapply(gff,function(x){
  exon_sequences <- unlist(lapply(x,function(y){
    chromosome <- as.character(seqnames(y)[1])
    n_chrom_name <- nchar(chromosome)
    chromosome_sequence <- genome[substr(names(genome),1,n_chrom_name) == chromosome,]

    current_seq <- subseq(chromosome_sequence,start=start(y),end=end(y))
    seq_junc <- toString(current_seq)
    return(seq_junc)
  }))
  exon_sequences <- paste(exon_sequences, collapse = '')
  if(as.character(strand(x))[1] == "-"){
    ## return reverse complement
    tmp <- DNAString(exon_sequences)
    tmp <- reverseComplement(tmp)
    exon_sequences <- toString(tmp)
  } 
  return(exon_sequences)
})

writeLines("dna sequences of isoforms:",log_file)
writeLines(paste(transcript_names,isoform_sequences,sep="\n"),log_file)
writeLines("\n",log_file)

## translate into AA sequence ####

isoform_aa <- lapply(isoform_sequences,function(x){
  x1 <- DNAString(x)
  x2 <- DNAString(substr(x,2,nchar(x)))
  x3 <- DNAString(substr(x,3,nchar(x)))
  
  aa1 <- toString(translate(x1))
  aa2 <- toString(translate(x2))
  aa3 <- toString(translate(x3))
  
  ## check for protein-coding sequence
  ## take the longest one
  ## jeweils erste codierende Sequenz wiedergeben, Beginn M, Ende *
  
  cds1 <- substr(aa1,regexpr("M",aa1)[1],regexpr("[*]",aa1)[1]-1)
  cds2 <- substr(aa2,regexpr("M",aa2)[1],regexpr("[*]",aa2)[1]-1)
  cds3 <- substr(aa3,regexpr("M",aa3)[1],regexpr("[*]",aa3)[1]-1)
  
  longest_cds <- max(c(length(cds1),length(cds2),length(cds3)))
  if(longest_cds == 1 ) return(cds1)
  if(longest_cds == 2 ) return(cds2)
  if(longest_cds == 3 ) return(cds3)
})

isoform_aa
writeLines("amino acid sequences of isoforms:",log_file)
writeLines(paste(transcript_names,isoform_aa,sep="\n"),log_file)
writeLines("\n",log_file)

## if equal for all input isoforms
## pairwise comparison
i <- 1
j <- 2
while(i < length(isoform_aa)){
  #print(paste0("i:",i))
  j <- i+1
  while(j <= length(isoform_aa)){
    #print(paste0("j:",j))
    if(isoform_aa[[i]] == isoform_aa[[j]]){
      print(paste0("Isoform ",i," (",transcript_names[i],") and isoform ",j," (",transcript_names[j],") have the same amino acid sequence."))
      writeLines(paste0("Isoform ",i," (",transcript_names[i],") and isoform ",j," (",transcript_names[j],") have the same amino acid sequence."),log_file)
    } 
    j <- j+1
  }
  i <- i+1
}
writeLines("\n",log_file)



if(length(transcript_names)==length(isoform_aa)){
  ffile <- file(paste0(gene_name,"_aa.fa"),"w")
  for(i in 1:length(isoform_aa)){
    if(nchar(isoform_aa[[i]][1]) > 0){
      writeLines(paste0(">",transcript_names[i]),ffile)
      writeLines(as.character(isoform_aa[i]),ffile)
    }
    else{
      print(paste0(transcript_names[i]," does not have an amino acid sequence"))
      writeLines(paste0(transcript_names[i]," does not have an amino acid sequence"),log_file)
    }
  }
  close(ffile)
}else print("number of transcript_names does not correspond to number of isoforms")


## check domains with interproscan ####
## run interproscan 
output_dir <- gene_name
input_fasta <- paste0(gene_name,"_aa.fa")
  
system(paste0("mkdir interproscan_",output_dir))
system(paste0("/home/patricia/bin/interproscan-5.23-62.0/interproscan.sh -i ",input_fasta," -d interproscan_",output_dir," -goterms -iprlookup"))

## check differences
### read tsv and check for different domains, pairwise
#tsv <- read_tsv(paste0("interproscan_",output_dir,"/",gene_name,"_aa.fa.tsv"),col_names=F,fill=T,na = c("", "NA"),quoted_na = TRUE)
tsv <- read.table(paste0("interproscan_",output_dir,"/",gene_name,"_aa.fa.tsv"),header=F,fill=T,sep="\t")

domains_per_isoform <- lapply(transcript_names,function(x){
  tmp <- tsv[tsv$V1 == as.character(x),]
  return(tmp)
})


## compare to the annotated isoform
domains_annotated_isoform <- unique(domains_per_isoform[[1]]$V13)
for(i in 2:length(domains_per_isoform)){
  if(nrow(domains_per_isoform[[i]]) > 0){
    print(paste0("compare protein domains in annotated isoform ",transcript_names[[1]]," with those in predicted isoform ",transcript_names[[i]]))
    writeLines(paste0("compare protein domains in annotated isoform ",transcript_names[[1]]," with those in predicted isoform ",transcript_names[[i]]),log_file)
    domains_i <- unique(as.character(domains_per_isoform[[i]]$V13))
    
    common <- as.character(domains_annotated_isoform[domains_annotated_isoform %in% domains_i])
    only_annot <- as.character(domains_annotated_isoform[!(domains_annotated_isoform %in% common)])
    only_i <- domains_i[!(domains_i %in% common)]
    
    if(length(only_annot) > 0){
      print(paste0(transcript_names[[i]]," misses the domain(s):"))
      print(only_annot)
      writeLines(paste0(transcript_names[[i]]," misses the domain(s):"),log_file)
      writeLines(only_annot,log_file)
    }
    if(length(only_i) > 0){
      print(paste0(transcript_names[[i]]," has additional domain(s):"))
      print(only_i)
      writeLines(paste0(transcript_names[[i]]," has additional domain(s):"),log_file)
      writeLines(only_i,log_file)
    }
    if(length(only_annot) == 0 & length(only_i) == 0){
      print("both include the same protein domains")
      writeLines("both include the same protein domains",log_file)
    } 
  }else{
    print(paste0(transcript_names[[i]]," does not have a valid amino acid sequence and thus no protein domains"))
    writeLines(paste0(transcript_names[[i]]," does not have a valid amino acid sequence and thus no protein domains"),log_file)
  }
}
close(log_file)

##TODO: compare all isoforms pairwise (there might be more than one annotated isoform in fungal annotations?!)


