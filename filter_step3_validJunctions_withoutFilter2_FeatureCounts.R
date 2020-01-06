library(Rsubread)
library(rtracklayer)
library(GenomicAlignments)

#install.packages('multicore',,'http://www.rforge.net/')
#library(multicore)
library(parallel)


location_as_analysis <- "/home/patricia/run_filter/A_fumigatus/KatrinLapp/"
location_bams <- "/home/patricia/run_filter/A_fumigatus/KatrinLapp/bams/"

setwd(location_as_analysis)

## import all as events #####
annot <- import.gff3("/home/patricia/run_filter/A_fumigatus/Aspergillus_fumigatusa1163.CADRE.31.gff3")
fasta_file <- "/home/patricia/run_filter/A_fumigatus/Aspergillus_fumigatusa1163.CADRE.31.dna.toplevel.fa"   ## still used later..
fasta <- import(fasta_file,format="fasta",type="DNA")

species <- "Afumigatus"


import.as.events <- function(condition){
  print(paste0(condition,":"))
  ## DEXSeq ####
  current_as <- NULL
  path_to_dexseq <- paste("DEXSeq/DEXSeq_",condition,"/results_",condition,".txt",sep="")
  if(file.exists(path_to_dexseq)){
    dexseq <- read.table(path_to_dexseq,sep="\t",header=TRUE)
    dexseq <- dexseq[dexseq$padj < 0.01,]
    dexseq <- unique(dexseq$groupID)
    print(paste("AS events DEXSeq: ",length(dexseq),sep=""))
    dexseq <- as.character(dexseq)
    
    ## filter 1: remove DEXSeq merge genes ####
    plus <- grep("[+]",dexseq)
    if(length(plus) > 0) dexseq <- dexseq[-plus]   
    current_as <- c(current_as,as.character(dexseq))
    print(paste("after removeing DEXSeq merged genes: ",length(current_as),sep=""))
  }
  else print("AS events DEXSeq: 0")  
  
  ## DiffSplice ####
  path_to_diffsplice <- paste("DiffSplice/DiffSplice2_",species,"_",condition,"/result/differential_transcription_filtered_mapped2genes.txt",sep="")
  if(file.exists(path_to_diffsplice)){
    diffsplice <- read.table(path_to_diffsplice,header=TRUE)
    diffsplice <- diffsplice[[1]]
    print(paste("AS events DiffSplice: ",length(diffsplice),sep=""))
    current_as <- c(current_as,as.character(diffsplice))
  }
  else print("AS events DiffSplice: 0")
  
  ## MATS ####
  import_mats_output <- function(path_to_mats_output,fdr=0.01){
    mats_files <- list.files(path_to_mats_output)
    mats_table <- NULL
    for(i in mats_files){
      current_table <- read.table(paste(path_to_mats_output,i,sep=""),header=TRUE)
      sub_table <- current_table[c("GeneID","chr","strand","PValue","FDR","IncLevel1","IncLevel2")]
      mats_table <- rbind(mats_table,sub_table)
    }
    mats_table <- mats_table[mats_table$FDR < fdr,]
    return(mats_table)
  }
  
  mats <- import_mats_output(paste("MATS/MATS_",species,"_",condition,"/MATS_output/",sep=""),0.01)
  mats <- mats$GeneID
  print(paste("AS events MATS: ",length(mats),sep=""))
  current_as <- c(current_as,as.character(mats))
  
  ## MISO ####
  read_length <- 51
  conds <- strsplit(as.character(condition),"vs",fixed=TRUE)
  cond_miso <- paste(species,"_",conds[[1]][1],".l",read_length,"_vs_",species,"_",conds[[1]][2],".l",read_length,sep="")
  path_to_miso <- paste("MISO/FASTMISO/comparisons/",cond_miso,"/bayes-factors/",cond_miso,".filtered.miso_bf",sep="")
  if(file.exists(path_to_miso)){
    miso <- read.table(path_to_miso,sep="\t",header=TRUE)
    ## already filtered by bayes factor, no additional filter necessary..
    miso <- miso$event_name
    print(paste("AS events MISO: ",length(miso),sep=""))
    current_as <- c(current_as,as.character(miso))  ##  
  }
  else print("AS events MISO: 0")
  
  ## SplicingCompass ####  
  path_to_splicingcompass <- paste("SplicingCompass/",condition,"/SplicingCompass_sigGenes_",condition,".txt",sep="")
  if(file.exists(path_to_splicingcompass)){
    splicingcompass <- read.table(path_to_splicingcompass,sep="\t",header=TRUE)
    
    splicingcompass <- splicingcompass[splicingcompass$tTestPairWiseAnglesPValHochberg < 0.01,]
    
    splicingcompass <- splicingcompass$gene_id
    splicingcompass <- gsub("mRNA","gene",splicingcompass)
    print(paste("AS events SplicingCompass: ",length(splicingcompass),sep=""))
    current_as <- c(current_as,as.character(splicingcompass))
  }
  else print("AS events SplicingCompass: 0")
  
  ## output list of all as events (unique)
  current_as <- unique(current_as)
  print(paste("all unique AS events: ",length(current_as),sep=""))
  return(current_as)
}

comparisons <- c("C30vsT30","C90vsT90","C180vsT180","T30vsT0","T90vsT0","T180vsT0","C30vsT0","C90vsT0","C180vsT0")

as_all <- lapply(comparisons,function(x){
  import.as.events(x)
})


as_all <- lapply(as_all,function(x){
  unlist(lapply(x,function(y){
    splittd <- strsplit(as.character(y),":")[[1]]
    len <- length(splittd)
    return(splittd[len])
  }))
})



exons <- annot[annot$type == "exon",]
parents <- as.character(exons$Parent)
exons$transcript_id <- unlist(lapply(parents,function(x){strsplit(as.character(x),":")[[1]][2]}))


transcripts <- annot[annot$type == "transcript",]
parents <- as.character(transcripts$Parent)
gene <- unlist(lapply(parents,function(x){strsplit(as.character(x),":")[[1]][2]}))
transcripts$gene_id <- gene

gene_ids <- lapply(exons$transcript_id,function(x){
  id <- transcripts[transcripts$transcript_id == as.character(x),]$gene_id
  if(length(id) == 0){return("None")
  }else{return(id)}
})
exons$gene_id <- unlist(gene_ids)
exons$gene <- unlist(gene_ids)


rm(list=ls(pattern="exons_as_all"))
gc()

exons_as_all <- lapply(as_all,function(x){
  exons_as <- lapply(x,function(y){
    position <- exons[exons$gene == y,]
    return(position)
  })
  exons_as <- lapply(exons_as,function(y){
    new_exons <- GRanges(seqnames=seqnames(y),ranges=ranges(y),strand=strand(y),source=y$source,type=y$type,ID=y$ID,Name=y$ID,
                         Parent=y$transcript_id,gene=y$gene)  ## changed from Parent=y$Parent
    return(new_exons)
  })
  return(exons_as)
})


## additionally to remove NAs from ID column
exons_as_all <- lapply(exons_as_all,function(x){
  n <- length(x)
  if(n > 0){
    for(i in 1:n){
      #print(i)
      current_gene <- x[[i]]
      if(length(current_gene) > 0){
        transcripts <- unlist(current_gene$Parent)
        n_exons <- length(transcripts)
        exon_id <- paste(transcripts,"-E",sep="")
        if(as.character(strand(current_gene)[1]) == "+") exon_id <- paste0(exon_id,1:n_exons)
        else exon_id <- paste0(exon_id,n_exons:1)
        x[[i]]$ID <- exon_id
        x[[i]]$Name <- exon_id
      }
    }
  }
  return(x)
})

## funtion define_introns: 
## input exons as GRanges object
define_introns <- function(exons_as){
  n <- length(exons_as)
  introns <- GRanges()
  exonsAndIntrons <- GRanges()
  
  for(i in 1:n){
    current_gene <- exons_as[[i]]
    ## split for each transcript, check separately
    transcripts <- unique(unlist(current_gene$Parent))
    for(j in transcripts){
      current_transcript <- current_gene[as.character(current_gene$Parent) == j,]
      
      ## don't consider transcripts without introns (just one exon..)
      current_intron <- NULL
      
      if(length(current_transcript) >= 2){
        current_intron <- gaps(current_transcript,start=min(start(current_transcript)))
      }
      
      n_intron <- length(current_intron)
      if(n_intron > 0){
        introns <- c(introns,current_intron)
        
        intron_id <- paste(j,"-I",sep="")
        if(as.character(strand(introns)[1]) == "+") intron_id <- paste0(intron_id,1:n_intron)
        else intron_id <- paste0(intron_id,n_intron:1)
        
        current_intron_ext <- GRanges(seqnames = seqnames(current_intron),ranges=ranges(current_intron),
                                      strand=strand(current_intron),source=factor(rep("calculated",n_intron)),
                                      type=factor(rep("intron",n_intron)),ID=intron_id,Name=intron_id,
                                      Parent=rep(as.character(j),n_intron),gene=as.character(rep(current_gene$gene[1],n_intron)))          ## changed from Parent=rep(CharacterList(j),n_intron)  
        
        current_exonsIntrons <- c(current_transcript,current_intron_ext)
        exonsAndIntrons <- c(exonsAndIntrons,current_exonsIntrons)  
        ## to have a collection of all genes
      }
      else{  ## no introns but its AS anyway
        exonsAndIntrons <- c(exonsAndIntrons,current_transcript)
      }
    }
  }
  return(exonsAndIntrons)
}

exonsAndIntrons_as_all <- mclapply(exons_as_all,function(x){
  if(length(x) > 0){
    exonsAndIntrons <- define_introns(x)
    print(paste("length:",length(exonsAndIntrons)))
    table(exonsAndIntrons$type) 
    print(paste("length_unique:",length(unique(exonsAndIntrons$gene))))
    return(exonsAndIntrons)
  }else return(list())
},mc.cores=4,mc.preschedule=FALSE)



setwd(location_bams)
filenames <- list.files(pattern="+.sorted.bam")
filenames <- filenames[-grep("bai",filenames)]

exonsAndIntrons_as_all_featureCounts <- lapply(exonsAndIntrons_as_all,function(x){
  x$id <- x$ID
  return(x)
})

annot_for_featureCounts <- lapply(exonsAndIntrons_as_all_featureCounts,function(x){
  if(length(x) > 0) return(createAnnotationFile(x))
  else return(NULL)
})

j <- 1
annot_for_featureCounts_all <- NULL
while(is.null(annot_for_featureCounts_all)){
  annot_for_featureCounts_all <- annot_for_featureCounts[[j]]
  j <- j+1
}
##Fehler hier drin..
while(j <= length(annot_for_featureCounts)){
  if(length(annot_for_featureCounts[[j]]) > 0) annot_for_featureCounts_all <- rbind(annot_for_featureCounts_all,annot_for_featureCounts[[j]])
  j <- j+1
}


objectnames <- c("bam_Af_D1C180","bam_Af_D1C30","bam_Af_D1C90","bam_Af_D1T0","bam_Af_D1T180","bam_Af_D1T30","bam_Af_D1T90",
                 "bam_Af_D2T180","bam_Af_D4C180","bam_Af_D4C30","bam_Af_D4C90","bam_Af_D4T0","bam_Af_D4T180","bam_Af_D4T30","bam_Af_D4T90",
                 "bam_Af_D5C180","bam_Af_D5C30","bam_Af_D5C90","bam_Af_D5T0","bam_Af_D5T30","bam_Af_D5T90",
                 "bam_Af_D6C180","bam_Af_D6C30","bam_Af_D6C90","bam_Af_D6T0","bam_Af_D6T180","bam_Af_D6T30","bam_Af_D6T90")

junctionnames <- c("junc_Af_D1C180","junc_Af_D1C30","junc_Af_D1C90","junc_Af_D1T0","junc_Af_D1T180","junc_Af_D1T30","junc_Af_D1T90",
                   "junc_Af_D2T180","junc_Af_D4C180","junc_Af_D4C30","junc_Af_D4C90","junc_Af_D4T0","junc_Af_D4T180","junc_Af_D4T30","junc_Af_D4T90",
                   "junc_Af_D5C180","junc_Af_D5C30","junc_Af_D5C90","junc_Af_D5T0","junc_Af_D5T30","junc_Af_D5T90",
                   "junc_Af_D6C180","junc_Af_D6C30","junc_Af_D6C90","junc_Af_D6T0","junc_Af_D6T180","junc_Af_D6T30","junc_Af_D6T90")

len_files <- length(filenames)

paired_end <- TRUE
for(i in 1:len_files){
  assign(objectnames[i],featureCounts(filenames[i],annot.ext=annot_for_featureCounts_all,isPairedEnd=paired_end))
}


rm(list=ls(pattern="annot"))
rm(list=ls(pattern="exons_as_all"))
rm(list=ls(pattern="exonsAndIntrons_as_all_featureCounts"))
rm(list=ls(pattern="annot_for_featureCounts_all"))
gc()

setwd(location_as_analysis)



## exonsAndIntrons: object; replicates: names of bam replicate objects in a vector
count_reads <- function(exonsAndIntrons,replicates_condition1,replicates_condition2){
  len_exIn <- length(exonsAndIntrons)
  exonsAndIntrons$count1 <- rep(0,len_exIn)  ## calculate read count value, condition1
  exonsAndIntrons$count1_sd <- rep(0,len_exIn)
  exonsAndIntrons$count2 <- rep(0,len_exIn)  ## calculate read count value, condition2
  exonsAndIntrons$count2_sd <- rep(0,len_exIn)
  
  for(i in 1:len_exIn){
    current_element <- as.character(exonsAndIntrons[i]$ID)
    
    counter <- lapply(replicates_condition1,function(z){
      current_replicate <- get(z)$counts
      sum <- current_replicate[current_element,]
      return(sum)
    })
    exonsAndIntrons[i]$count1 <- mean(unlist(counter))
    exonsAndIntrons[i]$count1_sd <- sd(unlist(counter))
    
    counter <- lapply(replicates_condition2,function(z){
      current_replicate <- get(z)$counts
      sum <- current_replicate[current_element,]
      return(sum)
    })
    exonsAndIntrons[i]$count2 <- mean(unlist(counter))
    exonsAndIntrons[i]$count2_sd <- sd(unlist(counter))
  }
  return(exonsAndIntrons)
}

exonsAndIntrons_as_all <- mcmapply(function(x,y){
  if(length(x) > 0){
    current_conds <- strsplit(as.character(y),"vs",fixed=TRUE)
    cond1 <- paste0(current_conds[[1]][1])
    cond2 <- paste0(current_conds[[1]][2])
    replicates_condition1 <- objectnames[grep(cond1,objectnames,fixed=TRUE)]
    replicates_condition2 <- objectnames[grep(cond2,objectnames,fixed=TRUE)]
    x <- count_reads(x,replicates_condition1,replicates_condition2)
  }
  return(x)
},exonsAndIntrons_as_all,comparisons,mc.cores=4,mc.preschedule=FALSE)

test_exonsAndIntrons_as_all <- exonsAndIntrons_as_all

## calculate rpkm values (Mortazavi2008) for each exon/intron
## rpkm = (10^9 * number of mapped reads)/(total number of mapped reads * length)
## exonsAndIntrons: object with count1 and count 2 values
rpkm <- function(exonsAndIntrons,replicates_condition1,replicates_condition2){
  len_exIn <- length(exonsAndIntrons)
  exonsAndIntrons$rpkm1 <- rep(0,len_exIn)  ## calculate rpkm value, condition1
  exonsAndIntrons$rpkm2 <- rep(0,len_exIn)  ## calculate rpkm value, condition2
  
  ## total read count of all replicates of one condition:
  count <- unlist(lapply(replicates_condition1,function(z){
    current_replicate <- get(z)$counts
    counter <- sum(current_replicate)
    return(counter)
  }))
  total_count_cond1 <- as.numeric(sum(count))
  
  count <- unlist(lapply(replicates_condition2,function(z){
    current_replicate <- get(z)$counts
    counter <- sum(current_replicate)
    return(counter)
  }))
  total_count_cond2 <- as.numeric(sum(count))
  
  for(i in 1:len_exIn){
    current_element <- exonsAndIntrons[i]
    exonsAndIntrons[i]$rpkm1 <- (10^9 * as.numeric(current_element$count1))/(total_count_cond1 * as.numeric(width(current_element)))
    exonsAndIntrons[i]$rpkm2 <- (10^9 * as.numeric(current_element$count2))/(total_count_cond2 * as.numeric(width(current_element)))
  }
  return(exonsAndIntrons)
}

exonsAndIntrons_as_all <- mcmapply(function(x,y){
  if(length(x) > 0){
    current_conds <- strsplit(as.character(y),"vs",fixed=TRUE)
    cond1 <- paste0(current_conds[[1]][1])
    cond2 <- paste0(current_conds[[1]][2])
    replicates_condition1 <- objectnames[grep(cond1,objectnames,fixed=TRUE)]
    replicates_condition2 <- objectnames[grep(cond2,objectnames,fixed=TRUE)]
    x <- rpkm(x,replicates_condition1,replicates_condition2) 
  }
  return(x)
},exonsAndIntrons_as_all,comparisons,mc.cores=4,mc.preschedule=FALSE)

rm(list=ls(pattern="bam_"))
gc()

test_exonsAndIntrons_as_all <- exonsAndIntrons_as_all


## filter3: exclude low read count ####
## exclude gene if read count of exons is lower than 10 reads (on average)
exclude_low_read_count <- function(exonsAndIntrons,threshold=20){
  ## if the coverage of all elements of a transcripts is lower than the threshold, then remove the transcript
  parents <- unlist(lapply(exonsAndIntrons$Parent,function(z){as.character(z[[1]])}))
  transcripts <- unique(parents)
  len_tr <- length(transcripts)
  exclude_transcript <- NULL ## collect transcripts that will be excluded
  for(i in 1:len_tr){
    #print(i)
    cat(i)
    current_transcript <- exonsAndIntrons[parents == transcripts[i],]
    current_transcript <- current_transcript[current_transcript$type == "exon",]  ## just exons
    if(length(current_transcript) > 0){
      if(mean(current_transcript$count1) <= threshold || mean(current_transcript$count2) <= threshold){
        exclude_transcript <- c(exclude_transcript,current_transcript$gene[1])
      }
    }
  }
  exonsAndIntrons_new <- unlist(lapply(exonsAndIntrons,function(k){
    match_exclude_transcript <- which(exclude_transcript == k$gene)
    if(length(match_exclude_transcript) > 0) return(NULL)
    else return(k)
  }))
  if(!is.null(exonsAndIntrons_new)) exonsAndIntrons_new <- do.call("c", exonsAndIntrons_new)
  return(exonsAndIntrons_new)
}

exonsAndIntrons_as_all <- mcmapply(function(x,y){
  if(length(x) > 0){
    cat(y)
    x <- exclude_low_read_count(x)
    print(paste0("unique genes ",y,":",length(unique(x$gene)))) 
  }
  return(x)
},exonsAndIntrons_as_all,comparisons,mc.cores=4,mc.preschedule=FALSE)

lapply(exonsAndIntrons_as_all,function(x){return(unique(x$gene))})
test_exonsAndIntrons_as_all <- exonsAndIntrons_as_all


## filter4: exclude low rpkm ####
## exclude by rpkm coverage
# at least one exon with rpkm > 5 (Aschoff2013, Trapnell2009), otherwise noise
##--> too hard, would exclude several identified AS events!, use 1
exclude_low_rpkm <- function(exonsAndIntrons,threshold=1){
  ## if the coverage of all elements of a transcripts is lower than the threshold, then remove the transcript
  parents <- unlist(lapply(exonsAndIntrons$Parent,function(z){as.character(z[[1]])}))
  transcripts <- unique(parents)
  len_tr <- length(transcripts)
  exclude_transcript <- NULL ## collect transcripts that will be excluded
  for(i in 1:len_tr){
    current_transcript <- exonsAndIntrons[as.character(exonsAndIntrons$Parent) == transcripts[i],]
    current_transcript <- current_transcript[current_transcript$type == "exon",]  ## just exons
    if(mean(current_transcript$rpkm1) <= threshold || mean(current_transcript$rpkm2) <= threshold) 
      exclude_transcript <- c(exclude_transcript,current_transcript$gene[1])
  }
  exonsAndIntrons_new <- unlist(lapply(exonsAndIntrons,function(k){
    match_exclude_transcript <- which(exclude_transcript == k$gene)
    if(length(match_exclude_transcript) > 0) return(NULL)
    else return(k)
  }))
  if(!is.null(exonsAndIntrons_new)) exonsAndIntrons_new <- do.call("c", exonsAndIntrons_new)
  return(exonsAndIntrons_new)
}

exonsAndIntrons_as_all <- mcmapply(function(x,y){
  if(length(x) > 0){
    x <- exclude_low_rpkm(x)
    print(paste0("unique genes ",y,":",length(unique(x$gene))))
  }
  return(x)
},exonsAndIntrons_as_all,comparisons,mc.cores=4,mc.preschedule=FALSE)

lapply(exonsAndIntrons_as_all,function(x){return(unique(x$gene))})
test_exonsAndIntrons_as_all <- exonsAndIntrons_as_all


setwd(location_bams)
if(paired_end){  ## for paired end data
  for(i in 1:len_files){
    assign(objectnames[i], readGAlignmentPairs(filenames[i]))
    assign(junctionnames[i],summarizeJunctions(get(objectnames[i])))
  }  
}else{  ## for single end data
  for(i in 1:len_files){
    assign(objectnames[i], readGAlignments(filenames[i]))
    assign(junctionnames[i],summarizeJunctions(get(objectnames[i])))
  }
} 

## just junctions are needed for further steps.
rm(list=ls(pattern="bam_"))
gc()

setwd(location_as_analysis)

## count_junctions
## exonsAndIntrons: object; replicates: names of junction replicate objects in a vector
count_junctions <- function(exonsAndIntrons,replicates_condition1,replicates_condition2){
  len_exIn <- length(exonsAndIntrons)
  nr_replicates <- length(replicates_condition1)
  exonsAndIntrons$junctions1 <- rep(0,len_exIn)  ## add number of spanning junctions, condition1
  exonsAndIntrons$junctions1_sd <- rep(0,len_exIn)  
  exonsAndIntrons$junctions2 <- rep(0,len_exIn)  ## add number of spanning junctions, condition2
  exonsAndIntrons$junctions2_sd <- rep(0,len_exIn) 
  
  introns <- exonsAndIntrons[exonsAndIntrons$type == "intron",]
  len_intr <- length(introns)
  for(i in 1:len_intr){
    current_intron <- introns[i]
    strand_info <- as.character(strand(current_intron))
    
    sum_up <- vector(mode="numeric",length=nr_replicates)
    for(j in 1:nr_replicates){
      current_replicate <- get(replicates_condition1[j])
      ov <- findOverlaps(current_intron,current_replicate)
      current_junction <- current_replicate[subjectHits(ov),]
      current_junction <- current_junction[ranges(current_junction) == ranges(current_intron),]
      
      if(length(current_junction) == 1){
        if(strand_info == "+") sum_up[j] <- current_junction$plus_score
        else{
          if(strand_info == "-") sum_up[j] <- current_junction$minus_score
          else sum_up[j] <- current_junction$score  ## strand == "*", so no information about strand
        }
      }
      if(length(current_junction) >= 2) writeLines(paste(i,",",j))
    }
    exonsAndIntrons[exonsAndIntrons$ID == as.character(current_intron$ID),]$junctions1 <- mean(sum_up)
    exonsAndIntrons[exonsAndIntrons$ID == as.character(current_intron$ID),]$junctions1_sd <- sd(sum_up)
    
    sum_up <- vector(mode="numeric",length=nr_replicates)
    for(j in 1:nr_replicates){
      current_replicate <- get(replicates_condition2[j])
      ov <- findOverlaps(current_intron,current_replicate)
      current_junction <- current_replicate[subjectHits(ov),]
      current_junction <- current_junction[ranges(current_junction) == ranges(current_intron),]
      
      if(length(current_junction) == 1){
        if(strand_info == "+") sum_up[j] <- current_junction$plus_score
        else{
          if(strand_info == "-") sum_up[j] <- current_junction$minus_score
          else sum_up[j] <- current_junction$score  ## strand == "*", so no information about strand
        }
      }
      if(length(current_junction) >= 2) writeLines(paste(i,",",j))
    }
    exonsAndIntrons[exonsAndIntrons$ID == as.character(current_intron$ID),]$junctions2 <- mean(sum_up)
    exonsAndIntrons[exonsAndIntrons$ID == as.character(current_intron$ID),]$junctions2_sd <- sd(sum_up)
  }
  return(exonsAndIntrons)
}
exonsAndIntrons_as_all <- mcmapply(function(x,y){
  if(length(x) > 0){
    current_conds <- strsplit(as.character(y),"vs",fixed=TRUE)
    cond1 <- paste0(current_conds[[1]][1])
    cond2 <- paste0(current_conds[[1]][2])
    replicates_condition1 <- junctionnames[grep(cond1,junctionnames,fixed=TRUE)]
    replicates_condition2 <- junctionnames[grep(cond2,junctionnames,fixed=TRUE)]
    x <- count_junctions(x,replicates_condition1,replicates_condition2)
  }
  return(x)
},exonsAndIntrons_as_all,comparisons,mc.cores=4,mc.preschedule=FALSE)

test_exonsAndIntrons_as_all <- exonsAndIntrons_as_all


## in every pattern of alternative splicing, the intron junction changes
## thus, another portion of junctions may indicate another splicing pattern
#install.packages("gtools", dependencies=TRUE, repos='http://cran.rstudio.com/')
library(gtools)
logFoldJunction <- function(exonsAndIntrons){
  len_exIn <- length(exonsAndIntrons)
  exonsAndIntrons$logFoldChange_junctions <- rep(0,len_exIn)
  introns <- exonsAndIntrons[exonsAndIntrons$type == "intron",]
  len_intr <- length(introns)
  for(i in 1:len_intr){
    current_intron <- introns[i]
    fchange <- foldchange(current_intron$junctions1,current_intron$junctions2)
    logfchange <- foldchange2logratio(fchange,base=2)
    exonsAndIntrons[exonsAndIntrons$ID == as.character(current_intron$ID),]$logFoldChange_junctions <- logfchange
  }
  return(exonsAndIntrons)
}

exonsAndIntrons_as_all <- mcmapply(function(x,y){
  if(length(x) > 0){
    x <- logFoldJunction(x)
  }
  return(x)
},exonsAndIntrons_as_all,comparisons,mc.cores=4,mc.preschedule=FALSE)

test_exonsAndIntrons_as_all <- exonsAndIntrons_as_all

#### filter5: exclude transcripts without any change in exon count / junction count ####
#source("/home/patricia/run_filter/filter_exclusion_excludeInvalicJunctions.R")
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome")
source("filter_exclusion_excludeInvalicJunctions.R")

exonsAndIntrons_as_all <- mcmapply(function(x,y){
  if(length(x) > 0){
    current_conds <- strsplit(as.character(y),"vs",fixed=TRUE)
    cond1 <- paste0(current_conds[[1]][1])
    cond2 <- paste0(current_conds[[1]][2])
    replicates_condition1 <- junctionnames[grep(cond1,junctionnames,fixed=TRUE)]
    replicates_condition2 <- junctionnames[grep(cond2,junctionnames,fixed=TRUE)]
    x <- exclude_no_change(x,replicates_condition1,replicates_condition2,fasta=fasta)
    print(paste0("unique genes ",y,":",length(unique(x$gene))))
  }
  return(x)
},exonsAndIntrons_as_all,comparisons,mc.cores=4,mc.preschedule=FALSE)

lapply(exonsAndIntrons_as_all,function(x){return(unique(x$gene))})
#test_exonsAndIntrons_as_all <- exonsAndIntrons_as_all


test <- mapply(function(x,y){
  if(length(x) > 0){
    print(paste0("unique genes ",y,":"))
    print(unique(x$gene))
    current <- paste0("exonsAndIntrons_as_withoutFilter2_FeatureCounts_",y)
    ffile <- file(paste0(current,"_filtered.txt"),"w")
    writeLines(paste("seqnames","start","end","strand","type","ID","gene","junctions1","junctions1_sd",
                     "junctions2","junctions2_sd","count1","count1_sd","count2","count2_sd",
                     "rpkm1","rpkm2","logFoldChange_junctions",sep="\t"),ffile)
    writeLines(paste(as.character(seqnames(x)),as.character(start(x)),
                     as.character(end(x)),as.character(strand(x)),
                     as.character(x$type),as.character(x$ID),
                     as.character(x$gene),
                     as.character(x$junctions1),as.character(x$junctions1_sd),
                     as.character(x$junctions2),as.character(x$junctions2_sd),
                     as.character(x$count1),as.character(x$count1_sd),
                     as.character(x$count2),as.character(x$count2_sd),
                     as.character(x$rpkm1),as.character(x$rpkm2),
                     as.character(x$logFoldChange_junctions),
                     sep="\t"),ffile)
    close(ffile)
  }
},exonsAndIntrons_as_all,comparisons)
