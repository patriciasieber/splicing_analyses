## statistics of gene annotations

setwd("/home/bioinformatik/Hacken2017/")
#gtf <- read.csv("annots_original/Homo_sapiens.GRCh38.85.chr.gtf",header=F,skip=5,sep="\t")
gtf <- read.csv("annots_original/Danio_rerio.GRCz10.release83.annotation.gtf",header=F,skip=5,sep="\t")

table(gtf$V3)

gtf_modified <- gtf[gtf$V3 != "Selenocysteine",]

ffile <- file("annots_mats/Danio_rerio.GRCz10.release83.annotation.modified.gtf","w")
writeLines(paste(gtf_modified$V1,gtf_modified$V2,gtf_modified$V3,gtf_modified$V4,gtf_modified$V5,
                 gtf_modified$V6,gtf_modified$V7,gtf_modified$V8,gtf_modified$V9,sep="\t"),ffile)
close(ffile)




## for mmu
setwd("/home/bioinformatik/Hacken2017/")
gtf <- read.csv("annots_original/Mus_musculus.GRCm38.87.chr.sorted.chr.gtf",header=F,skip=5,sep="\t")
table(gtf$V3)

gtf_modified <- gtf[gtf$V3 != "Selenocysteine",]

#chr_new <- as.character(gtf_modified$V1)
#chr_new <- paste("chr",chr_new,sep="")

ffile <- file("annots_mats/Mus_musculus.GRCm38.87.chr.sorted.chr.modified.gtf","w")
writeLines(paste(gtf_modified$V1,gtf_modified$V2,gtf_modified$V3,gtf_modified$V4,gtf_modified$V5,
                 gtf_modified$V6,gtf_modified$V7,gtf_modified$V8,gtf_modified$V9,sep="\t"),ffile)
close(ffile)




## for nfu:

gtf <- read.csv("annots_original/nothobranchius_furzeri_annotation.gff",header=F,skip=5,sep="\t")
table(gtf$V3)

gtf_modified <- gtf[gtf$V3 != "pre_miRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "mature_miRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "lncRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "ncRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "ribozyme",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "rRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "snoRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "snRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "sRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "tRNA",]


ffile <- file("nothobranchius_furzeri_annotation.modified.gtf","w")
writeLines(paste(gtf_modified$V1,gtf_modified$V2,gtf_modified$V3,gtf_modified$V4,gtf_modified$V5,
                 gtf_modified$V6,gtf_modified$V7,gtf_modified$V6,gtf_modified$V9,sep="\t"),ffile)
close(ffile)



gtf <- read.csv("annots_original/nothobranchius_furzeri_annotation.gff",header=F,skip=5,sep="\t")
table(gtf$V3)

gtf_modified <- gtf[gtf$V3 != "pre_miRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "mature_miRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "lncRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "ncRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "ribozyme",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "rRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "snoRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "snRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "sRNA",]
gtf_modified <- gtf_modified[gtf_modified$V3 != "tRNA",]

table(gtf_modified$V3)


ffile <- file("nothobranchius_furzeri_annotation.modified.gff3","w")
writeLines(paste(gtf_modified$V1,gtf_modified$V2,gtf_modified$V3,gtf_modified$V4,gtf_modified$V5,
                 gtf_modified$V6,gtf_modified$V7,gtf_modified$V6,gtf_modified$V9,sep="\t"),ffile)
close(ffile)

## prak27@dessert:~$ /mnt/prostlocal/programs/cufflinks/2.2.1/gffread nothobranchius_furzeri_annotation.modified.gff3 -T -o nothobranchius_furzeri_annotation.modified.gtf

setwd("/home/bioinformatik/Hacken2017/")
gtf <- read.csv("annots_original/Homo_sapiens.GRCh38.85.chr.sorted.nochr.gtf",header=F,skip=5,sep="\t")

table(gtf$V3)

gtf_modified <- gtf[gtf$V3 != "Selenocysteine",]

ffile <- file("annots_mats/Homo_sapiens.GRCh38.85.chr.sorted.nochr.modified.gtf","w")
writeLines(paste(gtf_modified$V1,gtf_modified$V2,gtf_modified$V3,gtf_modified$V4,gtf_modified$V5,
                 gtf_modified$V6,gtf_modified$V7,gtf_modified$V8,gtf_modified$V9,sep="\t"),ffile)
close(ffile)






# source("http://bioconductor.org/biocLite.R")
# biocLite("rtracklayer")

#library(rtracklayer)

#setwd("/media/patricia/UDISK/annotation_stats/L.corymbifera/")

## import gene annotation ####
## gff3 file
file_name <- "Mus_musculus.GRCm38.87.chr.sorted.gtf"


gff <- import.gff(file_name)

mrna <- gff[gff$type == "mRNA",]
gene <- unlist(lapply(mrna$Parent,function(x){
  return(as.character(x[[1]]))
}))
mrna$gene <- gene
mrna$transcript_id <- mrna$ID

exons <- gff[gff$type == "exon",]
tr <- unlist(lapply(as.character(exons$Parent),function(x){
  return(current <- as.character(x[[1]]))
}))
exons$transcript_id <- tr

gene_id <- unlist(lapply(exons$transcript_id,function(x){
  current <- mrna[mrna$transcript_id == as.character(x),]$gene
  if(length(current) == 1) return(current)
  else return("none")
}))
exons$gene <- gene_id


## exclude non-mrna genes
mrna_genes <- unique(mrna$gene)
mrna_exons <- lapply(mrna_genes,function(x){
  return(exons[exons$gene == x,])
})
mrna_exons <- do.call(c,mrna_exons) 

transcript_id <- as.character(mrna_exons$Parent)
mrna_exons$transcript_id <- transcript_id

## annotation stats ####
ffile <- file(paste(file_name,"_stats",sep=""),"w")
writeLines("stats for all annotated transcripts",ffile)

genes <- unique(mrna_exons$gene)
writeLines(paste("number of protein-coding genes:",length(genes),sep="\t"),ffile)

transcript_ids <- unique(mrna_exons$transcript_id)
writeLines(paste("number of protein-coding transcripts:",length(transcript_ids),sep="\t"),ffile)

## define number of exons per transcript
writeLines(paste("number of protein-coding exons:",length(mrna_exons),sep="\t"),ffile)

exons_per_transcript <- table(mrna_exons$transcript_id)
ept <- data.frame(exons_per_transcript)

## define exon length
exon_length <- width(mrna_exons)
writeLines(paste("exon length:","mean:",mean(exon_length),"standard deviation:",sd(exon_length),"min:",min(exon_length),"max:",max(exon_length),sep="\t"),ffile)


## define intron length (min, mean, max)
## exclude genes without introns
intron_table <- ept[ept$Freq > 1,]
max_freq <- max(intron_table$Freq)
writeLines(paste("maximum intron number of a transcript:",max_freq,sep="\t"),ffile)
tr_with_introns <- as.character(intron_table$Var1)
writeLines(paste("number of transcripts with more than one exon (intron-containing):",length(tr_with_introns),sep="\t"),ffile)
proportion <- length(tr_with_introns)/length(transcript_ids) * 100
writeLines(paste("proportion of intron-containing transcripts:",paste(proportion,"%",sep=" "),sep="\t"),ffile)

genes_with_introns <- unlist(lapply(tr_with_introns,function(x){
  return(mrna[mrna$transcript_id == x]$gene)
}))
proportion <- length(genes_with_introns)/length(genes) * 100
writeLines(paste("proportion of intron-containing genes:",paste(proportion,"%",sep=" "),sep="\t"),ffile)

intron_numbers <- intron_table$Freq
number_of_introns <- sum(intron_numbers) - length(intron_numbers)  # for each transcript one intron less than exons
writeLines(paste("number of introns:",number_of_introns,sep="\t"),ffile)

exons_with_introns<- lapply(genes_with_introns,function(x){
  return(mrna_exons[mrna_exons$gene == x,])
})
exons_with_introns <- do.call(c,exons_with_introns)

## define range between exons of each transcript --> introns
all_introns <- lapply(tr_with_introns,function(x){
  current_exons <- exons_with_introns[exons_with_introns$transcript_id == as.character(x),]
  return(gaps(current_exons,start=min(start(current_exons))))
})
all_introns <- do.call(c, all_introns)

intron_length <- width(all_introns)
writeLines(paste("intron length:","mean:",mean(intron_length),"standard deviation:",sd(intron_length),"min:",min(intron_length),"max:",max(intron_length),sep="\t"),ffile)
#length(intron_length)

writeLines("------------------------------------------------------",ffile)

## define AS patterns ####
writeLines("stats for genes with more than one transcript",ffile)

gene_table <- table(mrna$gene)
ttable <- data.frame(gene_table)

multiple_tr <- ttable[ttable$Freq > 1,]
writeLines(paste("number of genes with multiple transcripts:",nrow(multiple_tr),sep="\t"),ffile)

proportion <- nrow(multiple_tr)/length(genes) * 100
writeLines(paste("proportion of genes with multiple transcripts:",paste(proportion,"%",sep=" "),sep="\t"),ffile)

multiple_exons <- lapply(multiple_tr$Var1,function(x){
  return(mrna_exons[mrna_exons$gene == as.character(x),])
})
multiple_exons <- do.call(c,multiple_exons)


writeLines("------------------------------------------------------",ffile)
if(length(multiple_exons) == 0){
  writeLines("no alternative splicing patterns",ffile)
}else{
  writeLines("alternative splicing patterns",ffile)
  
  multiple_transcripts <- unique(multiple_tr$Var1)  ## genes with multiple transcripts, used for AS pattern identification
  len_multtr <- length(multiple_transcripts)
  
  altStartSite <- 0
  altTerminationSite <- 0
  alt5End <- 0
  alt3End <- 0
  mxe <- 0
  intrRetention <- 0
  exonSkipping <- 0
  
  for(i in 1:len_multtr){
    current_gene <- as.character(multiple_transcripts[i])
    current_all <- multiple_exons[multiple_exons$gene == as.character(current_gene),]
    
    all_tr <- table(current_all$transcript_id)
    all_tr <- data.frame(all_tr)
    all_tr <- as.character(all_tr$Var1)
    
    current_transcripts <- lapply(all_tr,function(x){
      output <- sort(current_all[current_all$transcript_id == x,])
      ## if two exons have the identical range, remove one
      disjoined <- disjoin(output)
      if(length(disjoined) != length(output)){
        ov <- findOverlaps(output,disjoined)
        hits <- data.frame(table(subjectHits(ov)))  
        ## we can have just multiple maps when a full exon is present multiple times
        if(any(hits$Freq > 1)){
          mult_hit <- hits[hits$Freq > 1,]$Var1
          toremove <- NULL
          for(j in 1:length(mult_hit)){
            current_hit <- mult_hit[j]
            hitsQuery <- queryHits(ov[subjectHits(ov) == current_hit,])
            toremove <- c(toremove,hitsQuery[2:length(hitsQuery)])
          }
          output <- output[-c(toremove),]
        }
        ## if different exons overlap, ? necessary ?
      }
      return(output)
    })
    
    strand <- as.character(strand(current_transcripts[[1]])[1])
    
    ## unique transcripts
    k <- 1
    l <- 2
    while(k < length(current_transcripts) & length(current_transcripts) > 1){
      if(length(current_transcripts[[k]]) == length(current_transcripts[[l]])){
        current_tr1 <- current_transcripts[[k]]
        current_tr2 <- current_transcripts[[l]]
        if(all(start(current_tr1) == start(current_tr2) & end(current_tr1) == end(current_tr2))){
          ## then transcript is duplicated, remove one
          current_transcripts[[k]] <- NULL
        }
        else{
          if(l == length(current_transcripts)){
            k <- k+1
            l <- k+1
          }
          else l <- l+1
        }
      }
      else{
        if(l == length(current_transcripts)){
          k <- k+1
          l <- k+1
        }
        else l <- l+1
      }
    }
    
    ## compare transcripts against each other, reference transcript against all the others
    ## TODO which is the reference? for now the first
    if(length(current_transcripts) >= 2){
      for(j in 2:length(current_transcripts)){
        
        ## Alternative transcription start Site
        ## different first exon
        if(strand == "+"){
          if(min(start(current_transcripts[[1]])) != min(start(current_transcripts[[j]]))){
            altStartSite <- altStartSite + 1
            writeLines(paste(current_gene,"Alternative Start site",sep="\t"),ffile)
          }
        }
        if(strand == "-"){
          if(max(end(current_transcripts[[1]])) != max(end(current_transcripts[[j]]))){
            altStartSite <- altStartSite + 1
            writeLines(paste(current_gene,"Alternative Start site",sep="\t"),ffile)
          }
        }
        
        ## Alternative Termination Site
        ## different last exon
        if(strand == "+"){
          if(max(end(current_transcripts[[1]])) != max(end(current_transcripts[[j]]))){
            altTerminationSite <- altTerminationSite + 1
            writeLines(paste(current_gene,"Alternative Termination Site",sep="\t"),ffile)
          }
        }
        if(strand == "-"){
          if(min(start(current_transcripts[[1]])) != min(start(current_transcripts[[j]]))){
            altTerminationSite <- altTerminationSite + 1
            writeLines(paste(current_gene,"Alternative Termination Site",sep="\t"),ffile)
          }
        }
        
        ## now consider just overlapping exons (exons more up- or downstream are already identified and don't need to be considered in the following steps)
        ov_transcripts <- findOverlaps(current_transcripts[[1]],current_transcripts[[j]])
        
        subset_current_transcripts <- current_transcripts
        subset_current_transcripts[[1]] <- subset_current_transcripts[[1]][queryHits(ov_transcripts),]
        subset_current_transcripts[[1]] <- unique(subset_current_transcripts[[1]])
        subset_current_transcripts[[j]] <- subset_current_transcripts[[j]][subjectHits(ov_transcripts),]
        subset_current_transcripts[[j]] <- unique(subset_current_transcripts[[j]])
        
        if(length(subset_current_transcripts[[1]]) == length(subset_current_transcripts[[j]])){
          ## we can have A3E, A5E, MXE (or a combination of ES and IR)
          n <- length(subset_current_transcripts[[1]])
          disjoined <- disjoin(c(subset_current_transcripts[[1]],subset_current_transcripts[[j]]))
          
          ## no alternative splicing event for length(disjoined) == length(subset_current_transcripts[[1]])
          ## we need to consider just the opposite case:
          if(length(disjoined) != n){
            ## check range of each exon
            for(k in 1:n){
              exon1 <- subset_current_transcripts[[1]][k]
              exon2 <- subset_current_transcripts[[j]][k]
              
              if(start(exon1) != start(exon2) && end(exon1) != end(exon2)){
                min_whole_transcript <- min(start(subset_current_transcripts[[1]]),start(subset_current_transcripts[[j]]))
                max_whole_transcript <- max(end(subset_current_transcripts[[1]]),end(subset_current_transcripts[[j]]))
                
                if(min_whole_transcript == min(start(exon1),start(exon2))){
                  ## alternative start/end was already considered, so consider just end
                  if(end(exon1) != end(exon2)){
                    if(strand == "+"){
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                    if(strand == "-"){
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                  }
                }
                if((max_whole_transcript == max(end(exon1),end(exon2)))){
                  ## alternative end/start was already considered, so consider just end
                  if(start(exon1) != start(exon2)){
                    if(strand == "+"){
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                    if(strand == "-"){
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                  }
                }
                if(min_whole_transcript != min(start(exon1),start(exon2)) && (max_whole_transcript != max(end(exon1),end(exon2)))){
                  ## alternative start and termination was already considered
                  
                  ov <- findOverlaps(exon1,exon2)
                  if(length(ov) == 0){
                    ## another exon, considered as MXE
                    ## Mutually Exclusive Exon
                    ## each transcript isoform with a unique set of exons (with at least one exon per set), the rest stays unchanged
                    mxe <- mxe + 1
                    writeLines(paste(current_gene,"Mutually Exclusive Exon",sep="\t"),ffile)
                  }else{
                    ## overlapping exons, so we consider them as A3'E and A5'E
                    alt3End <- alt3End + 1
                    alt5End <- alt5End + 1
                    writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                  }
                }
              }else{
                if(start(exon1) != start(exon2)){
                  min_whole_transcript <- min(start(subset_current_transcripts[[1]]),start(subset_current_transcripts[[j]]))
                  if(strand == "+"){
                    if(min_whole_transcript != min(start(exon1),start(exon2))){
                      ## the opposite case was considered already as alternative transcription start
                      
                      ## Alternative 3'End
                      ## exon has an alternative start
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                  }
                  if(strand == "-"){
                    if(min_whole_transcript != min(start(exon1),start(exon2))){
                      ## the opposite case was considered already as alternative termination
                      
                      ## Alternative 5'End
                      ## exon has an alternative end
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                  }
                }
                if(end(exon1) != end(exon2)){
                  max_whole_transcript <- max(end(subset_current_transcripts[[1]]),end(subset_current_transcripts[[j]]))
                  if(strand == "+"){
                    if(max_whole_transcript != max(end(exon1),end(exon2))){
                      ## the opposite case was considered already as alternative termination
                      ## Alternative 5'End
                      ## exon has an alternative start
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    } 
                  }
                  if(strand == "-"){
                    if(max_whole_transcript != max(end(exon1),end(exon2))){
                      ## the opposite case was considered already as alternative transcription start
                      ## Alternative 3'End
                      ## exon has an alternative end
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                  }
                }
              }
              ## sonst noch was ?? 
            } ## Z. 252
          } ## Z.250
        } ## Z.243
        
        if(length(subset_current_transcripts[[1]]) != length(subset_current_transcripts[[j]])){
          ## we have ES or/and IR events (and maybe additionally A3E, A5E, MXE..)
          ## use a boolean vector to identify equal ranges of exons
          if(length(subset_current_transcripts[[1]]) > length(subset_current_transcripts[[j]])){
            longer_subset_current_transcripts <- subset_current_transcripts[[1]]
            shorter_subset_current_transcripts <- subset_current_transcripts[[j]]
          } else{
            shorter_subset_current_transcripts <- subset_current_transcripts[[1]]
            longer_subset_current_transcripts <- subset_current_transcripts[[j]]
          }
          n <- max(length(subset_current_transcripts[[1]]),length(subset_current_transcripts[[j]]))
          
          range_overlap <- vector(length=n)
          for(k in 1:n){
            exon <- longer_subset_current_transcripts[k]
            ov <- findOverlaps(exon,shorter_subset_current_transcripts)
            
            ## create vector to identify IR and ES
            if(length(ov) == 1) range_overlap[k] <- TRUE
            else range_overlap[k] <- FALSE
          }
          
          ## Exon Skipping
          ## an additional exon
          if(sum(range_overlap) == n-1){
            exonSkipping <- exonSkipping + 1
            writeLines(paste(current_gene,"Exon Skipping",sep="\t"),ffile)
            
            ## check for additional A3'E / A5'E
            ## check for common start and end just for overlapping exons
            ov <- findOverlaps(shorter_subset_current_transcripts,longer_subset_current_transcripts)
            n_short <- length(shorter_subset_current_transcripts)
            
            for(l in 1:n_short){
              exon1 <- shorter_subset_current_transcripts[l]
              hits <- ov[queryHits(ov) == l,]
              exons2 <- longer_subset_current_transcripts[subjectHits(hits)]
              
              if(length(exons2) == 1){
                if(start(exon1) != start(exons2) && end(exon1) != end(exons2)){
                  ## we know already that they are overlapping, so we have both A3'E and A5'E
                  alt3End <- alt3End + 1
                  alt5End <- alt5End + 1
                  writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                  writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                }
                else{
                  if(start(exon1) != start(exons2)){
                    min_whole_transcript <- min(start(subset_current_transcripts[[1]]),start(subset_current_transcripts[[j]]))
                    if(strand == "+" && min_whole_transcript != min(start(exon1),start(exons2))){
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                    if(strand == "-" && min_whole_transcript != min(start(exon1),start(exons2))){
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                  }
                  if(end(exon1) != end(exons2)){
                    max_whole_transcript <- max(end(subset_current_transcripts[[1]]),end(subset_current_transcripts[[j]]))
                    if(strand == "+" && max_whole_transcript != max(end(exon1),end(exons2))){
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                    if(strand == "-" && max_whole_transcript != max(end(exon1),end(exons2))){
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                  }
                }
              }
              if(length(exons2) > 1){
                ## should not be the case..
                print(paste(current_gene,"; i =",i,"; j=",j,"; l=",l,sep=""))
              } 
            }
          }
          
          ## Intron Retention
          ## less exons with same start /end
          if(sum(range_overlap) == n){
            intrRetention <- intrRetention + 1
            writeLines(paste(current_gene,"Intron Retention",sep="\t"),ffile)
            
            ## check for additional A3'E / A5'E
            ## check for common start and end just for overlapping exons, exclude part of retained intron
            ov <- findOverlaps(shorter_subset_current_transcripts,longer_subset_current_transcripts)
            n_short <- length(shorter_subset_current_transcripts)
            for(l in 1:n_short){
              exon1 <- shorter_subset_current_transcripts[l]
              hits <- ov[queryHits(ov) == l,]
              exons2 <- longer_subset_current_transcripts[subjectHits(hits)]
              
              if(length(exons2) == 1){
                ## just an overlap with one exon, test for A3'E / A5'E normally
                if(start(exon1) != start(exons2) && end(exon1) != end(exons2)){
                  ## we know already that they are overlapping, so we have both A3'E and A5'E
                  alt3End <- alt3End + 1
                  alt5End <- alt5End + 1
                  writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                  writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                }
                else{
                  if(start(exon1) != start(exons2)){
                    min_whole_transcript <- min(start(subset_current_transcripts[[1]]),start(subset_current_transcripts[[j]]))
                    if(strand == "+" && min_whole_transcript != min(start(exon1),start(exons2))){
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                    if(strand == "-" && min_whole_transcript != min(start(exon1),start(exons2))){
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                  }
                  if(end(exon1) != end(exons2)){
                    max_whole_transcript <- max(end(subset_current_transcripts[[1]]),end(subset_current_transcripts[[j]]))
                    if(strand == "+" && max_whole_transcript != max(end(exon1),end(exons2))){
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                    if(strand == "-" && max_whole_transcript != max(end(exon1),end(exons2))){
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                  }
                }      
              }
              if(length(exons2) > 1){
                ## more than one exon overlap, so we compare just the first start and last end
                if(start(exon1) != min(start(exons2)) && end(exon1) != max(end(exons2))){
                  ## we know already that they are overlapping, so we have both A3'E and A5'E
                  alt3End <- alt3End + 1
                  alt5End <- alt5End + 1
                  writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                  writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                }
                else{
                  if(start(exon1) != min(start(exons2))){
                    min_whole_transcript <- min(start(subset_current_transcripts[[1]]),start(subset_current_transcripts[[j]]))
                    if(strand == "+" && min_whole_transcript != min(start(exon1),min(start(exons2)))){
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                    if(strand == "-" && min_whole_transcript != min(start(exon1),min(start(exons2)))){
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                  }
                  if(end(exon1) != max(end(exons2))){
                    max_whole_transcript <- max(end(subset_current_transcripts[[1]]),end(subset_current_transcripts[[j]]))
                    if(strand == "+" && max_whole_transcript != max(end(exon1),max(end(exons2)))){
                      alt5End <- alt5End + 1
                      writeLines(paste(current_gene,"Alternative 5'End",sep="\t"),ffile)
                    }
                    if(strand == "-" && max_whole_transcript != max(end(exon1),max(end(exons2)))){
                      alt3End <- alt3End + 1
                      writeLines(paste(current_gene,"Alternative 3'End",sep="\t"),ffile)
                    }
                  }
                }
              }
            }
            
            if(sum(range_overlap) < n-1){
              ## could be mxe..
              print(paste(current_gene,"; i=",i,"; j=",j,sep=""))
            } 
          }
        }
      } 
    } 
  }
  writeLines("------------------------------------------------------",ffile)
  writeLines("AS summary:",ffile)
  writeLines(paste("Alternative Start Site:",altStartSite,sep="\t"),ffile)
  writeLines(paste("Alternative Termination Site:",altTerminationSite,sep="\t"),ffile)
  writeLines(paste("Alternative 5'End:",alt5End,sep="\t"),ffile)
  writeLines(paste("Alternative 3'End:",alt3End,sep="\t"),ffile)
  writeLines(paste("Mutually Exclusive Exon:",mxe,sep="\t"),ffile)
  writeLines(paste("Exon Skipping:",exonSkipping,sep="\t"),ffile)
  writeLines(paste("Intron Retention:",intrRetention,sep="\t"),ffile)
}

close(ffile)
