library(BSgenome)

exclude_no_change <- function(exonsAndIntrons,replicates_condition1,replicates_condition2,fasta,threshold_junc=10,threshold_ratio=2){
   transcripts <- unique(as.character(exonsAndIntrons$Parent))
   len_tr <- length(transcripts)
   as_genes <- NULL ## collect all genes that show alternative splicing
   if(len_tr > 0){
      for(i in 1:len_tr){
       #print(i)
       #flush.console()
       cat(i)
       #current_transcript <- exonsAndIntrons[substr(exonsAndIntrons$ID,1,nchar(exonsAndIntrons$ID)-3) == transcripts[i],]
       current_transcript <- exonsAndIntrons[as.character(exonsAndIntrons$Parent) == transcripts[i],]
       gRange <- ranges(current_transcript)
       geneRange <- GRanges(seqnames(current_transcript)[1],IRanges(start=min(start(gRange)),end=max(end(gRange))),strand(current_transcript)[1])
       len_rep <- length(replicates_condition1)
      
      ### part1: compare ratio of annotated and not-annotated junctions ####
      ## to identify ES, alt. start/end, MXE, alt 5'/3' end
       junc_hits1 <- NULL
       junc_hits2 <- NULL
       for(j in 1:len_rep){
         current_replicate1 <- get(replicates_condition1[j])
         ov <- findOverlaps(geneRange,current_replicate1)
         current_hits <- current_replicate1[subjectHits(ov),]
         junc_hits1 <- c(junc_hits1,current_hits)
         
         current_replicate2 <- get(replicates_condition2[j])
         ov <- findOverlaps(geneRange,current_replicate2)
         current_hits <- current_replicate2[subjectHits(ov),]
         junc_hits2 <- c(junc_hits2,current_hits)
       }
       
      if(as.character(strand(current_transcript)[1]) == "+"){ ## plus strang einträge zählen
        junc_hits1 <- unlist(lapply(junc_hits1,function(x){
          x$score <- x$plus_score
          return(x)
        }))
        junc_hits2 <- unlist(lapply(junc_hits2,function(x){
          x$score <- x$plus_score
          return(x)
        }))
      }else{ ## minus strang einträge zählen
        junc_hits1 <- unlist(lapply(junc_hits1,function(x){
          x$score <- x$minus_score
          return(x)
        }))
        junc_hits2 <- unlist(lapply(junc_hits2,function(x){
          x$score <- x$minus_score
          return(x)
        }))
      }
      ## sum same ranges
      junc_hits1 <- do.call(c,junc_hits1)
      junc_hits2 <- do.call(c,junc_hits2)
      
      junc_counts1 <- GRanges()
      junc_counts2 <- GRanges()
      while(length(junc_hits1) > 0){
        current_hit <- junc_hits1[1]
        same_range <- junc_hits1[ranges(junc_hits1) == ranges(current_hit),]
        count <- sum(same_range$score)
        current_count <- GRanges(seqnames(current_hit),ranges(current_hit),as.character(strand(current_transcript)[1]),score=count)
        junc_counts1 <- c(junc_counts1,current_count)
        junc_hits1 <- junc_hits1[ranges(junc_hits1) != ranges(current_hit),]
      }
      while(length(junc_hits2) > 0){
        current_hit <- junc_hits2[1]
        same_range <- junc_hits2[ranges(junc_hits2) == ranges(current_hit),]
        count <- sum(same_range$score)
        current_count <- GRanges(seqnames(current_hit),ranges(current_hit),as.character(strand(current_transcript)[1]),score=count)
        junc_counts2 <- c(junc_counts2,current_count)
        junc_hits2 <- junc_hits2[ranges(junc_hits2) != ranges(current_hit),]
      }
      
      junc_counts1 <- junc_counts1[junc_counts1$score > threshold_junc,]
      junc_counts2 <- junc_counts2[junc_counts2$score > threshold_junc,]

      junc_ranges <- c(ranges(junc_counts1),ranges(junc_counts2))
      junc_ranges <- unique(junc_ranges)
      
      ## do not compare if junctions have the same position like the annotated introns
      introns <- current_transcript[current_transcript$type == "intron",]
      if(length(introns) > 0){
        for(k in 1:length(introns)){
          junc_ranges <- junc_ranges[junc_ranges != ranges(introns[k])]
        }
      }
      
      # check intron borders: GT..AG, otherwise don't consider them (assumed to be artefacts)
      if(length(junc_ranges) > 0){
        chromosome <- as.character(seqnames(current_transcript)[1])
        n_chrom_name <- nchar(chromosome)
        chromosome_sequence <- fasta[substr(names(fasta),1,n_chrom_name) == chromosome,]
        #chromosome_sequence <- fasta[names(fasta) == chromosome,]
        
        strand <- strand(current_transcript)[1]
        
        valid_junc_ranges <- NULL
        for(k in 1:length(junc_ranges)){
          current_ranges <- junc_ranges[k]
          current_seq <- subseq(chromosome_sequence,start=start(current_ranges),end=end(current_ranges))
          seq_junc <- toString(current_seq)
          width_junc <- width(current_seq)
          if(as.character(strand == "+")){
            if(substr(seq_junc,1,2) == "GT" & substr(seq_junc,width_junc-1,width_junc) == "AG"){
              valid_junc_ranges <- c(valid_junc_ranges,junc_ranges[k])
            }
          }else{
            if(substr(seq_junc,1,2) == "CT" & substr(seq_junc,width_junc-1,width_junc) == "AC"){ ## reverse complement
              valid_junc_ranges <- c(valid_junc_ranges,junc_ranges[k])
            }
          }
        }
        junc_ranges <- valid_junc_ranges
      }
      
      
      ## compare ratios between those and annotated ones
      if(length(junc_ranges) > 0){
        logFC_intron <- introns$logFoldChange_junctions
        if(all(!is.na(logFC_intron))){
          for(k in 1:length(junc_ranges)){
            current_ranges <- junc_ranges[[k]]
            current_junc_counts1 <- junc_counts1[ranges(junc_counts1) == current_ranges,]
            current_junc_counts2 <- junc_counts2[ranges(junc_counts2) == current_ranges,]
            if(length(current_junc_counts1) > 0 && length(current_junc_counts2) > 0){
              fc_nointron <- foldchange(current_junc_counts1$score,current_junc_counts2$score)
              logFC_nointron <- foldchange2logratio(fc_nointron,base=2)
              
              ## compare logFCs of introns and nointrons
              fc_of_fc <- foldchange(logFC_intron,logFC_nointron)
              logFC_fc <- foldchange2logratio(fc_of_fc)
              
              logFC_fc <- logFC_fc[!is.na(logFC_fc)]
              ## threshold > 1, < -1
              if(length(logFC_fc) > 0){
                if(any(abs(logFC_fc) >= threshold_ratio)) as_genes <- c(as_genes,current_transcript$gene[1])
              }
            }else{ ## junctions do not appear in both conditions, could be AS
              as_genes <- c(as_genes,current_transcript$gene[1])}
          }
        }
      }
      
      ### part2: compare ratio of exons (pairwise) ####
      ## to identify ES, alt. start / end, MXE
      exons <- current_transcript[current_transcript$type == "exon",]
      exon_ratio <- foldchange(exons$rpkm1,exons$rpkm2)
      logFC_exons <- foldchange2logratio(exon_ratio, base=2)
      
      if(length(logFC_exons) > 1 && !any(is.na(logFC_exons))){
        ## logFC have different signs (+ and -)
        if(sum(sign(logFC_exons)) != length(logFC_exons)) as_genes <- c(as_genes,current_transcript$gene[1])
        ## high difference of exon expression, or no expression for one condition (Inf)
        if(sum(sign(logFC_exons)) == length(logFC_exons) && any(diff(logFC_exons) > threshold_ratio || is.na(diff(logFC_exons)))) as_genes <- c(as_genes,current_transcript$gene[1]) 
      }
      
      ### part3: compare ratio of introns and junctions (pairwise) ####
      ## to identify IR
      intron_ratio <- foldchange(introns$rpkm1,introns$rpkm2)
      logFC_introns <- foldchange2logratio(intron_ratio, base=2)
      
      ## logFC have different signs (+ and -)
      ## if NA, then there are no junctions, or missing expression -> no IR
      if(!(any(is.na(introns$logFoldChange_junctions))||(any(is.na(logFC_introns))))){  
        #print(paste0("i:",i,",transcript:",current_transcript$gene[1]))
        if(sum(sign(logFC_introns)) != sum(sign(introns$logFoldChange_junctions))){ as_genes <- c(as_genes,current_transcript$gene[1])
        }else{
          difff <- abs(logFC_introns - introns$logFoldChange_junctions)
          for(l in 1:length(difff)){
            if(!is.na(difff[l]) && any(difff[l] > threshold_ratio)) as_genes <- c(as_genes,current_transcript$gene[1])
          }  
        }
      } 
     }
   }
  as_genes <- unique(as_genes)
  
  exonsAndIntrons_new <- unlist(lapply(as_genes,function(x){
    return(exonsAndIntrons[exonsAndIntrons$gene == x,])
  }))
  if(!is.null(exonsAndIntrons_new)) exonsAndIntrons_new <- do.call(c,exonsAndIntrons_new)
  
  return(exonsAndIntrons_new)
}
