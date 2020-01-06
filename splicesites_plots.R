## barplot for difference between used splicesites
library(colorspace)
setwd("/data/PhD/species_comparison_revision/")
library(RColorBrewer)
library(pheatmap)


## splice sites heatmap for all fungi and conditions ####
afum <- c(99.18,0.70,0.09,NA,NA,NA,99.14,0.72,0.10)
#afum293 <- c(99.15,0.71,0.10,99.18,0.70,0.09,98.75,0.92,0.22)
calb <- c(88.35,2.87,5.66,85.58,3.14,7.26,91.38,2.38,4.16)
calbk <- c(88.71,2.85,5.47,87.80,2.87,5.99,91.44,2.37,4.12)
cgla <- c(87.37,4.04,5.34,NA,NA,NA,87.37,4.04,5.34)
cglak <- c(88.42,3.75,4.60,90.15,3.15,3.99,86.40,4.41,5.41)
cpar <- c(99.889,0,0.1,NA,NA,NA,99.94,0,0.06)
cpark <- c(99.94,0,0.06,99.78,0.06,0.16,99.94,0,0.06)
cneo <- c(97.69,2.05,0.26,NA,NA,NA,97.75,2.04,0.21)
lcor <- c(96.78,2.97,0.24,96.78,2.98,0.24,96.78,2.97,0.24)
hcap <- c(96.23,2.13,1.64,96.23,2.13,1.64,NA,NA,NA)

#all <- c(afum,calbk,cglak,cpark,cneo,lcor,hcap)

all <- c(afum,calb,calbk,cgla,cglak,cpar,cpark,cneo,lcor,hcap)
matrix_all <- matrix(all,nrow=3)

others <- apply(matrix_all,2,function(x){
  return(100-sum(x))
})
matrix_all <- rbind(matrix_all,others)

rownames(matrix_all) <- c("GT-AG", "GC-AG", "AT-AC","others")

names_allss <- c("Af","Ca","CaK","Cg","CgK","Cp","CpK","Cn","Lc","Hc")
#names_allss <- c("Af","Ca","Cg","Cp","Cn","Lc","Hc")

names_allss <- unlist(lapply(names_allss,function(x){return(rep(x,3))}))
conditions <- c("control","infection","stress")
conditions <- rep(conditions,10)
#conditions <- rep(conditions,7)

all <- paste(names_allss,conditions,sep="_")

colnames(matrix_all) <- all

log_matrix_all <- log10(matrix_all)
log_matrix_all <- unlist(lapply(log_matrix_all,function(x){
  if(!is.finite(x)) x <- NA
  return(x)
}))
log_matrix_all <- matrix(log_matrix_all,nrow=4)
rownames(log_matrix_all) <- c("GT-AG", "GC-AG", "AT-AC","others")
colnames(log_matrix_all) <- all
# 
#heatmap(log_matrix_all, Rowv=NA, Colv=NA, col = brewer.pal(9, "Blues"), scale="column", margins=c(5,10))

png("/data/PhD/species_comparison/splicesites_log.png",width=1600,height=300)
pheatmap(log_matrix_all,cellwidth=40,cellheight=40,width=1800,height=600,
         cluster_rows=F,cluster_cols=F,
         scale="none",display_numbers=T,margins=c(10,20))
dev.off()

png("/data/PhD/species_comparison/splicesites.png",width=1600,height=300)
pheatmap(matrix_all,cellwidth=40,cellheight=40,width=1800,height=600,
         cluster_rows=F,cluster_cols=F,
         scale="none",display_numbers=T,margins=c(10,20))
#,color=brewer.pal(9, "Blues")
dev.off()


matrix_all_boolean <- as.vector(apply(matrix_all,2,function(x){
  if(is.na(sum(x))){ return(FALSE)
  }else return(TRUE)
}))
matrix_all_woNA <- matrix_all[,matrix_all_boolean]


png("/data/PhD/species_comparison/splicesites_phylogeny.png",height=300,width=1200)
pheatmap(matrix_all_woNA,cellwidth=40,cellheight=40,width=1800,height=600,
         cluster_rows=F,cluster_cols=T,
         scale="none",display_numbers=T,margins=c(10,20))
#color=brewer.pal(9, "Blues")
dev.off()


png("/data/PhD/species_comparison/splicesites_phylogeny.png",height=300,width=1200)
pheatmap(matrix_all_woNA,cellwidth=40,cellheight=40,width=1800,height=600,
         cluster_rows=F,cluster_cols=T,
         scale="none",display_numbers=T,margins=c(10,20))
dev.off()

# log_matrix_all <- log10(matrix_all_woNA)
# cols <- colnames(log_matrix_all)
# log_matrix_all <- unlist(lapply(log_matrix_all,function(x){
#   if(!is.finite(x)) x <- NA
#   return(x)
# }))

# log_matrix_all <- matrix(log_matrix_all,nrow=4)
# rownames(log_matrix_all) <- c("GT-AG", "GC-AG", "AT-AC","others")
# colnames(log_matrix_all) <- cols
# 
# png("/data/PhD/species_comparison/splicesites_phylogeny_log.png",width=1200,height=300)
# pheatmap(log_matrix_all,cellwidth=30,cellheight=40,width=1800,height=600,
#          cluster_rows=F,cluster_cols=T,
#          scale="none",display_numbers=T,margins=c(10,20))
# dev.off()

##mean values for each species over all conditions
names_unique <- unique(names_allss)
mean_values <- lapply(names_unique,function(x){
  current <- matrix_all[1,][grep(x,names(matrix_all[1,]))]
  return(mean(current,na.rm=T))
})
 

####################################
afum <- c(99.18,0.70,0.09,NA,NA,NA,99.14,0.72,0.10)
#afum293 <- c(99.15,0.71,0.10,99.18,0.70,0.09,98.75,0.92,0.22)
calb <- c(88.35,2.87,5.66,85.58,3.14,7.26,91.38,2.38,4.16)
calbk <- c(88.71,2.85,5.47,87.80,2.87,5.99,91.44,2.37,4.12)
cgla <- c(87.37,4.04,5.34,NA,NA,NA,87.37,4.04,5.34)
cglak <- c(88.42,3.75,4.60,90.15,3.15,3.99,86.40,4.41,5.41)
cpar <- c(99.889,0,0.1,NA,NA,NA,99.94,0,0.06)
cpark <- c(99.94,0,0.06,99.78,0.06,0.16,99.94,0,0.06)
cneo <- c(97.69,2.05,0.26,NA,NA,NA,97.75,2.04,0.21)
lcor <- c(96.78,2.97,0.24,96.78,2.98,0.24,96.78,2.97,0.24)
hcap <- c(96.23,2.13,1.64,96.23,2.13,1.64,NA,NA,NA)

all <- c(calbk,cpark,cglak,afum,hcap,cneo,lcor)
matrix_all <- matrix(all,nrow=3)

others <- apply(matrix_all,2,function(x){
  return(100-sum(x))
})
matrix_all <- rbind(matrix_all,others)

rownames(matrix_all) <- c("GT-AG", "GC-AG", "AT-AC","others")

names_allss <- c("Ca","Cp","Cg","Af","Hc","Cn","Lc")
names_allss <- unlist(lapply(names_allss,function(x){return(rep(x,3))}))
conditions <- c("control","infection","stress")
conditions <- rep(conditions,7)
all <- paste(names_allss,conditions,sep="_")

colnames(matrix_all) <- all

#png("/data/PhD/species_comparison/splicesites_log.png",width=1600,height=300)
svg("/data/PhD/species_comparison/splicesites_log.svg",width=18,height=4)
pheatmap(log_matrix_all,cellwidth=40,cellheight=40,width=1800,height=600,
         cluster_rows=F,cluster_cols=F,
         scale="none",display_numbers=T,margins=c(10,20))
dev.off()

library(pheatmap)
#png("/data/PhD/species_comparison/splicesites.png",width=1600,height=300)
svg("/data/PhD/species_comparison/splicesites.svg",width=15,height=4)
pheatmap(matrix_all,cellwidth=40,cellheight=40,width=1800,height=600,
         cluster_rows=F,cluster_cols=F,
         scale="none",display_numbers=T,margins=c(10,20))
#,color=brewer.pal(9, "Blues")
dev.off()



matrix_all_boolean <- as.vector(apply(matrix_all,2,function(x){
  if(is.na(sum(x))){ return(FALSE)
  }else return(TRUE)
}))
matrix_all_woNA <- matrix_all[,matrix_all_boolean]


png("/data/PhD/species_comparison/splicesites_phylogeny.png",height=300,width=1200)
svg("/data/PhD/species_comparison/splicesites_phylogeny.svg",width=15,height=4)
pheatmap(matrix_all_woNA,cellwidth=40,cellheight=40,width=1800,height=600,
         cluster_rows=F,cluster_cols=T,
         scale="none",display_numbers=T,margins=c(10,20))
#color=brewer.pal(9, "Blues")
dev.off()



paletteLength <- 20
myColor <- colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(paletteLength)  ## RdYlBu,RdBu,Spectral
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(matrix_all,na.rm=T), 1, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(matrix_all,na.rm=T)/paletteLength, max(matrix_all,na.rm=T), length.out=floor(paletteLength/2)))

svg("splicesites.svg",width=15,height=4)
pdf("splicesites.pdf",width=15,height=4,onefile = F)
pheatmap(matrix_all,cellwidth=40,cellheight=40,width=1800,height=600,
         cluster_rows=F,cluster_cols=F,color=myColor, breaks=myBreaks,
         scale="none",display_numbers=T,margins=c(10,20))
dev.off()
