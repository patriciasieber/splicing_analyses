library(RColorBrewer)
library(pheatmap)

## expressed genes /isoforms for all fungi and conditions ####
afum <- c(8384,11320,NA,NA,8424,11236)
calbk <- c(5624,5942,5673,5972,5683,5975)
cglak <- c(4717,4760,4627,4648,4839,4891)
cpark <- c(5743,5784,5127,5144,5718,5751)
cneo <- c(5866,6597,NA,NA,6415,7219)
lcor <- c(8853,10408,8144,9277,8760,10322)
hcap <- c(8984,11570,8317,10439,NA,NA)

all <- c(calbk,cpark,cglak,afum,hcap,cneo,lcor)


matrix_all <- matrix(all,nrow=2)
rownames(matrix_all) <- c("genes", "transcripts")


names_allss <- c("Ca","Cp","Cg","Af","Hc","Cn","Lc")

names_allss <- unlist(lapply(names_allss,function(x){return(rep(x,3))}))
conditions <- c("control","infection","stress")

conditions <- rep(conditions,7)

all <- paste(names_allss,conditions,sep="_")

colnames(matrix_all) <- all


pdf("/data/species_comparison/expression.pdf")
pheatmap(matrix_all,cellwidth=40,cellheight=40,
         cluster_rows=F,cluster_cols=F,
         scale="none",display_numbers=T,width=200,height=200)
#,color=brewer.pal(9, "Blues")
dev.off()


## ("CaK","CpK","CgK","Af","Hc","Cn","Lc")
number_tr <- c(6719,5897,5108,13879,17182,9650,14015) 
number_genes <- c(6158,5845,4977,10144,11216,7543,11350)

normalize_genes <- unlist(lapply(number_genes,function(x){return(rep(x,3))}))
normalize_tr <- unlist(lapply(number_tr,function(x){return(rep(x,3))}))


matrix_all_normalized <- matrix_all
matrix_all_normalized[1,] <- matrix_all_normalized[1,]/normalize_genes
matrix_all_normalized[2,] <- matrix_all_normalized[2,]/normalize_tr


svg("/data/PhD/species_comparison_revision/expression_normalized.svg",width=15,height=3)
pheatmap(matrix_all_normalized,cellwidth=40,cellheight=40,
         cluster_rows=F,cluster_cols=F,
         scale="none",display_numbers=T,width=400,height=400)
#,color=brewer.pal(9, "Blues")
dev.off()


matrix_all_boolean <- as.vector(apply(matrix_all_normalized,2,function(x){
  if(is.na(sum(x))){ return(FALSE)
  }else return(TRUE)
}))
matrix_all_woNA <- matrix_all_normalized[,matrix_all_boolean]


pdf("/data/species_comparison/expression_normalized_clustered.pdf",height=3,width=14)
pheatmap(matrix_all_woNA,cellwidth=40,cellheight=40,
         cluster_rows=F,cluster_cols=T,
         scale="none",display_numbers=T,width=400,height=400)
#,color=brewer.pal(9, "Blues")
dev.off()
