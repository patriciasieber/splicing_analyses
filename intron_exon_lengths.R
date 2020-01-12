setwd("/data/species_comparison/")

## intron lengths of all fungi
introns_all <- NULL ## TODO: fill!
names(introns_all) <- c("Sc","Ca","Cp","Cg","Af","Hc","Cn","Lc")

introns_all_rev <- rev(introns_all)

## exon lengths of all fungi

exons_all <- NULL ## TODO: fill!
names(exons_all) <- names(introns_all)

exons_all_rev <- rev(exons_all)

## read intron / exon lengths
introns_all_rev <- lapply(introns_all_rev,function(x){
  current_lengths <- readLines(x)
  current_lengths <- as.integer(current_lengths)
  return(current_lengths)
})

exons_all_rev <- lapply(exons_all_rev,function(x){
  current_lengths <- readLines(x)
  current_lengths <- as.integer(current_lengths)
  return(current_lengths)
})

svg("intron_lengths_boxplot.svg")
boxplot(introns_all_rev,horizontal=T,outline=F,
        xlab="Intron lengths",las=2)
dev.off()

svg("exon_lengths_boxplot.svg")
boxplot(exons_all_rev,horizontal=T,outline=F,
        xlab="Exon lengths",las=2)
dev.off()


summary_introns <- lapply(introns_all_rev,function(x){summary(x)})
