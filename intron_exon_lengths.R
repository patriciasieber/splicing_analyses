setwd("/data/PhD/species_comparison_revision/")

## intron lengths of all fungi
introns_afum <- "/data/PhD/species_comparison/A_fumigatus/stringtie_merge_g1.gtf_intronLength"
introns_calb <- "/data/PhD/species_comparison/C_albicans/stringtie_merge_g1.gtf_intronLength"
introns_calbk <- "/data/PhD/species_comparison/C_albicans_kaemmer/stringtie2_merge_g1.gtf_intronLength"
introns_cgla <- "/data/PhD/species_comparison/C_glabrata/stringtie_merge_g1.gtf_intronLength"
introns_cglak <- "/data/PhD/species_comparison/C_glabrata_kaemmer/stringtie2_merge_g1.gtf_intronLength"
introns_cpar <- "/data/PhD/species_comparison/C_parapsilosis/stringtie_merge_g1.gtf_intronLength"
introns_cpark <- "/data/PhD/species_comparison/C_parapsilosis_kaemmer/stringtie2_merge_g1.gtf_intronLength"
introns_cneo <- "/data/PhD/species_comparison/C_neoformans/stringtie_merge_g1.gtf_intronLength"
introns_lcor <- "/data/PhD/species_comparison/L_corymbifera/stringtie_merge_g1.gtf_intronLength"
introns_hcap <- "/data/PhD/species_comparison/H_capsulatum/stringtie_merge_g1.gtf_intronLength"

introns_scer <- "/data/PhD/species_comparison/S_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.31.gff3_intronLength"

introns_all <- list(introns_afum,introns_calb,introns_calbk,introns_cgla,introns_cglak,introns_cpar,introns_cpark,introns_cneo,introns_lcor,introns_hcap)
names(introns_all) <- c("Af","Ca","CaK","Cg","CgK","Cp","CpK","Cn","Lc","Hc")

introns_all <- list(introns_scer,introns_calbk,introns_cpark,introns_cglak,introns_afum,introns_hcap,introns_cneo,introns_lcor)
names(introns_all) <- c("Sc","Ca","Cp","Cg","Af","Hc","Cn","Lc")

introns_all_rev <- rev(introns_all)

## exon lengths of all fungi
exons_afum <- "/data/PhD/species_comparison/A_fumigatus/stringtie_merge_g1.gtf_exonLength"
exons_calb <- "/data/PhD/species_comparison/C_albicans/stringtie_merge_g1.gtf_exonLength"
exons_calbk <- "/data/PhD/species_comparison/C_albicans_kaemmer/stringtie2_merge_g1.gtf_exonLength"
exons_cgla <- "/data/PhD/species_comparison/C_glabrata/stringtie_merge_g1.gtf_exonLength"
exons_cglak <- "/data/PhD/species_comparison/C_glabrata_kaemmer/stringtie2_merge_g1.gtf_exonLength"
exons_cpar <- "/data/PhD/species_comparison/C_parapsilosis/stringtie_merge_g1.gtf_exonLength"
exons_cpark <- "/data/PhD/species_comparison/C_parapsilosis_kaemmer/stringtie2_merge_g1.gtf_exonLength"
exons_cneo <- "/data/PhD/species_comparison/C_neoformans/stringtie_merge_g1.gtf_exonLength"
exons_lcor <- "/data/PhD/species_comparison/L_corymbifera/stringtie_merge_g1.gtf_exonLength"
exons_hcap <- "/data/PhD/species_comparison/H_capsulatum/stringtie_merge_g1.gtf_exonLength"

exons_scer <- "/data/PhD/species_comparison/S_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.31.gff3_exonLength"


exons_all <- list(exons_afum,exons_calb,exons_calbk,exons_cgla,exons_cglak,exons_cpar,exons_cpark,exons_cneo,exons_lcor,exons_hcap)
names(exons_all) <- c("Af","Ca","CaK","Cg","CgK","Cp","CpK","Cn","Lc","Hc")

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
#png("intron_lengths_boxplot.png")
boxplot(introns_all_rev,horizontal=T,outline=F,
        xlab="Intron lengths",las=2)
dev.off()

#svg("exon_lengths_boxplot.svg")
png("exon_lengths_boxplot.png")
boxplot(exons_all_rev,horizontal=T,outline=F,
        xlab="Exon lengths",las=2)
dev.off()


summary_introns <- lapply(introns_all_rev,function(x){summary(x)})
