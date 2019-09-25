library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)



# function_for_extracting_ranges_for_ambiguous_genes ----------------------


df_ls2_vs_ns3_smad2enriched %>% filter(geneId %in% ambi_ls2_vs_ns3, str_detect(string = annotation, pattern = "Promoter")) %>% dplyr::select(geneId) %>% distinct(geneId) %>% unlist(use.names = FALSE) %>% is.element(el = ., set = ambi_ls2_vs_ns3) ### sanity check for filtering strategy

range_extract <- function(ambi_list, df){ # ambi_list = ambiguous gene vector, df = annotation dataframe
  ## extracts subset of annotation dataframe corresponding to ambiguous genes
  df_ambi_subset <- df %>% 
    filter(geneId %in% ambi_list, str_detect(string = annotation, pattern = "Promoter"))
  
  # creates ranges based on data from the ambiguous gene data frame subset
  
  ranges <- GRanges(seqnames = df_ambi_subset$seqnames, ranges = IRanges(start = df_ambi_subset$start, df_ambi_subset$end), strand = df_ambi_subset$strand, GENENAME = df_ambi_subset$GENENAME, SYMBOL = df_ambi_subset$SYMBOL, log10_LR = df_ambi_subset$log10_LR,geneId = df_ambi_subset$geneId, distanceToTSS = df_ambi_subset$distanceToTSS)
  
  return(ranges)
}



# creating_granges_list_from_dataframes ------------------------------------

## first set of granges is from Smad2 peaks, second set is from Smad3 peaks

ambi_ranges_ls2_vs_ns3 <- GRangesList(range_extract(ambi_list = ambi_ls2_vs_ns3, df = df_ls2_vs_ns3_smad2enriched), range_extract(ambi_list = ambi_ls2_vs_ns3, df = df_ls2_vs_ns3_smad3enriched))



mcols(ambi_ranges_ls2_vs_ns3) <- c("Smad2", "Smad3")

# naming granges list elements for easier access and identification

names(ambi_ranges_ls2_vs_ns3) <- c("Smad2", "Smad3")


ambi_ranges_ls2_vs_ns3$Smad2

ambi_ranges_ls2_vs_ns3$Smad3

### providing genome info for all the ranges in both the grange list elements


genome(ambi_ranges_ls2_vs_ns3) <- "hg38"


# processing_ranges -------------------------------------------------------


# creating annotation tracks for Gviz from individual grange list elements

smad2_track <- AnnotationTrack(range = ambi_ranges_ls2_vs_ns3$Smad2, name="Smad2")

smad3_track <- AnnotationTrack(range = ambi_ranges_ls2_vs_ns3$Smad3, name="Smad3")

# creating genome axis track

gtrack <- GenomeAxisTrack()

itrack <- IdeogramTrack(chromosome = "chr21", genome = "hg38")


# finding overlaps (note: none were found)

findOverlaps(query = ambi_ranges_ls2_vs_ns3$Smad2, ambi_ranges_ls2_vs_ns3$Smad3)


#### plotting ranges

seqnames(ambi_ranges_ls2_vs_ns3) # common chromosomes are 6, 9, 10 and 21


plotTracks(trackList = c(gtrack, smad2_track, smad3_track, IdeogramTrack(chromosome = "chr6", genome = "hg38")), chromosome = "chr6")

ambi_ranges_ls2_vs_ns3$Smad2[seqnames(ambi_ranges_ls2_vs_ns3$Smad2) == "chr6"] %>% 
  mcols() %>% 
  as_tibble() %>%
  mutate(TF = "Smad2") %>% 
full_join(x = ., y = (ambi_ranges_ls2_vs_ns3$Smad3[seqnames(ambi_ranges_ls2_vs_ns3$Smad3) == "chr6"] %>% mcols() %>% as_tibble()) %>% mutate(TF = "Smad3"))



plotTracks(trackList = c(gtrack, smad2_track, smad3_track, IdeogramTrack(chromosome = "chr9", genome = "hg38")), chromosome = "chr9")

ambi_ranges_ls2_vs_ns3$Smad2[seqnames(ambi_ranges_ls2_vs_ns3$Smad2) == "chr9"] %>% 
  mcols() %>% 
  as_tibble() %>%
  mutate(TF = "Smad2") %>% 
  full_join(x = ., y = (ambi_ranges_ls2_vs_ns3$Smad3[seqnames(ambi_ranges_ls2_vs_ns3$Smad3) == "chr9"] %>% mcols() %>% as_tibble()) %>% mutate(TF = "Smad3"))



plotTracks(trackList = c(gtrack, smad2_track, smad3_track, IdeogramTrack(chromosome = "chr10", genome = "hg38")), chromosome = "chr10")

ambi_ranges_ls2_vs_ns3$Smad2[seqnames(ambi_ranges_ls2_vs_ns3$Smad2) == "chr10"] %>% 
  mcols() %>% 
  as_tibble() %>%
  mutate(TF = "Smad2") %>% 
  full_join(x = ., y = (ambi_ranges_ls2_vs_ns3$Smad3[seqnames(ambi_ranges_ls2_vs_ns3$Smad3) == "chr10"] %>% mcols() %>% as_tibble()) %>% mutate(TF = "Smad3"))



plotTracks(trackList = c(gtrack, smad2_track, smad3_track, IdeogramTrack(chromosome = "chr21", genome = "hg38")), chromosome = "chr21")


ambi_ranges_ls2_vs_ns3$Smad2[seqnames(ambi_ranges_ls2_vs_ns3$Smad2) == "chr21"] %>% 
  mcols() %>% 
  as_tibble() %>%
  mutate(TF = "Smad2") %>% 
  full_join(x = ., y = (ambi_ranges_ls2_vs_ns3$Smad3[seqnames(ambi_ranges_ls2_vs_ns3$Smad3) == "chr21"] %>% mcols() %>% as_tibble()) %>% mutate(TF = "Smad3"))

plotTracks(trackList = c(gtrack, smad2_track, smad3_track, IdeogramTrack(chromosome = "chr21", genome = "hg38")), chromosome = "chr21", from = 8.19999*10^6, to = 8.261*10^6)



plotTracks(trackList = c(gtrack, smad2_track, smad3_track, IdeogramTrack(chromosome = "chr21", genome = "hg38")), chromosome = "chr21", from = 8.42*10^6, to = 8.455*10^6) 




### exporting range plots

pdf("chr21_8.44mb_to_8.45mb.pdf") 

plotTracks(trackList = c(gtrack, smad2_track, smad3_track), chromosome = "chr21", from = 8.44*10^6, to = 8.45*10^6) 

dev.off()

### plotting ranges 2

pdf("chr21_8.44mb_to_8.45mb.pdf") 

plotTracks(trackList = c(gtrack, smad2_track, smad3_track), chromosome = "21") 

dev.off() 



plotTracks(trackList = c((ambi_ranges_ls2_vs_ns3$Smad2[mcols(ambi_ranges_ls2_vs_ns3$Smad2)$GENENAME == "RNA, 5.8S ribosomal N5"] %>% AnnotationTrack(range = ., name = "Smad2")), (ambi_ranges_ls2_vs_ns3$Smad3[mcols(ambi_ranges_ls2_vs_ns3$Smad3)$GENENAME == "RNA, 5.8S ribosomal N5"] %>% AnnotationTrack(range = ., name = "Smad3")), gtrack), chromosome = "chr21")

export.bed(ambi_ranges_ls2_vs_ns3$Smad2, con = "ambi_ranges_ls2_vs_ns3_Smad2.bed")

export.bed(ambi_ranges_ls2_vs_ns3$Smad3, con = "ambi_ranges_ls2_vs_ns3_Smad3.bed")



#####

ambi_ranges_ls2_vs_ls3 <- GRangesList(range_extract(ambi_list = ambi_ls2_vs_ls3, df = df_ls2_vs_ls3_smad2enriched), range_extract(ambi_list = ambi_ls2_vs_ls3, df = df_ls2_vs_ls3_smad3enriched))


mcols(ambi_ranges_ls2_vs_ls3) <- c("Smad2", "Smad3")

names(ambi_ranges_ls2_vs_ls3) <- c("Smad2", "Smad3")



####


mcols(ambi_ranges_ls2_vs_ns3$Smad2)[mcols(ambi_ranges_ls2_vs_ns3$Smad2)$GENENAME == "RNA, 5.8S ribosomal N5", "distanceToTSS"]


mcols(ambi_ranges_ls2_vs_ns3$Smad3)[mcols(ambi_ranges_ls2_vs_ns3$Smad3)$GENENAME == "RNA, 5.8S ribosomal N5", "distanceToTSS"]


# getting-bed-files-for-motif-analysis -------------------------------------------------


ambi_ranges_ls2_vs_ns3$Smad2 %>% as.data.frame() %>% as_tibble()




