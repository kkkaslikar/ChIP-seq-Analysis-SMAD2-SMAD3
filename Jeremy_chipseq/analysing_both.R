
library(tidyverse)
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(clusterProfiler)
library(Gviz)

bed_tr_ls2_vs_ns3_smad2enriched <- read_tsv("./treated/treated-diffpeak_results/diff_lap-smad2_vs_lap-smad3_c3.0_cond1.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

tr_ls2_vs_ns3_smad2enriched <- GRanges(seqnames = bed_tr_ls2_vs_ns3_smad2enriched$chrom, ranges = IRanges(start = bed_tr_ls2_vs_ns3_smad2enriched$start, end = bed_tr_ls2_vs_ns3_smad2enriched$end), log10_LR = bed_tr_ls2_vs_ns3_smad2enriched$log10_likelihood_ratio)


bed_tr_ls2_vs_ns3_smad3enriched <- read_tsv("./treated/treated-diffpeak_results/diff_lap-smad2_vs_lap-smad3_c3.0_cond2.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

tr_ls2_vs_ns3_smad3enriched <- GRanges(seqnames = bed_tr_ls2_vs_ns3_smad3enriched$chrom, ranges = IRanges(start = bed_tr_ls2_vs_ns3_smad3enriched$start, end = bed_tr_ls2_vs_ns3_smad3enriched$end), log10_LR = bed_tr_ls2_vs_ns3_smad3enriched$log10_likelihood_ratio)

bed_tr_ls2_vs_ns3_bothenriched <- read_tsv("./treated/treated-diffpeak_results/diff_lap-smad2_vs_lap-smad3_c3.0_common.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

tr_ls2_vs_ns3_bothenriched <- GRanges(seqnames = bed_tr_ls2_vs_ns3_bothenriched$chrom, ranges = IRanges(start = bed_tr_ls2_vs_ns3_bothenriched$start, end = bed_tr_ls2_vs_ns3_bothenriched$end), log10_LR = bed_tr_ls2_vs_ns3_bothenriched$log10_likelihood_ratio)

bed_tr_ls2 <- read_tsv("./treated/treated-lap-smad2/lap-smad2_peaks.narrowPeak", col_names = c("chrom", "start", "end", "name"))

tr_ls2 <- GRanges(seqnames = bed_tr_ls2$chrom, ranges = IRanges(start = bed_tr_ls2$start, end = bed_tr_ls2$end))

bed_tr_ns3 <- read_tsv("./treated/treated-native-smad3/native-smad3_peaks.narrowPeak", col_names = c("chrom", "start", "end", "name"))

tr_ns3 <- GRanges(seqnames = bed_tr_ns3$chrom, ranges = IRanges(start = bed_tr_ns3$start, end = bed_tr_ns3$end))

TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene ### loading tx database for annotation



ap_tr_ls2_vs_ns3_smad2enriched <- annotatePeak(peak = tr_ls2_vs_ns3_smad2enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_tr_ls2_vs_ns3_smad3enriched <- annotatePeak(peak = tr_ls2_vs_ns3_smad3enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_tr_ls2_vs_ns3_bothenriched <- annotatePeak(peak = tr_ls2_vs_ns3_bothenriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_tr_ls2 <- annotatePeak(peak = tr_ls2, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_tr_ns3 <- annotatePeak(peak = tr_ns3, TxDb = TxDb, annoDb = "org.Hs.eg.db")

df_tr_ls2_vs_ns3_smad2enriched <- AnnotationDbi::as.data.frame(ap_tr_ls2_vs_ns3_smad2enriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, SYMBOL, geneId, ENSEMBL, everything())

df_tr_ls2_vs_ns3_smad3enriched <- AnnotationDbi::as.data.frame(ap_tr_ls2_vs_ns3_smad3enriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, SYMBOL, geneId, ENSEMBL, everything())

df_tr_ls2_vs_ns3_bothenriched <- AnnotationDbi::as.data.frame(ap_tr_ls2_vs_ns3_bothenriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, SYMBOL, geneId, ENSEMBL, everything())

df_tr_ls2 <- AnnotationDbi::as.data.frame(ap_tr_ls2) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, SYMBOL, geneId, ENSEMBL, everything())

df_tr_ns3 <- AnnotationDbi::as.data.frame(ap_tr_ns3) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, SYMBOL, geneId, ENSEMBL, everything())

GenomicRanges::findOverlaps(query = tr_ls2, subject = tr_ls2_vs_ns3_bothenriched)

tr_ls2[46]

tr_ls2_vs_ns3_bothenriched[2]

hits_1 <- GenomicRanges::findOverlaps(query = tr_ns3, subject = tr_ls2_vs_ns3_smad3enriched)

tr_ns3[hits_1@from]


####

track_tr_ls2_vs_ns3_smad2enriched <- AnnotationTrack(range = tr_ls2_vs_ns3_smad2enriched, name="ls2_vs_ns3_smad2enriched")

track_tr_ls2_vs_ns3_smad3enriched <- AnnotationTrack(range = tr_ls2_vs_ns3_smad3enriched, name="ls2_vs_ns3_smad3enriched")

track_tr_ls2_vs_ns3_bothenriched <- AnnotationTrack(range = tr_ls2_vs_ns3_bothenriched, name="ls2_vs_ns3_bothenriched")

track_tr_ls2 <- AnnotationTrack(range = tr_ls2, name= "lap-smad2")

track_tr_ns3 <- AnnotationTrack(range = tr_ns3, name="native-smad3")


gtrack <- GenomeAxisTrack()



plotTracks(trackList = c(gtrack, track_tr_ls2_vs_ns3_smad2enriched, track_tr_ls2_vs_ns3_smad3enriched), chromosome = "chr1")

#####

prom_genes_tr_ls2_vs_ns3_smad2enriched <- df_tr_ls2_vs_ns3_smad2enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_tr_ls2 <- df_tr_ls2 %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_tr_ls2_vs_ns3_smad3enriched <- df_tr_ls2_vs_ns3_smad3enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_tr_ns3 <- df_tr_ns3%>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_tr_ls2_vs_ns3_bothenriched <- df_tr_ls2_vs_ns3_bothenriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

#####

intersect((prom_genes_tr_ls2 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE)), (prom_genes_tr_ls2_vs_ns3_smad2enriched%>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE)))

intersect((prom_genes_tr_ls2 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE)), (prom_genes_tr_ls2_vs_ns3_bothenriched %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE)))


prom_genes_tr_ls2

prom_genes_tr_ls2_vs_ns3_smad2enriched

prom_genes_tr_ls2_vs_ns3_bothenriched

prom_genes_tr_ns3

prom_genes_tr_ls2_vs_ns3_smad3enriched


####

ego_prom_genes_tr_ns3 <- enrichGO(gene = (prom_genes_tr_ns3 %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)


ego_prom_genes_tr_ns3@result %>% 
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID) %>% View()

#####

x <- intersect((prom_genes_tr_ls2 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE)), (prom_genes_tr_ns3 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE)))

y <- intersect((prom_genes_tr_ls2 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE)), (prom_genes_tr_ls2_vs_ns3_smad3enriched %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE)))

intersect(x, y)

setdiff(x, y)

setdiff(y, x)





prom_genes_tr_ns3 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE) %>% length()

prom_genes_tr_ls2_vs_ns3_smad3enriched %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE) %>% length()

####

setdiff(prom_genes_tr_ns3 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE), prom_genes_tr_ls2_vs_ns3_smad3enriched %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE))

setdiff( prom_genes_tr_ls2_vs_ns3_smad3enriched %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE), prom_genes_tr_ns3  %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE))

z <- setdiff(prom_genes_tr_ns3 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE), prom_genes_tr_ls2_vs_ns3_smad3enriched %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE))

ego_prom_genes_tr_ns3 <- enrichGO(gene = prom_genes_tr_ns3 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)

ego_prom_genes_tr_ns3 <- enrichGO(gene = prom_genes_tr_ns3 %>% dplyr::select(GENENAME) %>% unlist(use.names = FALSE), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)

save.image(file = "analysing_both.RData")
