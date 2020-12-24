library(tidyverse)
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(clusterProfiler)

bed_ls2_vs_ls3_smad2enriched <- read_tsv("./bamfiles/diff_lap-smad2_vs_lap-smad3_c3.0_cond1.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ls3_smad2enriched <- GRanges(seqnames = bed_ls2_vs_ls3_smad2enriched$chrom, ranges = IRanges(start = bed_ls2_vs_ls3_smad2enriched$start, end = bed_ls2_vs_ls3_smad2enriched$end), log10_LR = bed_ls2_vs_ls3_smad2enriched$log10_likelihood_ratio)

ranges(ls2_vs_ls3_smad2enriched)

bed_ls2_vs_ls3_smad3enriched <- read_tsv("./bamfiles/diff_lap-smad2_vs_lap-smad3_c3.0_cond2.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ls3_smad3enriched <- GRanges(seqnames = bed_ls2_vs_ls3_smad3enriched$chrom, ranges = IRanges(start = bed_ls2_vs_ls3_smad3enriched$start, end = bed_ls2_vs_ls3_smad3enriched$end), log10_LR = bed_ls2_vs_ls3_smad3enriched$log10_likelihood_ratio)

bed_ls2_vs_ls3_bothenriched <- read_tsv("./bamfiles/diff_lap-smad2_vs_lap-smad3_c3.0_common.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ls3_bothenriched <- GRanges(seqnames = bed_ls2_vs_ls3_bothenriched$chrom, ranges = IRanges(start = bed_ls2_vs_ls3_bothenriched$start, end = bed_ls2_vs_ls3_bothenriched$end), log10_LR = bed_ls2_vs_ls3_bothenriched$log10_likelihood_ratio)

TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene ### loading tx database for annotation

ap_ls2_vs_ls3_smad2enriched <- annotatePeak(peak = ls2_vs_ls3_smad2enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_ls2_vs_ls3_smad3enriched <- annotatePeak(peak = ls2_vs_ls3_smad3enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_ls2_vs_ls3_bothenriched <- annotatePeak(peak = ls2_vs_ls3_bothenriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

### Creating tibbles from peak annotation data

df_ls2_vs_ls3_smad2enriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ls3_smad2enriched) %>% as_tibble()

df_ls2_vs_ls3_smad3enriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ls3_smad3enriched) %>% as_tibble()

df_ls2_vs_ls3_bothenriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ls3_bothenriched) %>% as_tibble()

save.image(file = "differential_binding.RData") 

df_ls2_vs_ls3_smad2enriched <- df_ls2_vs_ls3_smad2enriched %>% dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())

df_ls2_vs_ls3_smad3enriched <- df_ls2_vs_ls3_smad3enriched %>% dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())

df_ls2_vs_ls3_bothenriched <- df_ls2_vs_ls3_bothenriched %>% dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())

### Filtering only to those peaks located near gene promoters and deriving a unique (i.e., devoid of duplicates) vector of genes associated with these peaks for each condition

prom_genes_ls2_vs_ls3_smad2enriched <- df_ls2_vs_ls3_smad2enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()


prom_genes_ls2_vs_ls3_smad3enriched <-  df_ls2_vs_ls3_smad3enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()


prom_genes_ls2_vs_ls3_bothenriched <-  df_ls2_vs_ls3_bothenriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()



# lap-smad2-vs-lap-smad3-setops-on-gene-vectors --------------------------------------------------


intersect(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_smad3enriched$GENENAME)

intersect(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_bothenriched$GENENAME)

intersect(prom_genes_ls2_vs_ls3_smad3enriched$GENENAME, prom_genes_ls2_vs_ls3_bothenriched$GENENAME)

setdiff(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_smad3enriched$GENENAME)

setdiff(prom_genes_ls2_vs_ls3_smad3enriched$GENENAME, prom_genes_ls2_vs_ls3_smad2enriched$GENENAME)



# lap-smad2-vs-native-smad3-creating-GRanges -----------------------------------------------


bed_ls2_vs_ns3_smad2enriched <- read_tsv("./bamfiles/diff_lap-smad2_vs_native-smad3_c3.0_cond1.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ns3_smad2enriched <- GRanges(seqnames = bed_ls2_vs_ns3_smad2enriched$chrom, ranges = IRanges(start = bed_ls2_vs_ns3_smad2enriched$start, end = bed_ls2_vs_ns3_smad2enriched$end), log10_LR = bed_ls2_vs_ns3_smad2enriched$log10_likelihood_ratio)

ranges(ls2_vs_ns3_smad2enriched)

bed_ls2_vs_ns3_smad3enriched <- read_tsv("./bamfiles/diff_lap-smad2_vs_native-smad3_c3.0_cond2.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ns3_smad3enriched <- GRanges(seqnames = bed_ls2_vs_ns3_smad3enriched$chrom, ranges = IRanges(start = bed_ls2_vs_ns3_smad3enriched$start, end = bed_ls2_vs_ns3_smad3enriched$end), log10_LR = bed_ls2_vs_ns3_smad3enriched$log10_likelihood_ratio)

bed_ls2_vs_ns3_bothenriched <- read_tsv("./bamfiles/diff_lap-smad2_vs_native-smad3_c3.0_common.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ns3_bothenriched <- GRanges(seqnames = bed_ls2_vs_ns3_bothenriched$chrom, ranges = IRanges(start = bed_ls2_vs_ns3_bothenriched$start, end = bed_ls2_vs_ns3_bothenriched$end), log10_LR = bed_ls2_vs_ns3_bothenriched$log10_likelihood_ratio)


# lap-smad2-vs-native-smad3-annotating-peaks ------------------------------

TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene ### loading tx database for annotation


ap_ls2_vs_ns3_smad2enriched <- annotatePeak(peak = ls2_vs_ns3_smad2enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_ls2_vs_ns3_smad3enriched <- annotatePeak(peak = ls2_vs_ns3_smad3enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_ls2_vs_ns3_bothenriched <- annotatePeak(peak = ls2_vs_ns3_bothenriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

df_ls2_vs_ns3_smad2enriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ns3_smad2enriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())

df_ls2_vs_ns3_smad3enriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ns3_smad3enriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())

df_ls2_vs_ns3_bothenriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ns3_bothenriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())


# lap-smad2-vs-native-smad3-extracting-promoter-gene-lists ----------------


prom_genes_ls2_vs_ns3_smad2enriched <- df_ls2_vs_ns3_smad2enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_ls2_vs_ns3_smad3enriched <-  df_ls2_vs_ns3_smad3enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_ls2_vs_ns3_bothenriched <-  df_ls2_vs_ns3_bothenriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()



# lap-smad2-vs-native-smad3-setops-on-gene-vectors --------------------------------------

intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_smad3enriched$GENENAME)

intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME)

intersect(prom_genes_ls2_vs_ns3_smad3enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME)

save.image("differential_binding.RData")


# some-exploration --------------------------------------------------------

df_ls2_vs_ls3_smad2enriched %>% filter(GENENAME %in% intersect(prom_genes_ls2_vs_ls3_smad2enriched, prom_genes_ls2_vs_ls3_smad3enriched)) %>% dplyr::select(GENENAME, log10_LR)

### Removing common gene annotations occurring across smad2-enriched, smad3-enriched, and both-enriched, from the smad2-enriched and smad3-enriched gene lists

prom_genes_ls2_vs_ls3_bothenriched
prom_genes_ls2_vs_ls3_smad2enriched
prom_genes_ls2_vs_ls3_smad3enriched

intersect(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_bothenriched$GENENAME)

intersect(prom_genes_ls2_vs_ls3_smad3enriched$GENENAME, prom_genes_ls2_vs_ls3_bothenriched$GENENAME)

intersect(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_smad3enriched$GENENAME)

union(intersect(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_bothenriched$GENENAME), intersect(prom_genes_ls2_vs_ls3_smad3enriched$GENENAME, prom_genes_ls2_vs_ls3_bothenriched$GENENAME)) %>% 
  union(., intersect(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_smad3enriched$GENENAME))


ambi_ls2_vs_ls3 <- union(intersect(prom_genes_ls2_vs_ls3_smad2enriched$geneId, prom_genes_ls2_vs_ls3_bothenriched$geneId), intersect(prom_genes_ls2_vs_ls3_smad3enriched$geneId, prom_genes_ls2_vs_ls3_bothenriched$geneId)) %>% 
  union(., intersect(prom_genes_ls2_vs_ls3_smad2enriched$geneId, prom_genes_ls2_vs_ls3_smad3enriched$geneId))

prom_genes_ls2_vs_ls3_smad2enriched %>% filter(!geneId %in% ambi_ls2_vs_ls3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)

prom_genes_ls2_vs_ls3_smad3enriched %>% filter(!geneId %in% ambi_ls2_vs_ls3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)

ego_prom_genes_ls2_vs_ls3_smad2enriched <- enrichGO(gene = (prom_genes_ls2_vs_ls3_smad2enriched %>% filter(!geneId %in% ambi_ls2_vs_ls3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)

ego_prom_genes_ls2_vs_ls3_smad3enriched <- enrichGO(gene = (prom_genes_ls2_vs_ls3_smad3enriched %>% filter(!geneId %in% ambi_ls2_vs_ls3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)

####

prom_genes_ls2_vs_ns3_bothenriched
prom_genes_ls2_vs_ns3_smad2enriched
prom_genes_ls2_vs_ns3_smad3enriched

intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME)

intersect(prom_genes_ls2_vs_ns3_smad3enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME)

intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_smad3enriched$GENENAME)

union(intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME), intersect(prom_genes_ls2_vs_ns3_smad3enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME)) %>% 
  union(., intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_smad3enriched$GENENAME))


ambi_ls2_vs_ns3 <- union(intersect(prom_genes_ls2_vs_ns3_smad2enriched$geneId, prom_genes_ls2_vs_ns3_bothenriched$geneId), intersect(prom_genes_ls2_vs_ns3_smad3enriched$geneId, prom_genes_ls2_vs_ns3_bothenriched$geneId)) %>% 
  union(., intersect(prom_genes_ls2_vs_ns3_smad2enriched$geneId, prom_genes_ls2_vs_ns3_smad3enriched$geneId))

prom_genes_ls2_vs_ns3_smad2enriched %>% filter(!geneId %in% ambi_ls2_vs_ns3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)

prom_genes_ls2_vs_ns3_smad3enriched %>% filter(!geneId %in% ambi_ls2_vs_ns3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)

ego_prom_genes_ls2_vs_ns3_smad2enriched <- enrichGO(gene = (prom_genes_ls2_vs_ns3_smad2enriched %>% filter(!geneId %in% ambi_ls2_vs_ns3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)

ego_prom_genes_ls2_vs_ns3_smad3enriched <- enrichGO(gene = (prom_genes_ls2_vs_ns3_smad3enriched %>% filter(!geneId %in% ambi_ls2_vs_ns3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)


save.image("differential_binding.RData")

###

ego_prom_genes_ls2_vs_ls3_smad2enriched@result %>% 
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID) %>% View()

ego_prom_genes_ls2_vs_ls3_smad3enriched@result %>%
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID) %>% View()


####

ego_prom_genes_ls2_vs_ns3_smad2enriched@result %>% 
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID) %>% View()

ego_prom_genes_ls2_vs_ns3_smad3enriched@result %>%
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID) %>% View()




####

ego_prom_genes_ls2_vs_ls3_smad3enriched@result %>%
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(geneID) %>% unlist(use.names = FALSE) %>% head(1) %>% str_split(string = ., pattern = "\\/", simplify = FALSE) %>% unlist(use.names = FALSE) %>% select(x = org.Hs.eg.db, keys = ., keytype = "SYMBOL", columns = c("GENENAME"))