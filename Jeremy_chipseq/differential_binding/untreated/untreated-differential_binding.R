library(tidyverse)
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(clusterProfiler)

file_list <- list.files(full.names = TRUE, recursive = TRUE)

file_list[str_detect(string = file_list, pattern = "cond.+.bed")]


### creating granges

bed_ls2_vs_ls3_smad2enriched <- read_tsv("./differential_binding/untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_lap-smad3_c3.0_cond1.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ls3_smad2enriched <- GRanges(seqnames = bed_ls2_vs_ls3_smad2enriched$chrom, ranges = IRanges(start = bed_ls2_vs_ls3_smad2enriched$start, end = bed_ls2_vs_ls3_smad2enriched$end), log10_LR = bed_ls2_vs_ls3_smad2enriched$log10_likelihood_ratio)


bed_ls2_vs_ls3_smad3enriched <- read_tsv("./differential_binding/untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_lap-smad3_c3.0_cond2.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ls3_smad3enriched <- GRanges(seqnames = bed_ls2_vs_ls3_smad3enriched$chrom, ranges = IRanges(start = bed_ls2_vs_ls3_smad3enriched$start, end = bed_ls2_vs_ls3_smad3enriched$end), log10_LR = bed_ls2_vs_ls3_smad3enriched$log10_likelihood_ratio)

bed_ls2_vs_ls3_bothenriched <- read_tsv("./differential_binding/untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_lap-smad3_c3.0_common.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ls3_bothenriched <- GRanges(seqnames = bed_ls2_vs_ls3_bothenriched$chrom, ranges = IRanges(start = bed_ls2_vs_ls3_bothenriched$start, end = bed_ls2_vs_ls3_bothenriched$end), log10_LR = bed_ls2_vs_ls3_bothenriched$log10_likelihood_ratio)

### annotating peaks

TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene ### loading tx database for annotation

ap_ls2_vs_ls3_smad2enriched <- annotatePeak(peak = ls2_vs_ls3_smad2enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_ls2_vs_ls3_smad3enriched <- annotatePeak(peak = ls2_vs_ls3_smad3enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_ls2_vs_ls3_bothenriched <- annotatePeak(peak = ls2_vs_ls3_bothenriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

### Creating tibbles from peak annotation data

df_ls2_vs_ls3_smad2enriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ls3_smad2enriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())

df_ls2_vs_ls3_smad3enriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ls3_smad3enriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())

df_ls2_vs_ls3_bothenriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ls3_bothenriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())

### Extracting unique list of genes containing promoter-proximal peaks

prom_genes_ls2_vs_ls3_smad2enriched <- df_ls2_vs_ls3_smad2enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_ls2_vs_ls3_smad3enriched <-  df_ls2_vs_ls3_smad3enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()


prom_genes_ls2_vs_ls3_bothenriched <-  df_ls2_vs_ls3_bothenriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_ls2_vs_ls3_smad2enriched

prom_genes_ls2_vs_ls3_smad3enriched

prom_genes_ls2_vs_ls3_bothenriched


intersect(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_smad3enriched$GENENAME) %>% # finding common peaks between LAP-Smad2 and LAP-Smad3
  tibble::enframe(name = "No.", value = "Promoter-proximal peak annotations") %>% 
  tidyr::drop_na()

intersect(prom_genes_ls2_vs_ls3_smad2enriched$GENENAME, prom_genes_ls2_vs_ls3_bothenriched$GENENAME) %>% 
  tibble::enframe(name = "No.", value = "Promoter-proximal peak annotations") %>% 
  tidyr::drop_na()

intersect(prom_genes_ls2_vs_ls3_smad3enriched$GENENAME, prom_genes_ls2_vs_ls3_bothenriched$GENENAME) %>% 
  tibble::enframe(name = "No.", value = "Promoter-proximal peak annotations") %>% 
  tidyr::drop_na()

#### Collecting a list of "ambiguous" genes

ambi_ls2_vs_ls3 <- dplyr::union(intersect(prom_genes_ls2_vs_ls3_smad2enriched$geneId, prom_genes_ls2_vs_ls3_bothenriched$geneId), intersect(prom_genes_ls2_vs_ls3_smad3enriched$geneId, prom_genes_ls2_vs_ls3_bothenriched$geneId)) %>% 
  dplyr::union(., intersect(prom_genes_ls2_vs_ls3_smad2enriched$geneId, prom_genes_ls2_vs_ls3_smad3enriched$geneId))

ego_prom_genes_ls2_vs_ls3_smad2enriched <- enrichGO(gene = (prom_genes_ls2_vs_ls3_smad2enriched %>% filter(!geneId %in% ambi_ls2_vs_ls3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)

ego_prom_genes_ls2_vs_ls3_smad3enriched <- enrichGO(gene = (prom_genes_ls2_vs_ls3_smad3enriched %>% filter(!geneId %in% ambi_ls2_vs_ls3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)

ego_prom_genes_ls2_vs_ls3_smad2enriched@result %>% 
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID)


ego_prom_genes_ls2_vs_ls3_smad3enriched@result %>% 
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID)

################### vs native

bed_ls2_vs_ns3_smad2enriched <- read_tsv("./differential_binding/untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_native-smad3_c3.0_cond1.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ns3_smad2enriched <- GRanges(seqnames = bed_ls2_vs_ns3_smad2enriched$chrom, ranges = IRanges(start = bed_ls2_vs_ns3_smad2enriched$start, end = bed_ls2_vs_ns3_smad2enriched$end), log10_LR = bed_ls2_vs_ns3_smad2enriched$log10_likelihood_ratio)

ranges(ls2_vs_ns3_smad2enriched)

bed_ls2_vs_ns3_smad3enriched <- read_tsv("./differential_binding/untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_native-smad3_c3.0_cond2.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ns3_smad3enriched <- GRanges(seqnames = bed_ls2_vs_ns3_smad3enriched$chrom, ranges = IRanges(start = bed_ls2_vs_ns3_smad3enriched$start, end = bed_ls2_vs_ns3_smad3enriched$end), log10_LR = bed_ls2_vs_ns3_smad3enriched$log10_likelihood_ratio)

bed_ls2_vs_ns3_bothenriched <- read_tsv("./differential_binding/untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_native-smad3_c3.0_common.bed", skip = 1, col_names = c("chrom", "start", "end", "name", "log10_likelihood_ratio"))

ls2_vs_ns3_bothenriched <- GRanges(seqnames = bed_ls2_vs_ns3_bothenriched$chrom, ranges = IRanges(start = bed_ls2_vs_ns3_bothenriched$start, end = bed_ls2_vs_ns3_bothenriched$end), log10_LR = bed_ls2_vs_ns3_bothenriched$log10_likelihood_ratio)


### Annotating peaks

TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene ### loading tx database for annotation

ap_ls2_vs_ns3_smad2enriched <- annotatePeak(peak = ls2_vs_ns3_smad2enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_ls2_vs_ns3_smad3enriched <- annotatePeak(peak = ls2_vs_ns3_smad3enriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

ap_ls2_vs_ns3_bothenriched <- annotatePeak(peak = ls2_vs_ns3_bothenriched, TxDb = TxDb, annoDb = "org.Hs.eg.db")

### Creating tibbles from peak annotation data

df_ls2_vs_ns3_smad2enriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ns3_smad2enriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())


df_ls2_vs_ns3_smad3enriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ns3_smad3enriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())


df_ls2_vs_ns3_bothenriched <- AnnotationDbi::as.data.frame(ap_ls2_vs_ns3_bothenriched) %>% as_tibble() %>% 
  dplyr::select(seqnames, start, end, width, annotation, GENENAME, log10_LR, SYMBOL, geneId, ENSEMBL, everything())


### Extracting unique list of genes containing promoter-proximal peaks

prom_genes_ls2_vs_ns3_smad2enriched <- df_ls2_vs_ns3_smad2enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()

prom_genes_ls2_vs_ns3_smad3enriched <-  df_ls2_vs_ns3_smad3enriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()


prom_genes_ls2_vs_ns3_bothenriched <-  df_ls2_vs_ns3_bothenriched %>% filter(str_detect(string = annotation, pattern = "Promoter")) %>% 
  dplyr::select(GENENAME, geneId) %>% 
  distinct(geneId, .keep_all = TRUE) %>% 
  tidyr::drop_na()


prom_genes_ls2_vs_ns3_smad2enriched

prom_genes_ls2_vs_ns3_smad3enriched

prom_genes_ls2_vs_ns3_bothenriched

### Setops

#### Genes common to both the Smad2-differentially enriched peak set and the Smad3-differentially enriched peak set

intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_smad3enriched$GENENAME) %>% # finding common peaks between LAP-Smad2 and LAP-Smad3
  tibble::enframe(name = "No.", value = "Promoter-proximal peak annotations") %>% 
  tidyr::drop_na()

#### Genes common to both the Smad2-differentially enriched peak set and the non-differentially expressed peak set

intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME) %>% # finding common peaks between LAP-Smad2 and LAP-Smad3
  tibble::enframe(name = "No.", value = "Promoter-proximal peak annotations") %>% 
  tidyr::drop_na()

#### Genes common to both the Smad3-differentially enriched peak set and the non-differentially expressed peak set

intersect(prom_genes_ls2_vs_ns3_smad3enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME) %>% # finding common peaks between LAP-Smad2 and LAP-Smad3
  tibble::enframe(name = "No.", value = "Promoter-proximal peak annotations") %>% 
  tidyr::drop_na()

#### Collecting a list of "ambiguous" genes


dplyr::union(intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_bothenriched$GENENAME), 
             intersect(prom_genes_ls2_vs_ns3_smad3enriched$GENENAME,prom_genes_ls2_vs_ns3_bothenriched$GENENAME)) %>% 
  dplyr::union(., intersect(prom_genes_ls2_vs_ns3_smad2enriched$GENENAME, prom_genes_ls2_vs_ns3_smad3enriched$GENENAME)) %>% 
  tibble::enframe(name = "No.", value = "'Ambiguous' genes") %>% 
  tidyr::drop_na()

ambi_ls2_vs_ns3 <- dplyr::union(intersect(prom_genes_ls2_vs_ns3_smad2enriched$geneId, prom_genes_ls2_vs_ns3_bothenriched$geneId), intersect(prom_genes_ls2_vs_ns3_smad3enriched$geneId, prom_genes_ls2_vs_ns3_bothenriched$geneId)) %>% 
  dplyr::union(., intersect(prom_genes_ls2_vs_ns3_smad2enriched$geneId, prom_genes_ls2_vs_ns3_smad3enriched$geneId))

ambi_ls2_vs_ns3

### GO term enrichment analysis


#### GO terms enriched in Smad2-exclusive differentially enriched peaks

prom_genes_ls2_vs_ns3_smad2enriched %>% filter(!geneId %in% ambi_ls2_vs_ns3) # only one gene = smad2


ego_prom_genes_ls2_vs_ns3_smad2enriched <- enrichGO(gene = (prom_genes_ls2_vs_ns3_smad2enriched %>% filter(!geneId %in% ambi_ls2_vs_ns3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)


ego_prom_genes_ls2_vs_ns3_smad2enriched@result %>% 
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID)

prom_genes_ls2_vs_ns3_smad3enriched %>% filter(!geneId %in% ambi_ls2_vs_ns3)


ego_prom_genes_ls2_vs_ns3_smad3enriched <- enrichGO(gene = (prom_genes_ls2_vs_ns3_smad3enriched %>% filter(!geneId %in% ambi_ls2_vs_ns3) %>% dplyr::select(geneId) %>% unlist(use.names = FALSE)), org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)



ego_prom_genes_ls2_vs_ns3_smad3enriched@result %>%
  as_tibble() %>% 
  arrange(p.adjust, desc(Count) , desc(GeneRatio))  %>% 
  filter(p.adjust < 0.05) %>%  
  dplyr::select(Description, p.adjust, GeneRatio, BgRatio, geneID) %>% View()

save(list = ls(), file = "untreated-differential_binding.RData")

save(list = c("ego_prom_genes_ls2_vs_ls3_smad2enriched", "ego_prom_genes_ls2_vs_ls3_smad3enriched", "ego_prom_genes_ls2_vs_ns3_smad2enriched", "ego_prom_genes_ls2_vs_ns3_smad3enriched"), file  = "untreated-differential_binding_go_term.RData")