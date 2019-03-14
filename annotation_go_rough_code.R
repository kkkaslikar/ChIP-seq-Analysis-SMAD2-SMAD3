
# Creating GRanges from MACS2 output ----------------------------------


library(tidyverse)
library(GenomicRanges)




## Takes a filename as string input, and if it has suffix ".narrowPeak", converts it to a GRanges object containing a q_value metadata column. Requires tidyverse and GRanges packages. 

narrow_to_granges <- function(file){stopifnot(file.exists(file), endsWith(x = file, suffix = ".narrowPeak")) ## checking for existence of file and correct extension
  
  bed <- as_tibble(read_tsv(file = file, col_names = c("chrom", "start", "end", "name", "score", "strand", "signal_value", "p_value", "q_value", "peak"))) 
  
  
  gr <- GRanges(seqnames = bed$chrom, ranges = IRanges(start = bed$start, end = bed$end), strand = if_else(condition = bed$strand == "+"|  bed$strand == "-", true = bed$strand, false = "*"), q_value = bed$q_value)
  
  return(gr)}



### To read in multiple .narrowPeak files and create a GRanges list of GRanges objects with the input of respective files.   

files <- dir(path = ".", pattern = "\\.narrowPeak$") ## specify files in character vector

grl <- lapply(files, narrow_to_granges)




# Annotate Peaks in GRanges -----------------------------------------------

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)


TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene

v <- str_remove_all(string = files, pattern = "\\.narrowPeak$") ## vector of variable names for granges in list

### Assign individual granges within grl to variable names in v. 

for(i in seq_along(grl)){
  assign(x = v[i], value = grl[[i]])
}

ap_v <- paste0("AP_", v) ## variable names for annotated peak objects


### Create annotated peak objects and assign them to names in ap_v
for(i in seq_along(grl)){
  assign(x = ap_v[i], value = annotatePeak(peak = grl[[i]], TxDb = TxDb, annoDb = "org.Hs.eg.db"))
}


ap_v_df <- paste0("df_", ap_v) ## variable names for data frames derived from annotated peaks


### Create tibbles from annotated peak objects listed in ap_v and assign them to names in ap_v_df

for(i in seq_along(ap_v)){
  assign(x = ap_v_df[i], value = AnnotationDbi::as.data.frame(eval(parse(text = ap_v[i]))) %>% as_tibble) 
}


ap_v_prom_df <- paste0("prom_df_", ap_v)

for(i in seq_along(ap_v_df)){
  assign(x = ap_v_prom_df[i], value = parse(text = ap_v_df[i]) %>% 
           eval() %>% 
           filter(str_detect(annotation, "Promoter")))
}

ap_v_gen_df <- paste0("gen_df_", ap_v) 

for(i in seq_along(ap_v_df)){
  assign(x = ap_v_gen_df[i], value = parse(text = ap_v_df[i]) %>% 
           eval() %>% 
           filter(!str_detect(annotation, "Distal Intergenic|Downstream")))
}




plotAnnoPie(AP_SMAD2_abInput_treated_peaks)
title(main = "AP_SMAD2_abInput_treated_peaks", line = -2, adj = 0)

plotAnnoPie(AP_SMAD2_abInput_untreated_peaks)
title(main = "AP_SMAD2_abInput_untreated_peaks", line = -2, adj = 0)

plotAnnoPie(AP_SMAD3_LAP_untreated_peaks)
title(main = "AP_SMAD3_LAP_untreated_peaks", line = -2, adj = 0)

plotAnnoPie(AP_SMAD3_LAP_treated_peaks)
title(main = "AP_SMAD3_LAP_treated_peaks", line = -2, adj = 0)



annobar_SMAD2_abInput_untreated_peaks <- plotAnnoBar(AP_SMAD2_abInput_untreated_peaks) + ggtitle(label = "SMAD2 Untreated Peaks") + theme(plot.title = element_text(hjust = 0.5))


annobar_SMAD2_abInput_treated_peaks <- plotAnnoBar(AP_SMAD2_abInput_untreated_peaks) + ggtitle(label = "SMAD2 Treated Peaks") + theme(plot.title = element_text(hjust = 0.5))

annobar_SMAD3_LAP_untreated_peaks <- plotAnnoBar(AP_SMAD2_abInput_untreated_peaks) + ggtitle(label = "SMAD3 Untreated Peaks") + theme(plot.title = element_text(hjust = 0.5))


annobar_SMAD3_LAP_treated_peaks <- plotAnnoBar(AP_SMAD2_abInput_untreated_peaks) + ggtitle(label = "SMAD3 Treated Peaks") + theme(plot.title = element_text(hjust = 0.5))


gridExtra::grid.arrange(annobar_SMAD2_abInput_untreated_peaks, annobar_SMAD2_abInput_treated_peaks, annobar_SMAD3_LAP_untreated_peaks, annobar_SMAD3_LAP_treated_peaks)





# GO-term Analysis Trial --------------------------------------------------


library(clusterProfiler)


ego_v <- paste0("ego_", v) ### Vector of variable names for go term enrichment objects

### Performing go term analysis using the clusterprofiler package and assigning the results to the variables in ego_v

for(i in seq_along(ego_v)){ 
  ego <- enrichGO(gene = eval(parse(text = (ap_v_df[i])))[["geneId"]], OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)
  if(!nrow(as.data.frame(ego)) == 0) 
    assign(x = ego_v[i], value = ego)}






### Creating dotplots of all the annotated peaks 

## No enriched terms for SMAD2 untreated peaks. 


dotplot(ego_SMAD2_abInput_treated_peaks) + ggtitle(label = "Enriched Biological Process GO Term in SMAD2 Treated Peaks") + theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text.y = element_text(size = 10))

dotplot(ego_SMAD3_LAP_untreated_peaks) + ggtitle(label = "Enriched Biological Process GO Term in SMAD3 Untreated Peaks") + theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text.y = element_text(size = 10))

dotplot(ego_SMAD3_LAP_treated_peaks) + ggtitle(label = "Enriched Biological Process GO Terms in SMAD3 Treated Peaks") + theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text.y = element_text(size = 10))





### 

prom_ego_v <- ego_v <- paste0("prom_ego_", v) ### Vector of variable names for go term enrichment objects limited to promoter peaks



for(i in seq_along(prom_ego_v)){ 
  ego <- enrichGO(gene = eval(parse(text = (ap_v_prom_df[i])))[["geneId"]], OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)
  if(!nrow(as.data.frame(ego)) == 0) 
    assign(x = prom_ego_v[i], value = ego)}





### 

dotplot(prom_ego_SMAD2_abInput_untreated_peaks) + ggtitle(label = "Enriched BP GO Terms in Untreated SMAD2 Promoter Peaks") + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_text(size = 7))

dotplot(prom_ego_SMAD3_LAP_untreated_peaks) + ggtitle(label = "Enriched BP GO Terms in Untreated SMAD3 Promoter Peaks") + theme(plot.title = element_text(hjust = 0.8, size = 15), axis.text.y = element_text(size = 7))

dotplot(prom_ego_SMAD3_LAP_treated_peaks) + ggtitle(label = "Enriched BP GO Term in Treated SMAD3 Promoter Peaks") + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_text(size = 7))


### 

gen_ego_v <- paste0("gen_ego_", v) ### Vector of variable names for go term enrichment objects 


for(i in seq_along(gen_ego_v)){ 
  ego <- enrichGO(gene = eval(parse(text = (ap_v_gen_df[i])))[["geneId"]], OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)
  if(!nrow(as.data.frame(ego))==0) 
    assign(x = gen_ego_v[i], value = ego)}



dotplot(gen_ego_SMAD2_abInput_treated_peaks) + ggtitle(label = "Enriched BP GO Terms in Treated SMAD2 Genic Peaks") + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_text(size = 7))

dotplot(gen_ego_SMAD3_LAP_untreated_peaks )+ ggtitle(label = "Enriched BP GO Terms in Untreated SMAD3 Genic Peaks") + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_text(size = 7))

dotplot(gen_ego_SMAD3_LAP_treated_peaks )+ ggtitle(label = "Enriched BP GO Terms in Treated SMAD3 Genic Peaks") + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_text(size = 7))




anti_join(gen_ego_SMAD2_abInput_treated_peaks@result %>% filter(p.adjust<0.05) %>% dplyr::select(Description), ego_SMAD2_abInput_treated_peaks@result %>% filter(p.adjust<0.05) %>% dplyr::select(Description)) %>% rename(Description = "Different Terms")


ego_SMAD2_abInput_treated_peaks@result %>% as_tibble() %>% filter(p.adjust<0.05)

