
# Quick access links

This section provides quick access to the analysis pages mentioned in the subsequent sections, all collected in one place.

* [Preliminary peak annotation and GO term analysis](https://kkkaslikar.github.io/ChIP-seq-Analysis-SMAD2-SMAD3/Jeremy_chipseq/peak_annotation_go_term_analysis.nb.html)
* [GO term analysis for TGFβ-treated differential peaks](https://kkkaslikar.github.io/ChIP-seq-Analysis-SMAD2-SMAD3/Jeremy_chipseq/differential_binding/treated/treated-differential_binding.nb.html)
* [Code for performing the differential binding analysis](https://kkkaslikar.github.io/ChIP-seq-Analysis-SMAD2-SMAD3/Jeremy_chipseq/differential_binding_pipeline.html) 


# Introduction

Smad2 and Smad3 are both DNA-binding proteins that are activated by TGFβ signalling, and have differing downstream targets and phenotypic effects. The predominance of either the Smad2 or the Smad3 pathway over the other can determine the overall effect of signalling TGFβ . The aim of this project was to determine the difference between the Smad2 and the Smad3 axes of the TGFβ pathway. ChIP-seq reads from antibodies targeting Smad2 and Smad3 were available, and had already been aligned using bowtie2, giving BAM files, which were used as a starting point for this analysis.

ChIP-seq had been performed for Smad2 and Smad3 before and after TGF-beta treatment. In addition, two different antibodies had been used to perform the ChIP-seq for Smad3: an antibody targeting the GFP portion of the Smad3 fusion protein, and an antibody targeting an epitope found within native Smad3. Hence, in total, there were at least six separate conditions, which have been listed in the table below:

| Sample ID                    | Antibody Type                                           | TGFβ Treatment | Sample Reference Name  |
| ---------------------------- | ------------------------------------------------------- | -------------- | ---------------------- |
| FCC6CLFACXX-wHAPSI016408-113 | anti-GFP (targeting GFP portion of Smad2)               | Untreated      | LAP Smad2 Untreated    |
| FCC6CLFACXX-wHAPSI016409-133 | Anti-GFP (targeting GFP portion of Smad2)               | Treated        | LAP Smad2 Treated      |
| S3xT_GFP                     | Anti-GFP (targeting GFP portion of Smad3)               | Untreated      | LAP Smad3 Untreated    |
| S3o_GFP                      | Anti-GFP (targeting GFP portion of Smad3)               | Treated        | LAP Smad3 Treated      |
| MDAxT_Smad3                  | Anti-Smad3 (targeting epitope in native Smad3 protein ) | Untreated      | Native Smad3 Untreated |
| MDAo_Smad3                   | Anti-Smad3 (targeting epitope in native Smad3 protein ) | Treated        | Native Smad3 Untreated |


It was determined that identifying the genes in which the ChIP-seq peaks fell, and performing a GO term analysis for these genes would be a good method to find out the differences between the Smad2 and Smad3 axes.

# Preliminary analysis of peak annotation and GO terms

As a preliminary analysis, peaks were first called for the treated and untreated LAP Smad2 and LAP Smad3 data using MACS2, in order to better explore the output of the peak annotation and to get a sense of the associated Gene Ontology terms. Please refer to [this link](https://kkkaslikar.github.io/ChIP-seq-Analysis-SMAD2-SMAD3/Jeremy_chipseq/peak_annotation_go_term_analysis.nb.html) for this analysis.

# Differential binding analysis

We decided that performing the analysis on *differentially*-called peaks would be a better way of appreciating differences between the Smad2 and Smad3 axes, rather than just looking at regularly-called peaks, since Smad2 and Smad3 share some binding sites. A differential peak would indicate more promiscous binding of one over the other at a particular location, and might indicate that a particular gene is more influenced by, for example, Smad3 rather than Smad2. Hence, such a gene would be part of the Smad3 axis rather than the Smad2 axis.

Differential binding analysis was performed on both the treated and untreated Smad2 and Smad3 data, for both the native as well as the GFP-targeted Smad3 antibody. The code used for the differential binding, as well as some of the output, is presented at [this link](https://kkkaslikar.github.io/ChIP-seq-Analysis-SMAD2-SMAD3/Jeremy_chipseq/differential_binding_pipeline.html)

## Analysis of differential peaks

Priority was given to the analysis of differential peaks from treated ChIP-seq samples, since our objective was to understand the difference between the Smad2 and Smad3 binding sites *after* TGFβ treatment. The rendered R notebook for this analysis is given at [this link](https://kkkaslikar.github.io/ChIP-seq-Analysis-SMAD2-SMAD3/Jeremy_chipseq/differential_binding/treated/treated-differential_binding.nb.html)
