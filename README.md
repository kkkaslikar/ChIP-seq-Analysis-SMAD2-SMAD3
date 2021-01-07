
# Introduction

This repository contains the code and some output files associated with the my project which attempted to find the difference between the Smad2 and the Smad3 axes of the TGFβ pathway.

ChIP-seq had been performed for Smad2 and Smad3 before and after TGF-beta treatment. Peak-calling and differential peak-calling was been performed on ChIP-seq peaks.

To see the a more detailed explanation of the analysis, as well as the rendered Rmarkdown notebooks, please visit the rendered GitHub pages for this repository at [this link](https://kkkaslikar.github.io/ChIP-seq-Analysis-SMAD2-SMAD3/)

# Files included in the repository

Due to size constraints on Github, many of the intermediate output files have not been uploaded to this repository. My intention was primarily to convey essential information such as the parameters and options used for the commands in the pipeline, as well as the R code used to perform the analysis (which is present in the form of Rmarkdown files). Hence, I have not uploaded very large data files, such as the BAM files for the aligned ChIP-seq reads, as well as the bedgraph pileup files used for differential expression analysis. Instead, in the relevant directories where these files are missing (as well as in some other directories), there will be files called labelled "localtree.txt" (derived using the Linux `tree` utility), which will include the filename of the missing files. This is done only for reference, so as to allow better interpretation of the code which refers to those files, especially with regards to the location of the file with respect to the directory in which a particular shell command is being invoked. 
