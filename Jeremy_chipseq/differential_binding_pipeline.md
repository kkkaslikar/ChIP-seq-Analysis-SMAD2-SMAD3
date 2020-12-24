

[TOC]



# Metadata

| Sample ID                    |  Antibody Type          | TGF-beta Treatment |
|------------------------------|-------------------------|--------------------|
| FCC6CLFACXX-wHAPSI016408-113 | LAP Smad2 untreated     | Untreated          |
| FCC6CLFACXX-wHAPSI016409-133 | LAP Smad2 treated       | Treated            |
| S3o_GFP                      | LAP Smad3 treated       | Treated            |
| S3xT_GFP                     | LAP Smad3 untreated     | Untreated          |
| MDAo_Smad3                   | Native Smad3 treated    | Treated            |
| MDAxT_Smad3                  |  Native Smad3 untreated | Untreated          |


# Creating a directory for storing files for differential binding

```bash
 mkdir differential_binding
```

# Treated sample differential peak-calling

## Getting uniform extension size for running callpeak using `predictd`

`macs2 predictd -i ./bamfiles/MDA-treated_Smad3.bam`

`macs2 predictd -i ./bamfiles/LAP-treated_Smad3.bam`

`macs2 predictd -i ./bamfiles/FCC6CLFACXX-wHAPSI016409-133_L5_1.bam`

### Output for predictd

#### MDA-treated_Smad3.bam


```
tag size is determined as 50 bps
INFO  @ Thu, 05 Sep 2019 14:21:33: # tag size = 50
INFO  @ Thu, 05 Sep 2019 14:21:33: # total tags in alignment file: 42588774
INFO  @ Thu, 05 Sep 2019 14:21:33: # Build Peak Model...
INFO  @ Thu, 05 Sep 2019 14:21:33: #2 looking for paired plus/minus strand peaks...
INFO  @ Thu, 05 Sep 2019 14:21:57: #2 number of paired peaks: 31741
INFO  @ Thu, 05 Sep 2019 14:21:57: start model_add_line...
INFO  @ Thu, 05 Sep 2019 14:22:00: start X-correlation...
INFO  @ Thu, 05 Sep 2019 14:22:00: end of X-cor
INFO  @ Thu, 05 Sep 2019 14:22:00: # finished!
INFO  @ Thu, 05 Sep 2019 14:22:00: # predicted fragment length is 311 bps
INFO  @ Thu, 05 Sep 2019 14:22:00: # alternative fragment length(s) may be 311,364 bps
INFO  @ Thu, 05 Sep 2019 14:22:00: # Generate R script for model : predictd
```

#### LAP-treated_Smad3.bam

```
INFO  @ Thu, 05 Sep 2019 14:27:09: tag size is determined as 50 bps
INFO  @ Thu, 05 Sep 2019 14:27:09: # tag size = 50
INFO  @ Thu, 05 Sep 2019 14:27:09: # total tags in alignment file: 39195161
INFO  @ Thu, 05 Sep 2019 14:27:09: # Build Peak Model...
INFO  @ Thu, 05 Sep 2019 14:27:09: #2 looking for paired plus/minus strand peaks...
INFO  @ Thu, 05 Sep 2019 14:27:34: #2 number of paired peaks: 78611
INFO  @ Thu, 05 Sep 2019 14:27:34: start model_add_line...
INFO  @ Thu, 05 Sep 2019 14:27:36: start X-correlation...
INFO  @ Thu, 05 Sep 2019 14:27:36: end of X-cor
INFO  @ Thu, 05 Sep 2019 14:27:36: # finished!
INFO  @ Thu, 05 Sep 2019 14:27:36: # predicted fragment length is 195 bps
INFO  @ Thu, 05 Sep 2019 14:27:36: # alternative fragment length(s) may be 195 bps
INFO  @ Thu, 05 Sep 2019 14:27:36: # Generate R script for model : predictd
```

#### FCC6CLFACXX-wHAPSI016409-133_L5_1.bam

```
INFO  @ Thu, 05 Sep 2019 14:35:20: tag size is determined as 49 bps
INFO  @ Thu, 05 Sep 2019 14:35:20: # tag size = 49
INFO  @ Thu, 05 Sep 2019 14:35:20: # total tags in alignment file: 101380397
INFO  @ Thu, 05 Sep 2019 14:35:20: # Build Peak Model...
INFO  @ Thu, 05 Sep 2019 14:35:20: #2 looking for paired plus/minus strand peaks...
INFO  @ Thu, 05 Sep 2019 14:35:34: #2 number of paired peaks: 6218
INFO  @ Thu, 05 Sep 2019 14:35:34: start model_add_line...
INFO  @ Thu, 05 Sep 2019 14:35:35: start X-correlation...
INFO  @ Thu, 05 Sep 2019 14:35:35: end of X-cor
INFO  @ Thu, 05 Sep 2019 14:35:35: # finished!
INFO  @ Thu, 05 Sep 2019 14:35:35: # predicted fragment length is 150 bps
INFO  @ Thu, 05 Sep 2019 14:35:35: # alternative fragment length(s) may be 150,492 bps
INFO  @ Thu, 05 Sep 2019 14:35:35: # Generate R script for model : predictd
```



Based on the predictd output, the uniform extension size to be used for peak calling (fed as the `--extsize` argument to `callpeak` ) was calculated to be the average of the predicted fragment lengths of the 3 bam files:

(195 + 150 + 311) / 3 = 218


## Calling peaks to generate BedGraph files

`macs2 callpeak -B -t bamfiles/MDA-treated_Smad3.bam -c bamfiles/input.bam -n treated-native-smad3 --nomodel --extsize 218`

`macs2 callpeak -B -t bamfiles/LAP-treated_Smad3.bam -c bamfiles/input.bam -n treated-lap-smad3 --nomodel --extsize 218`

`macs2 callpeak -B -t bamfiles/FCC6CLFACXX-wHAPSI016409-133_L5_1.bam -c bamfiles/input.bam -n treated-lap-smad2 --nomodel --extsize 218`

## Moving peak-calling results to another directory

```bash
mkdir differential_binding/treated
mv treated-lap-smad2/ differential_binding/treated/
mv treated-lap-smad3/ differential_binding/treated/
mv treated-native-smad3 differential_binding/treated/
```


## Tags after filtering

### lap-smad2 tags after filtering

tags after filtering in treatment: 43345452
tags after filtering in control: 31588061

### native-smad3 tags after filtering

tags after filtering in treatment: 19110302
tags after filtering in control: 31588061

### lap-smad3 tags after filtering

tags after filtering in treatment: 10967473
tags after filtering in control: 31588061


## Calling differential peaks

### lap-smad2 vs lap-smad3 differential peaks



```bash
macs2 bdgdiff --t1 differential_binding/treated/treated-lap-smad2/lap-smad2_treat_pileup.bdg --c1 differential_binding/treated/treated-lap-smad2/lap-smad2_control_lambda.bdg --t2 differential_binding/treated/treated-lap-smad3/lap-smad3_treat_pileup.bdg --c2 differential_binding/treated/treated-lap-smad2/lap-smad3_control_lambda.bdg --d1 31588061 --d2 10967473 --o-prefix differential_binding/treated/diffpeak_results/diff_lap-smad2_vs_lap-smad3
```



### lap-smad2 vs native-smad3 differential peaks


```bash
macs2 bdgdiff --t1 differential_binding/treated/treated-lap-smad2/lap-smad2_treat_pileup.bdg --c1 differential_binding/treated/treated-lap-smad2/lap-smad2_control_lambda.bdg --t2 differential_binding/treated/treated-native-smad3/native-smad3_treat_pileup.bdg --c2 differential_binding/treated/treated-native-smad3/native-smad3_control_lambda.bdg --d1 31588061 --d2 19110302 --o-prefix differential_binding/treated/diffpeak_results/diff_lap-smad2_vs_native-smad3
```

---



# Untreated sample differential peak-calling
	
## Getting uniform extension size for running callpeak using `predictd`

macs2 predictd -i ./bamfiles/LAP-untreated_Smad3.bam |tee -a predictd.txt
macs2 predictd -i ./bamfiles/MDA-untreated_Smad3.bam |tee -a predictd.txt
macs2 predictd -i ./bamfiles/FCC6CLFACXX-wHAPSI016408-113_L5_1.bam |tee -a predictd.txt

### Output for predictd

#### MDA-untreated_Smad3.bam


```
tag size is determined as 50 bps
INFO  @ Wed, 23 Oct 2019 20:50:22: # tag size = 50
INFO  @ Wed, 23 Oct 2019 20:50:22: # total tags in alignment file: 39023219
INFO  @ Wed, 23 Oct 2019 20:50:22: # Build Peak Model...
INFO  @ Wed, 23 Oct 2019 20:50:22: #2 looking for paired plus/minus strand peaks...
INFO  @ Wed, 23 Oct 2019 20:50:46: #2 number of paired peaks: 28780
INFO  @ Wed, 23 Oct 2019 20:50:46: start model_add_line...
INFO  @ Wed, 23 Oct 2019 20:50:47: start X-correlation...
INFO  @ Wed, 23 Oct 2019 20:50:47: end of X-cor
INFO  @ Wed, 23 Oct 2019 20:50:47: # finished!
INFO  @ Wed, 23 Oct 2019 20:50:47: # predicted fragment length is 258 bps
INFO  @ Wed, 23 Oct 2019 20:50:47: # alternative fragment length(s) may be 168,258 bps
INFO  @ Wed, 23 Oct 2019 20:50:47: # Generate R script for model : predictd
```

#### LAP-untreated_Smad3.bam

```
tag size is determined as 50 bps
INFO  @ Wed, 23 Oct 2019 20:27:46: # tag size = 50
INFO  @ Wed, 23 Oct 2019 20:27:46: # total tags in alignment file: 34934067
INFO  @ Wed, 23 Oct 2019 20:27:46: # Build Peak Model...
INFO  @ Wed, 23 Oct 2019 20:27:46: #2 looking for paired plus/minus strand peaks...
INFO  @ Wed, 23 Oct 2019 20:28:26: #2 number of paired peaks: 84007
INFO  @ Wed, 23 Oct 2019 20:28:26: start model_add_line...
INFO  @ Wed, 23 Oct 2019 20:28:29: start X-correlation...
INFO  @ Wed, 23 Oct 2019 20:28:29: end of X-cor
INFO  @ Wed, 23 Oct 2019 20:28:29: # finished!
INFO  @ Wed, 23 Oct 2019 20:28:29: # predicted fragment length is 177 bps
INFO  @ Wed, 23 Oct 2019 20:28:29: # alternative fragment length(s) may be 118,177 bps
INFO  @ Wed, 23 Oct 2019 20:28:29: # Generate R script for model : predictd
```

#### FCC6CLFACXX-wHAPSI016408-113_L5_1.bam

```
tag size is determined as 49 bps
INFO  @ Wed, 23 Oct 2019 21:12:07: # tag size = 49
INFO  @ Wed, 23 Oct 2019 21:12:07: # total tags in alignment file: 61101120
INFO  @ Wed, 23 Oct 2019 21:12:07: # Build Peak Model...
INFO  @ Wed, 23 Oct 2019 21:12:07: #2 looking for paired plus/minus strand peaks...
INFO  @ Wed, 23 Oct 2019 21:12:28: #2 number of paired peaks: 2348
INFO  @ Wed, 23 Oct 2019 21:12:28: start model_add_line...
INFO  @ Wed, 23 Oct 2019 21:12:30: start X-correlation...
INFO  @ Wed, 23 Oct 2019 21:12:31: end of X-cor
INFO  @ Wed, 23 Oct 2019 21:12:31: # finished!
INFO  @ Wed, 23 Oct 2019 21:12:31: # predicted fragment length is 48 bps
INFO  @ Wed, 23 Oct 2019 21:12:31: # alternative fragment length(s) may be 48,216,244,295,381,458,541 bps
INFO  @ Wed, 23 Oct 2019 21:12:31: # Generate R script for model : predictd
```



Based on the predictd output, the uniform extension size to be used for peak calling (fed as the `--extsize` argument to `callpeak` ) was calculated to be the average of the predicted fragment lengths of the 3 bam files. Since the predicted fragment length for  FCC6CLFACXX-wHAPSI016408-113_L5_1.bam was deemed to be too low,  the first of the alternative predicted fragment lengths, which is 216, was used instead.

(216 + 177 + 258) / 3 = 217



## Calling peaks to generate BedGraph files

`macs2 callpeak -B -t bamfiles/MDA-untreated_Smad3.bam -c bamfiles/input.bam -n untreated-native-smad3 --nomodel --extsize 217`

`macs2 callpeak -B -t bamfiles/LAP-untreated_Smad3.bam -c bamfiles/input.bam -n untreated-lap-smad3 --nomodel --extsize 217`

`macs2 callpeak -B -t bamfiles/FCC6CLFACXX-wHAPSI016408-113_L5_1.bam -c bamfiles/input.bam -n untreated-lap-smad2 --nomodel --extsize 217`

## Moving peak-calling results to another directory

```bash
mkdir differential_binding/treated
mv untreated-lap-smad2/ differential_binding/untreated/
mv untreated-lap-smad3/ differential_binding/untreated/
mv untreated-native-smad3 differential_binding/untreated/
```

## Tags after filtering

### lap-smad2 tags after filtering

tags after filtering in treatment: 46760825
tags after filtering in control: 31588061

### native-smad3 tags after filtering

tags after filtering in treatment: 13426831
tags after filtering in control: 31588061

### lap-smad3 tags after filtering

tags after filtering in treatment: 6675026
tags after filtering in control: 31588061


## Calling differential peaks

### lap-smad2 vs lap-smad3 differential peaks



```bash
macs2 bdgdiff --t1 differential_binding/untreated/untreated-lap-smad2/untreated-lap-smad2_treat_pileup.bdg --c1 differential_binding/untreated/untreated-lap-smad2/untreated-lap-smad2_control_lambda.bdg --t2 differential_binding/untreated/untreated-lap-smad3/untreated-lap-smad3_treat_pileup.bdg --c2 differential_binding/untreated/untreated-lap-smad3/untreated-lap-smad3_control_lambda.bdg --d1 31588061 --d2 31588061 --o-prefix differential_binding/untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_lap-smad3
```

### lap-smad2 vs native-smad3 differential peaks



```bash
macs2 bdgdiff --t1 differential_binding/untreated/untreated-lap-smad2/untreated-lap-smad2_treat_pileup.bdg --c1 differential_binding/untreated/untreated-lap-smad2/untreated-lap-smad2_control_lambda.bdg --t2 differential_binding/untreated/untreated-native-smad3/untreated-native-smad3_treat_pileup.bdg --c2 differential_binding/untreated/untreated-native-smad3/untreated-native-smad3_control_lambda.bdg --d1 31588061 --d2 13426831 --o-prefix differential_binding/untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_native-smad3
```
