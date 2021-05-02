############################################################################################################
############################################################################################################
### Diversity analysis for James et al. 2021, Virus Evolution, veab041, https://doi.org/10.1093/ve/veab041
### Author: Chase W. Nelson, cnelson@amnh.org, Academia Sinica, American Museum of Natural History

library(boot)
library(tidyverse)
library(patchwork)
library(RColorBrewer)


############################################################################################################
############################################################################################################
### FILTER VCFs

# IN: /SARS-CoV-2-South-Africa/vcfs_raw
# run_VCF_Filter_v3.py 5 .02 .05 .95 > run_VCF_Filter_v3.out
# NB: there are 109 VCF files (samples) here

# Placed filtered VCF files in: /SARS-CoV-2-South-Africa/vcfs_filtered/

# BEFORE
# IN: /SARS-CoV-2-South-Africa/vcfs_raw
# DO: ggrep -Ev "^#" *.vcf | gsed -E "s/:/\t/" > vcfs_raw_pooled.tsv

# IN: /SARS-CoV-2-South-Africa/vcfs_filtered
# DO: ggrep -Ev "^#" *.vcf | gsed -E "s/:/\t/" > vcfs_filtered_pooled.tsv


############################################################################################################
############################################################################################################
### IMPORT DATA BEFORE & AFTER FILTERING

### RAW
vcfs_raw_pooled <- read_tsv("/SARS-CoV-2-South-Africa/vcfs_raw_pooled.tsv", 
                            col_names = c("sample", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

# Extract AF
(vcfs_raw_pooled$AF <- vcfs_raw_pooled$INFO)
(vcfs_raw_pooled$AF <- str_replace(string = vcfs_raw_pooled$AF, pattern = "^.+AF=", replacement = ""))
(vcfs_raw_pooled$AF <- str_replace(string = vcfs_raw_pooled$AF, pattern = ";.+$", replacement = ""))
(vcfs_raw_pooled$AF <- as.numeric(vcfs_raw_pooled$AF))

# sample: remove ".lf.vcf" from sample name
vcfs_raw_pooled$sample <- str_replace(string = vcfs_raw_pooled$sample, pattern = ".lf.vcf", replacement = "")

# Add factor for raw
vcfs_raw_pooled$status <- "raw"


### FILTERED
vcfs_filtered_pooled <- read_tsv("/SARS-CoV-2-South-Africa/vcfs_filtered_pooled.tsv", 
                            col_names = c("sample", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

# Extract AF
(vcfs_filtered_pooled$AF <- vcfs_filtered_pooled$INFO)
(vcfs_filtered_pooled$AF <- str_replace(string = vcfs_filtered_pooled$AF, pattern = "^.+AF=", replacement = ""))
(vcfs_filtered_pooled$AF <- str_replace(string = vcfs_filtered_pooled$AF, pattern = ";.+$", replacement = ""))
(vcfs_filtered_pooled$AF <- as.numeric(vcfs_filtered_pooled$AF))

# sample: remove ".lf.vcf_filtered" from sample name
vcfs_filtered_pooled$sample <- str_replace(string = vcfs_filtered_pooled$sample, pattern = ".lf_filtered.vcf", replacement = "")

# Add factor for raw
vcfs_filtered_pooled$status <- "filtered"


### SUMMARIZE
# Combine
vcfs_ALL <- rbind(vcfs_raw_pooled, vcfs_filtered_pooled)

# Calculate MAF
(vcfs_ALL$MAF <- vcfs_ALL$AF)
(vcfs_ALL[vcfs_ALL$MAF > 0.5, ]$MAF <- (1 - vcfs_ALL[vcfs_ALL$MAF > 0.5, ]$MAF))

# MAF non-zero
vcfs_ALL$MAF_nonzero <- vcfs_ALL$MAF
(vcfs_ALL[vcfs_ALL$MAF_nonzero == 0, ]$MAF_nonzero <- NA)

# Summarize
(vcfs_ALL_SUMMARY <- vcfs_ALL %>% 
  group_by(sample, status) %>%
  summarize(
    count = n(),
    min_MAF = min(MAF_nonzero, na.rm = TRUE),
    median_MAF = median(MAF_nonzero, na.rm = TRUE),
    mean_MAF = mean(MAF_nonzero, na.rm = TRUE),
    max_MAF = max(MAF_nonzero, na.rm = TRUE)
  ))

# Summarize NONZERO MAF
(vcfs_ALL_SUMMARY_NONZEROMAF <- filter(vcfs_ALL, is.finite(MAF_nonzero)) %>% 
    group_by(sample, status) %>%
    summarize(
      count = n(),
      min_MAF = min(MAF_nonzero, na.rm = TRUE),
      median_MAF = median(MAF_nonzero, na.rm = TRUE),
      mean_MAF = mean(MAF_nonzero, na.rm = TRUE),
      max_MAF = max(MAF_nonzero, na.rm = TRUE)
    ))

# Inspect
#filter(vcfs_ALL, sample == "K002189_S23") # it's because ALL PASSINF are fixed non-REF

## Separate VCF RECORDS by raw/filtered
vcfs_ALL_RAW <- filter(vcfs_ALL, status == "raw")
vcfs_ALL_FILTERED <- filter(vcfs_ALL, status == "filtered")

vcfs_ALL_RAW_NONZEROMAF <- filter(vcfs_ALL, status == "raw",  is.finite(MAF_nonzero))
vcfs_ALL_FILTERED_NONZEROMAF <- filter(vcfs_ALL, status == "filtered",  is.finite(MAF_nonzero))

# Number of minor allele variable sites prefiltering
nrow(filter(vcfs_ALL_RAW_NONZEROMAF)) # 14642
nrow(filter(vcfs_ALL_FILTERED_NONZEROMAF)) # 2129

## Separate VCF SUMMARIES by raw/filtered, NON ZERO MAF ONLY
vcfs_ALL_SUMMARY_RAW_NONZEROMAF <- filter(vcfs_ALL_SUMMARY_NONZEROMAF, status == "raw")
vcfs_ALL_SUMMARY_FILTERED_NONZEROMAF <- filter(vcfs_ALL_SUMMARY_NONZEROMAF, status == "filtered")

# Number of minor allele variable sites PER SAMPLE prefiltering
summary(vcfs_ALL_SUMMARY_RAW_NONZEROMAF$count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#7.0    25.0    65.0   134.3   190.0  1078.0 

summary(vcfs_ALL_SUMMARY_FILTERED_NONZEROMAF$count)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    3.00    8.00   22.41   39.50   90.00 



############################################################################################################
############################################################################################################
### SNPGENIE ON UNIX

# IN: /SARS-CoV-2-South-Africa/vcfs_filtered
# DO: snpgenie.pl --vcfformat=3 --fastafile=../NC_045512.2.fasta --gtffile=../WuhanHu1_v4.gtf 2>&1 | tee SNPGenie.out

### IMPORT DATA
codon_results <- read_tsv("/SARS-CoV-2-South-Africa/vcfs_filtered/SNPGenie_Results/codon_results.txt")
#View(codon_results)

# Make nsp12 one thing
codon_results[codon_results$product %in% c("nsp12_1", "nsp12_2"), ]$product <- "nsp12"

# Change nsp11 before site 13468 as non-overlapping (it's same frame as nsp12)
codon_results[codon_results$product == "nsp11" & codon_results$site < (13468 - 2), ]$num_overlap_ORF_nts <- 0
#View(filter(codon_results, product == "nsp11"))

# Rename product column
names(codon_results)[names(codon_results) == "product"] <- "product_segment"

# Rename sample and remove extension
names(codon_results)[names(codon_results) == "file"] <- "sample"
codon_results$sample <- str_replace(string = codon_results$sample, pattern = ".lf_filtered.vcf", replacement = "")

# Get unique product names
(uniq_product_segments <- unique(codon_results$product_segment))

# Number the codons for each product segment

# Obtain unique products and sites
(productSegment_uniqSites <- codon_results %>%
    group_by(product_segment, site) %>%
    summarise(
      count = n()
    )) # automatically sorts by product and site

# Make sure they're all 109
sum(productSegment_uniqSites$count != 109) # 0

# remove count column
productSegment_uniqSites <- dplyr::select(productSegment_uniqSites, -count)

# add codon_num column
productSegment_uniqSites$codon_num <- NA

# loop each product and site
for (this_product in unique(productSegment_uniqSites$product_segment)) {
  #this_product <- "ORF6"
  (num_uniq_sites <- nrow(productSegment_uniqSites[productSegment_uniqSites$product_segment == this_product, ]))
  productSegment_uniqSites[productSegment_uniqSites$product_segment == this_product, ]$codon_num <- 1:num_uniq_sites
}

# Examine
#View(productSegment_uniqSites)

# JOIN codon numbers
codon_results <- left_join(x = codon_results, y = productSegment_uniqSites, by = c("product_segment", "site"))
codon_results <- dplyr::select(codon_results, sample, product, product_segment, codon_num, everything())
#View(codon_results)

# check codon numbers
codon_results %>%
  group_by(product_segment) %>%
  summarise(
    highest_codon_num = max(codon_num)
  ) # confirmed

#product highest_codon_num
#product_segment highest_codon_num
#1 E                              76
#2 M                             223
#3 N                             420
#4 nsp1                          180
#5 nsp10                         139
#6 nsp11                          13
#7 nsp12                         932
#8 nsp13                         601
#9 nsp14                         527
#10 nsp15                         346
## … with 23 more rows

# Manually rename products to group different segments of same product together
codon_results$product <- codon_results$product_segment
codon_results[codon_results$product == "nsp1", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp2", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp3", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp4", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp5", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp6", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp7", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp8", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp9", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp10", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp11", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp12", ]$product <- "ORF1ab"
#codon_results[codon_results$product == "nsp12_1", ]$product <- "ORF1ab"
#codon_results[codon_results$product == "nsp12_2", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp13", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp14", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp15", ]$product <- "ORF1ab"
codon_results[codon_results$product == "nsp16", ]$product <- "ORF1ab"

# Get unique product grouped names
(uniq_products <- unique(codon_results$product))

# Number the codons for each product

# Obtain unique products and sites
(product_uniqSites <- codon_results %>%
  group_by(product, site) %>%
  summarise(
    count = n()
  )) # automatically sorts by product and site

# Make sure they're all 109
sum(product_uniqSites$count != 109)
product_uniqSites[product_uniqSites$count != 109, ]
#product  site count
#<chr>   <dbl> <int>
#1 ORF1ab  13442   218
#2 ORF1ab  13445   218
#3 ORF1ab  13448   218
#4 ORF1ab  13451   218
#5 ORF1ab  13454   218
#6 ORF1ab  13457   218
#7 ORF1ab  13460   218
#8 ORF1ab  13463   218
#9 ORF1ab  13466   218
# That's okay, because these are shared by both nsp11 and nsp12

# remove count column
product_uniqSites <- dplyr::select(product_uniqSites, -count)

# add codon_num_ORF column
product_uniqSites$codon_num_ORF <- NA

# loop each product and site
for (this_product in unique(product_uniqSites$product)) {
  (num_uniq_sites <- nrow(product_uniqSites[product_uniqSites$product == this_product, ]))
  product_uniqSites[product_uniqSites$product == this_product, ]$codon_num_ORF <- 1:num_uniq_sites
}

# Examine
#View(product_uniqSites)

# JOIN codon numbers
codon_results <- left_join(x = codon_results, y = product_uniqSites, by = c("product", "site"))
codon_results <- dplyr::select(codon_results, sample, product, product_segment, codon_num, codon_num_ORF, everything())
#View(codon_results)

# check codon numbers
codon_results %>%
  group_by(product) %>%
  summarise(
    highest_codon_num = max(codon_num_ORF)
  ) # confirmed


################################################################################
### piN/piS by SAMPLE (whole-genome), NOL regions only
(codon_results_NOL_bySample <- filter(codon_results, num_overlap_ORF_nts == 0) %>%
    group_by(sample) %>%
    summarise(
      N_diffs = sum(N_diffs),
      S_diffs = sum(S_diffs),
      N_sites = sum(N_sites),
      S_sites = sum(S_sites)
    ))

codon_results_NOL_bySample$piN <- (codon_results_NOL_bySample$N_diffs / codon_results_NOL_bySample$N_sites)
codon_results_NOL_bySample$piS <- (codon_results_NOL_bySample$S_diffs / codon_results_NOL_bySample$S_sites)
codon_results_NOL_bySample$piNpiS <- codon_results_NOL_bySample$piN / codon_results_NOL_bySample$piS

# test piN==piS
(codon_results_NOL_bySample_MEANPIN <- mean(codon_results_NOL_bySample$piN, na.rm = TRUE)) # 0.0001361455
(codon_results_NOL_bySample_MEANPIS <- mean(codon_results_NOL_bySample$piS, na.rm = TRUE)) # 0.000151611
(codon_results_NOL_bySample_MEANPINPIS <- codon_results_NOL_bySample_MEANPIN / codon_results_NOL_bySample_MEANPIS) # 0.8979925
(codon_results_NOL_bySample_WILCOXP <- wilcox.test(codon_results_NOL_bySample$piN, codon_results_NOL_bySample$piS, paired = TRUE)) # p-value = 0.1971


###############################################################################
### IMPORT OUTBREAK INFO
sample_metadata <- read_tsv("/SARS-CoV-2-South-Africa/samples_metadata_purged_CWN.tsv")
#View(sample_metadata)

names(sample_metadata)[names(sample_metadata) == "seqIds_VCF"] <- "sample"
codon_results_NOL_bySample <- left_join(x = codon_results_NOL_bySample, y = sample_metadata, by = "sample")
#View(codon_results_NOL_bySample)

# RESULTS BY OUTBREAK
(codon_results_NOL_bySample$pi <- codon_results_NOL_bySample$piN + codon_results_NOL_bySample$piS)
summary(codon_results_NOL_bySample$pi)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 1.435e-05 5.970e-05 2.878e-04 5.721e-04 1.510e-03 

# SAVE
#write_tsv(codon_results_NOL_bySample, "/SARS-CoV-2-South-Africa/codon_results_NOL_bySample.tsv")

# SUMMARIZE BY OUTBREAK
(codon_results_NOL_bySample_byOutbreak <- codon_results_NOL_bySample %>%
  group_by(outbreak) %>%
  summarise(
    count = n(),
    
    min_pi = min(pi, na.rm = TRUE),
    median_pi = median(pi, na.rm = TRUE),
    mean_pi = mean(pi, na.rm = TRUE),
    max_pi = max(pi, na.rm = TRUE),
    sd_pi = sd(pi, na.rm = TRUE),
    
    min_piN = min(piN, na.rm = TRUE),
    median_piN = median(piN, na.rm = TRUE),
    mean_piN = mean(piN, na.rm = TRUE),
    max_piN = max(piN, na.rm = TRUE),
    sd_piN = sd(piN, na.rm = TRUE),
    
    min_piS = min(piS, na.rm = TRUE),
    median_piS = median(piS, na.rm = TRUE),
    mean_piS = mean(piS, na.rm = TRUE),
    max_piS = max(piS, na.rm = TRUE),
    sd_piS = sd(piS, na.rm = TRUE)
  ))

codon_results_NOL_bySample_byOutbreak$mean_piNpiS <- codon_results_NOL_bySample_byOutbreak$mean_piN / codon_results_NOL_bySample_byOutbreak$mean_piS
codon_results_NOL_bySample_byOutbreak
#View(codon_results_NOL_bySample_byOutbreak)
#outbreak count min_pi median_pi  mean_pi  max_pi    sd_pi min_piN median_piN  mean_piN  max_piN   sd_piN min_piS median_piS  mean_piS  max_piS   sd_piS mean_piNpiS
#<chr>    <int>  <dbl>     <dbl>    <dbl>   <dbl>    <dbl>   <dbl>      <dbl>     <dbl>    <dbl>    <dbl>   <dbl>      <dbl>     <dbl>    <dbl>    <dbl>       <dbl>
#1 CH1         35      0 0.000572  0.000508 0.00151 0.000359       0  0.000269  0.000246  0.000654 0.000168       0 0.000289   0.000262  0.000856 0.000202       0.938
#2 CH3         74      0 0.0000350 0.000184 0.00131 0.000323       0  0.0000188 0.0000842 0.000573 0.000147       0 0.00000550 0.0000993 0.000749 0.000182       0.848

# pi, CH1 vs. CH3
# CH1 pi = 0.0005080586
# CH3 pi = 0.0001835596
wilcox.test(filter(codon_results_NOL_bySample, outbreak == "CH1")$pi,
            filter(codon_results_NOL_bySample, outbreak == "CH3")$pi,
            paired = FALSE)
#Wilcoxon rank sum test with continuity correction
#data:  filter(codon_results_NOL_bySample, outbreak == "CH1")$pi and filter(codon_results_NOL_bySample, outbreak == "CH3")$pi
#W = 1948.5, p-value = 2.211e-05
#alternative hypothesis: true location shift is not equal to 0

# piN, CH1 vs. CH3
# CH1 piN = 0.000246
# CH3 piN = 0.0000842
wilcox.test(filter(codon_results_NOL_bySample, outbreak == "CH1")$piN,
            filter(codon_results_NOL_bySample, outbreak == "CH3")$piN,
            paired = FALSE)
#Wilcoxon rank sum test with continuity correction
#data:  filter(codon_results_NOL_bySample, outbreak == "CH1")$piN and filter(codon_results_NOL_bySample, outbreak == "CH3")$piN
#W = 2001, p-value = 4.415e-06
#alternative hypothesis: true location shift is not equal to 0

# piS, CH1 vs. CH3
# CH1 piS = 0.000262
# CH3 piS = 0.0000994
wilcox.test(filter(codon_results_NOL_bySample, outbreak == "CH1")$piS,
            filter(codon_results_NOL_bySample, outbreak == "CH3")$piS,
            paired = FALSE)
#Wilcoxon rank sum test with continuity correction
#data:  filter(codon_results_NOL_bySample, outbreak == "CH1")$piS and filter(codon_results_NOL_bySample, outbreak == "CH3")$piS
#W = 1925, p-value = 2.787e-05
#alternative hypothesis: true location shift is not equal to 0

# BUT RATIOS EXTREMELY SIMILAR:
# CH1 piN/piS = 0.938
# CH3 piN/piS = 0.847

# CH1 piN vs. piS
# piN = 0.000246
# piS = 0.000262
wilcox.test(filter(codon_results_NOL_bySample, outbreak == "CH1")$piN,
            filter(codon_results_NOL_bySample, outbreak == "CH1")$piS,
            paired = TRUE)
#Wilcoxon signed rank test with continuity correction
#data:  filter(codon_results_NOL_bySample, outbreak == "CH1")$piN and filter(codon_results_NOL_bySample, outbreak == "CH1")$piS
#V = 213, p-value = 0.345

# CH3 piN vs. piS
# piN = 0.000246
# piS = 0.000262
wilcox.test(filter(codon_results_NOL_bySample, outbreak == "CH3")$piN,
            filter(codon_results_NOL_bySample, outbreak == "CH3")$piS,
            paired = TRUE)
#Wilcoxon signed rank test with continuity correction
#data:  filter(codon_results_NOL_bySample, outbreak == "CH3")$piN and filter(codon_results_NOL_bySample, outbreak == "CH3")$piS
#V = 885, p-value = 0.4017
#alternative hypothesis: true location shift is not equal to 0


### Let's make sure there isn't a laboratory artifact
wilcox.test(filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH")$piN,
            filter(codon_results_NOL_bySample, originating_lab != "NHLS-IALCH")$piN,
            paired = FALSE)
#Wilcoxon rank sum test with continuity correction
#data:  filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH")$piN and filter(codon_results_NOL_bySample, originating_lab != "NHLS-IALCH")$piN
#W = 643, p-value = 0.00023
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH")$piS,
            filter(codon_results_NOL_bySample, originating_lab != "NHLS-IALCH")$piS,
            paired = FALSE)
#Wilcoxon rank sum test with continuity correction
#data:  filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH")$piS and filter(codon_results_NOL_bySample, originating_lab != "NHLS-IALCH")$piS
#W = 738, p-value = 0.00189
#alternative hypothesis: true location shift is not equal to 0

# So we DO need to double check this

(codon_results_NOL_bySample_byLaboratory <- codon_results_NOL_bySample %>%
  group_by(originating_lab, outbreak) %>%
  summarise(
    count = n(),
    
    min_pi = min(pi, na.rm = TRUE),
    median_pi = median(pi, na.rm = TRUE),
    mean_pi = mean(pi, na.rm = TRUE),
    max_pi = max(pi, na.rm = TRUE),
    sd_pi = sd(pi, na.rm = TRUE),
    
    min_piN = min(piN, na.rm = TRUE),
    median_piN = median(piN, na.rm = TRUE),
    mean_piN = mean(piN, na.rm = TRUE),
    max_piN = max(piN, na.rm = TRUE),
    sd_piN = sd(piN, na.rm = TRUE),
    
    min_piS = min(piS, na.rm = TRUE),
    median_piS = median(piS, na.rm = TRUE),
    mean_piS = mean(piS, na.rm = TRUE),
    max_piS = max(piS, na.rm = TRUE),
    sd_piS = sd(piS, na.rm = TRUE)
  ))

#originating_lab      outbreak count   min_pi median_pi  mean_pi   max_pi   sd_pi min_piN median_piN mean_piN max_piN  sd_piN min_piS median_piS mean_piS max_piS  sd_piS
#<chr>                <chr>    <int>    <dbl>     <dbl>    <dbl>    <dbl>   <dbl>   <dbl>      <dbl>    <dbl>   <dbl>   <dbl>   <dbl>      <dbl>    <dbl>   <dbl>   <dbl>
#1 AMPATH-DBN           CH1         25 0        0.000572  0.000483 0.00151  4.06e-4 0.       0.000271   2.32e-4 6.54e-4 1.87e-4 0.      0.000288    2.51e-4 8.56e-4 2.26e-4
#2 Molecular Diagnosti… CH1          2 0.000230 0.000341  0.000341 0.000452 1.57e-4 1.88e-4  0.000195   1.95e-4 2.01e-4 9.29e-6 2.84e-5 0.000146    1.46e-4 2.64e-4 1.67e-4
#3 Molecular Diagnosti… CH1          3 0.000352 0.000449  0.000548 0.000844 2.61e-4 1.86e-4  0.000241   3.01e-4 4.76e-4 1.54e-4 1.66e-4 0.000208    2.47e-4 3.68e-4 1.07e-4
#4 NHLS-IALCH           CH1          5 0.000568 0.000669  0.000677 0.000829 1.03e-4 2.40e-4  0.000274   3.02e-4 4.69e-4 9.53e-5 3.27e-4 0.000360    3.75e-4 4.73e-4 6.04e-5
#5 NHLS-IALCH           CH3         74 0        0.0000350 0.000184 0.00131  3.23e-4 0.       0.0000188  8.42e-5 5.73e-4 1.47e-4 0.      0.00000550  9.93e-5 7.49e-4 1.82e-4

wilcox.test(filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH", outbreak == "CH1")$piN,
            filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH", outbreak == "CH3")$piN,
            paired = FALSE)
#Wilcoxon rank sum test with continuity correction
#data:  filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH",  and filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH",     outbreak == "CH1")$piN and     outbreak == "CH3")$piN
#W = 334, p-value = 0.002678
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH", outbreak == "CH1")$piS,
            filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH", outbreak == "CH3")$piS,
            paired = FALSE)
#Wilcoxon rank sum test with continuity correction
#data:  filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH",  and filter(codon_results_NOL_bySample, originating_lab == "NHLS-IALCH",     outbreak == "CH1")$piS and     outbreak == "CH3")$piS
#W = 331, p-value = 0.00226
#alternative hypothesis: true location shift is not equal to 0

# GREAT! So, even true for samples from the same laboratory.


###############################################################################
### SUMMARIZE RESULTS BY CODON MEANS, NOL ONLY
unique(codon_results$product_segment) # includes nsp11
(codon_results_NOL_byProductCodon <- filter(codon_results, num_overlap_ORF_nts == 0) %>%
    #group_by(product, codon_num) %>%
   group_by(product_segment, codon_num) %>%
    summarise(
      N_sites = mean(N_sites, na.rm = TRUE),
      S_sites = mean(S_sites, na.rm = TRUE),
      N_diffs = mean(N_diffs, na.rm = TRUE),
      S_diffs = mean(S_diffs, na.rm = TRUE)
    ))
unique(codon_results_NOL_byProductCodon$product_segment) # DOES include nsp11 because we correctly updated the OLG sites
names(codon_results_NOL_byProductCodon)[names(codon_results_NOL_byProductCodon) == "product_segment"] <- "product"
unique(codon_results_NOL_byProductCodon$product) # DOES include nsp11 because we correctly updated the OLG sites


###############################################################################
# BOOTSTRAP PROCESS

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
codon_results_NOL_byProductCodon$num_defined_seqs <- 6

### SUMMARIZE RESULTS BY GENE
(codon_results_NOL_byProductCodon_summary <- codon_results_NOL_byProductCodon %>%
    group_by(product) %>%
    summarise(
      N_sites = sum(N_sites, na.rm = TRUE),
      S_sites = sum(S_sites, na.rm = TRUE),
      N_diffs = sum(N_diffs, na.rm = TRUE),
      S_diffs = sum(S_diffs, na.rm = TRUE)
    )
)

codon_results_NOL_byProductCodon_summary$dN <- codon_results_NOL_byProductCodon_summary$N_diffs / codon_results_NOL_byProductCodon_summary$N_sites
codon_results_NOL_byProductCodon_summary$dS <- codon_results_NOL_byProductCodon_summary$S_diffs / codon_results_NOL_byProductCodon_summary$S_sites
codon_results_NOL_byProductCodon_summary$dNdS <- codon_results_NOL_byProductCodon_summary$dN / codon_results_NOL_byProductCodon_summary$dS
#View(codon_results_NOL_byProductCodon_summary)

(codon_results_NOL_byProductCodon_summary_LONG <- codon_results_NOL_byProductCodon_summary %>%
    pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value"))


####################################################################################################
### BOOTSTRAP EACH STATISTIC

############################################################################################################
# *BASIC* BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT
dNdS_diff_boot_fun <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
  
  # Function for dN
  dN_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    return(dN)
  }
  
  # Function for dN
  dS_function <- function(D, indices) {
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    return(dS)
  }
  
  # Function for dN - dS
  dN_m_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_m_dS <- dN - dS
    return(dN_m_dS)
  }
  
  # Function for dN/dS
  dN_over_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_over_dS <- dN / dS
    return(dN_over_dS)
  }
  
  # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
  
  (dN <- sum(as.vector(codon_results[ , paste0(numerator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(numerator, "_sites")]), na.rm = TRUE))
  (dS <- sum(as.vector(codon_results[ , paste0(denominator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(denominator, "_sites")]), na.rm = TRUE))
  (dNdS <- dN / dS)
  
  # Run the BOOTSTRAPS
  # boot dN
  (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
  (dN <- boot_dN$t0)
  (boot_dN_SE <- sd(boot_dN$t))
  
  # boot dS
  (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dS <- boot_dS$t0)
  (boot_dS_SE <- sd(boot_dS$t))
  
  # boot dN - dS
  (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_m_dS <- boot_dN_m_dS$t0)
  (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
  (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
  (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
  
  # boot dN/dS
  (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_over_dS <- boot_dN_over_dS$t0)
  (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
  (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
  (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
  
  ### NEW: ASL (acheived significance level)
  boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0) # 345
  boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0) # 0
  boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0) # 655
  ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
               boot_dN_m_dS_SE, boot_dN_m_dS_P, 
               boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
               sep = "\t"))
}


############################################################################################################
### ANALYSIS VARIABLES: same as before
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 10000
NCPUS <- 6


############################################################################################################
### INITIALIZE DATA FRAME: intrahost
intrahost_results_bootstrap <- data.frame(gene_name = character(),
                                          num_bootstraps = integer(),
                                          min_defined_codons = integer(),
                                          num_codons = integer(),
                                          N_sites = numeric(),
                                          S_sites = numeric(),
                                          N_diffs = numeric(),
                                          S_diffs = numeric(),
                                          num_replicates = integer(),
                                          dN = numeric(),
                                          dS = numeric(),
                                          dNdS = numeric(),
                                          dN_m_dS = numeric(),
                                          boot_dN_SE = numeric(),
                                          boot_dS_SE = numeric(),
                                          boot_dN_over_dS_SE = numeric(),
                                          boot_dN_over_dS_P = numeric(),
                                          boot_dN_m_dS_SE = numeric(),
                                          P_value = numeric(),
                                          boot_dN_gt_dS_count = integer(), 
                                          boot_dN_eq_dS_count = integer(), 
                                          boot_dN_lt_dS_count = integer(), 
                                          ASL_dN_gt_dS_P = numeric(), 
                                          ASL_dN_lt_dS_P = numeric())


### LOOP EACH DATA SUBSET (4 minutes with 6 CPUs)
for (this_gene in sort(unique(codon_results_NOL_byProductCodon$product))) {
  
    # Filter by gene, frame, and minimum number of defined codons
    this_data <- filter(codon_results_NOL_byProductCodon, product == this_gene, num_defined_seqs >= MIN_DEFINED_CODONS) 
    
    if(nrow(this_data) >= 2) {
      # LEADING SUMMARY COLUMNS:
      N_sites <- sum(this_data$N_sites, na.rm = T)
      S_sites <- sum(this_data$S_sites, na.rm = T)
      N_diffs <- sum(this_data$N_diffs, na.rm = T)
      S_diffs <- sum(this_data$S_diffs, na.rm = T)
      
      summary_data <- paste(nrow(this_data),
                            N_sites, S_sites, 
                            N_diffs, S_diffs, 
                            sep = "\t")
      
      # BOOTSTRAP THE ONE RATIO
      boot_dNdS <- dNdS_diff_boot_fun(this_data, 'N', 'S', NBOOTSTRAPS, NCPUS)
      
      # RECORD HEADER
      boot_vector_names <- c('num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                             'num_replicates', 
                             'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                             'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
      
      boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
      
      # Add names
      names(boot_dNdS_vector) <- boot_vector_names
      
      # Prepare additional rows
      intrahost_results_bootstrap_ADDITION <- data.frame(gene_name = this_gene,
                                                         num_bootstraps = NBOOTSTRAPS,
                                                         min_defined_codons = MIN_DEFINED_CODONS,
                                                         num_codons = as.integer(boot_dNdS_vector['num_codons']),
                                                         N_sites = as.numeric(boot_dNdS_vector['N_sites']),
                                                         S_sites = as.numeric(boot_dNdS_vector['S_sites']),
                                                         N_diffs = as.numeric(boot_dNdS_vector['N_diffs']),
                                                         S_diffs = as.numeric(boot_dNdS_vector['S_diffs']),
                                                         num_replicates = as.integer(boot_dNdS_vector['num_replicates']),
                                                         dN = as.numeric(boot_dNdS_vector['dN']),
                                                         dS = as.numeric(boot_dNdS_vector['dS']),
                                                         dNdS = as.numeric(boot_dNdS_vector['dNdS']),
                                                         dN_m_dS = as.numeric(boot_dNdS_vector['dN_m_dS']),
                                                         boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE']),
                                                         boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE']),
                                                         boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE']),
                                                         boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P']),
                                                         boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE']),
                                                         P_value = as.numeric(boot_dNdS_vector['P_value']),
                                                         boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count']),
                                                         boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count']),
                                                         boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count']),
                                                         ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P']),
                                                         ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P']))
      
      # Add the 4 new rows to results
      intrahost_results_bootstrap <- rbind(intrahost_results_bootstrap, intrahost_results_bootstrap_ADDITION)
    }
}
# warnings: NAs introduced by coercion (OKAY)

names(intrahost_results_bootstrap) <- c('gene_name',
                                        'num_bootstraps', 'min_defined_codons', 'num_codons', 
                                        'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                        'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                        'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                        'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')


#View(intrahost_results_bootstrap)

### Manual 2-sided ASL P-value
intrahost_results_bootstrap$P_ALS <- NA
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P < intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P < intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_gt_dS_P
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P > intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P > intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_lt_dS_P
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$P_ALS) & intrahost_results_bootstrap$P_ALS == 0, ]$P_ALS <- 1 / NBOOTSTRAPS

### FDR
intrahost_results_bootstrap$Q_ASL_BH <- p.adjust(intrahost_results_bootstrap$P_ALS, method = "BH")
intrahost_results_bootstrap$Q_Z_BH <- p.adjust(intrahost_results_bootstrap$P_value, method = "BH")

### SAVE
#write_tsv(intrahost_results_bootstrap, "/SARS-CoV-2-South-Africa/intrahost_results_bootstrap_products.tsv")

#View(intrahost_results_bootstrap)


###############################################################################
### SUMMARIZE RESULTS BY CODON MEANS, NOL ONLY, NOW BY OUTBREAK!

# JOIN outbreak metadata
codon_results <- left_join(x = codon_results, y = sample_metadata, by = "sample")

# SUMMARIZE
(codon_results_NOL_byOutbreakProductCodon <- filter(codon_results, num_overlap_ORF_nts == 0) %>%
    group_by(outbreak, product_segment, codon_num) %>%
    summarise(
      N_sites = mean(N_sites, na.rm = TRUE),
      S_sites = mean(S_sites, na.rm = TRUE),
      N_diffs = mean(N_diffs, na.rm = TRUE),
      S_diffs = mean(S_diffs, na.rm = TRUE)
    ))
names(codon_results_NOL_byOutbreakProductCodon)[names(codon_results_NOL_byOutbreakProductCodon) == "product_segment"] <- "product"


###############################################################################
# BOOTSTRAP PROCESS

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
codon_results_NOL_byOutbreakProductCodon$num_defined_seqs <- 6

### SUMMARIZE RESULTS BY OUTBREAK AND GENE
(codon_results_NOL_byOutbreakProductCodon_summary <- codon_results_NOL_byOutbreakProductCodon %>%
    group_by(outbreak, product) %>%
    summarise(
      N_sites = sum(N_sites, na.rm = TRUE),
      S_sites = sum(S_sites, na.rm = TRUE),
      N_diffs = sum(N_diffs, na.rm = TRUE),
      S_diffs = sum(S_diffs, na.rm = TRUE)
    )
)

codon_results_NOL_byOutbreakProductCodon_summary$dN <- codon_results_NOL_byOutbreakProductCodon_summary$N_diffs / codon_results_NOL_byOutbreakProductCodon_summary$N_sites
codon_results_NOL_byOutbreakProductCodon_summary$dS <- codon_results_NOL_byOutbreakProductCodon_summary$S_diffs / codon_results_NOL_byOutbreakProductCodon_summary$S_sites
codon_results_NOL_byOutbreakProductCodon_summary$dNdS <- codon_results_NOL_byOutbreakProductCodon_summary$dN / codon_results_NOL_byOutbreakProductCodon_summary$dS
#View(codon_results_NOL_byOutbreakProductCodon_summary)

(codon_results_NOL_byOutbreakProductCodon_summary_LONG <- codon_results_NOL_byOutbreakProductCodon_summary %>%
    pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value"))


####################################################################################################
### BOOTSTRAP EACH STATISTIC: *SAME* BOOTSTRAP FUNCTION (dN - dS) as above

############################################################################################################
### ANALYSIS VARIABLES: same as before
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 10000
NCPUS <- 6


############################################################################################################
### INITIALIZE DATA FRAME: intrahost
intrahost_byOutbreak_results_bootstrap <- data.frame(gene_name = character(),
                                          num_bootstraps = integer(),
                                          min_defined_codons = integer(),
                                          num_codons = integer(),
                                          N_sites = numeric(),
                                          S_sites = numeric(),
                                          N_diffs = numeric(),
                                          S_diffs = numeric(),
                                          num_replicates = integer(),
                                          dN = numeric(),
                                          dS = numeric(),
                                          dNdS = numeric(),
                                          dN_m_dS = numeric(),
                                          boot_dN_SE = numeric(),
                                          boot_dS_SE = numeric(),
                                          boot_dN_over_dS_SE = numeric(),
                                          boot_dN_over_dS_P = numeric(),
                                          boot_dN_m_dS_SE = numeric(),
                                          P_value = numeric(),
                                          boot_dN_gt_dS_count = integer(), 
                                          boot_dN_eq_dS_count = integer(), 
                                          boot_dN_lt_dS_count = integer(), 
                                          ASL_dN_gt_dS_P = numeric(), 
                                          ASL_dN_lt_dS_P = numeric(), 
                                          outbreak = character()) ## ADD OUTBREAK


### LOOP EACH DATA SUBSET (4 minutes with 6 CPUs)
for (this_outbreak in sort(unique(codon_results_NOL_byOutbreakProductCodon$outbreak))) {
  
  for (this_gene in sort(unique(codon_results_NOL_byOutbreakProductCodon$product))) {
    
    # Filter by gene, frame, and minimum number of defined codons
    this_data <- filter(codon_results_NOL_byOutbreakProductCodon, product == this_gene, outbreak == this_outbreak, num_defined_seqs >= MIN_DEFINED_CODONS) 
    
    if(nrow(this_data) >= 2) {
      # LEADING SUMMARY COLUMNS:
      N_sites <- sum(this_data$N_sites, na.rm = T)
      S_sites <- sum(this_data$S_sites, na.rm = T)
      N_diffs <- sum(this_data$N_diffs, na.rm = T)
      S_diffs <- sum(this_data$S_diffs, na.rm = T)
      
      summary_data <- paste(nrow(this_data),
                            N_sites, S_sites, 
                            N_diffs, S_diffs, 
                            sep = "\t")
      
      # BOOTSTRAP THE ONE RATIO
      boot_dNdS <- dNdS_diff_boot_fun(this_data, 'N', 'S', NBOOTSTRAPS, NCPUS)
      
      # RECORD HEADER
      boot_vector_names <- c('num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                             'num_replicates', 
                             'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                             'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
      
      boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
      
      # Add names
      names(boot_dNdS_vector) <- boot_vector_names
      
      # Prepare additional rows
      intrahost_byOutbreak_results_bootstrap_ADDITION <- data.frame(gene_name = this_gene,
                                                                    num_bootstraps = NBOOTSTRAPS,
                                                                    min_defined_codons = MIN_DEFINED_CODONS,
                                                                    num_codons = as.integer(boot_dNdS_vector['num_codons']),
                                                                    N_sites = as.numeric(boot_dNdS_vector['N_sites']),
                                                                    S_sites = as.numeric(boot_dNdS_vector['S_sites']),
                                                                    N_diffs = as.numeric(boot_dNdS_vector['N_diffs']),
                                                                    S_diffs = as.numeric(boot_dNdS_vector['S_diffs']),
                                                                    num_replicates = as.integer(boot_dNdS_vector['num_replicates']),
                                                                    dN = as.numeric(boot_dNdS_vector['dN']),
                                                                    dS = as.numeric(boot_dNdS_vector['dS']),
                                                                    dNdS = as.numeric(boot_dNdS_vector['dNdS']),
                                                                    dN_m_dS = as.numeric(boot_dNdS_vector['dN_m_dS']),
                                                                    boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE']),
                                                                    boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE']),
                                                                    boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE']),
                                                                    boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P']),
                                                                    boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE']),
                                                                    P_value = as.numeric(boot_dNdS_vector['P_value']),
                                                                    boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count']),
                                                                    boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count']),
                                                                    boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count']),
                                                                    ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P']),
                                                                    ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P']),
                                                                    outbreak = this_outbreak)
      
      # Add the 4 new rows to results
      intrahost_byOutbreak_results_bootstrap <- rbind(intrahost_byOutbreak_results_bootstrap, intrahost_byOutbreak_results_bootstrap_ADDITION)
    }
  }
} # end this outbreak
# warnings: NAs introduced by coercion (OKAY)

names(intrahost_byOutbreak_results_bootstrap) <- c('gene_name',
                                        'num_bootstraps', 'min_defined_codons', 'num_codons', 
                                        'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                        'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                        'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                        'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P',
                                        "outbreak") ## ADDED THIS


#View(intrahost_byOutbreak_results_bootstrap)

### Manual 2-sided ASL P-value
intrahost_byOutbreak_results_bootstrap$P_ALS <- NA
intrahost_byOutbreak_results_bootstrap[! is.na(intrahost_byOutbreak_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_byOutbreak_results_bootstrap$ASL_dN_gt_dS_P < intrahost_byOutbreak_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_byOutbreak_results_bootstrap[! is.na(intrahost_byOutbreak_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_byOutbreak_results_bootstrap$ASL_dN_gt_dS_P < intrahost_byOutbreak_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_gt_dS_P
intrahost_byOutbreak_results_bootstrap[! is.na(intrahost_byOutbreak_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_byOutbreak_results_bootstrap$ASL_dN_gt_dS_P > intrahost_byOutbreak_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_byOutbreak_results_bootstrap[! is.na(intrahost_byOutbreak_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_byOutbreak_results_bootstrap$ASL_dN_gt_dS_P > intrahost_byOutbreak_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_lt_dS_P
intrahost_byOutbreak_results_bootstrap[! is.na(intrahost_byOutbreak_results_bootstrap$P_ALS) & intrahost_byOutbreak_results_bootstrap$P_ALS == 0, ]$P_ALS <- 1 / NBOOTSTRAPS

### FDR
intrahost_byOutbreak_results_bootstrap$Q_ASL_BH <- p.adjust(intrahost_byOutbreak_results_bootstrap$P_ALS, method = "BH")
intrahost_byOutbreak_results_bootstrap$Q_Z_BH <- p.adjust(intrahost_byOutbreak_results_bootstrap$P_value, method = "BH")

### SAVE
#write_tsv(intrahost_byOutbreak_results_bootstrap, "/SARS-CoV-2-South-Africa/intrahost_byOutbreak_results_bootstrap_products.tsv")

#View(intrahost_byOutbreak_results_bootstrap)



####################################################################################################
# RELOAD ALL
intrahost_results_bootstrap <- read_tsv("/SARS-CoV-2-South-Africa/intrahost_results_bootstrap_products.tsv")
intrahost_results_bootstrap$outbreak <- "Combined" # add this now
intrahost_byOutbreak_results_bootstrap <- read_tsv("/SARS-CoV-2-South-Africa/intrahost_byOutbreak_results_bootstrap_products.tsv")

# COMBINE
intrahost_results_bootstrap <- rbind(intrahost_results_bootstrap, intrahost_byOutbreak_results_bootstrap)

### CONVERT TO LONG
intrahost_results_bootstrap_LONG <- intrahost_results_bootstrap %>%
  pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value")
# so the 2 measures of d (dN and dS) are on different lines

#View(intrahost_results_bootstrap_LONG)


####################################################################################################
### PLOT WHOLE GENES, VARIOUS RATIOS
intrahost_results_bootstrap_LONG$d_measure <- factor(intrahost_results_bootstrap_LONG$d_measure, levels = c('dN', 'dS'))

# Add error bars
intrahost_results_bootstrap_LONG$d_SE_min <- NA
intrahost_results_bootstrap_LONG$d_SE_max <- NA
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_min <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_value - intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_max <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_value + intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_min <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_value - intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_max <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_value + intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_SE_min < 0 & ! is.na(intrahost_results_bootstrap_LONG$d_SE_min), ]$d_SE_min <- 0 # OK if fail because none negative

# significance for Q-value
intrahost_results_bootstrap_LONG$significance <- NA
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$Q_Z_BH) & intrahost_results_bootstrap_LONG$Q_Z_BH < 0.05, ]$significance <- '*'
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$Q_Z_BH) & intrahost_results_bootstrap_LONG$Q_Z_BH < 0.01, ]$significance <- '**'
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$Q_Z_BH) & intrahost_results_bootstrap_LONG$Q_Z_BH < 0.001, ]$significance <- '***'

# ORDER AND NAME THE GENES
unique(intrahost_results_bootstrap_LONG$gene_name)
unique(codon_results$product_segment)

#gene_ids_sorted <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
gene_ids_sorted <- c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp11", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", 
                     "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")

#gene_names_sorted <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
gene_names_sorted <- c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp11", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", 
                       "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")

intrahost_results_bootstrap_LONG$gene_name <- factor(intrahost_results_bootstrap_LONG$gene_name, 
                                                     levels = gene_ids_sorted, 
                                                     labels = gene_names_sorted)

intrahost_results_bootstrap_LONG$outbreak <- factor(intrahost_results_bootstrap_LONG$outbreak,
                                                    levels = c("Combined", "CH1", "CH3"))

# For horizontal lines
horizontal_lines_data <- filter(intrahost_results_bootstrap_LONG, significance %in% c('*', '**', '***')) %>%
  group_by(gene_name) %>% # frame
  summarise(
    y_value = max(d_SE_max)
  )

intrahost_error_bar_colors <- dplyr::select(intrahost_results_bootstrap_LONG, outbreak, gene_name, d_measure, d_value, d_SE_min, d_SE_max) # <-- CHANGE AS NEEDED
intrahost_error_bar_colors$this_color <- NA
intrahost_error_bar_colors[intrahost_error_bar_colors$d_measure == 'dN', ]$this_color <- 'pink'
intrahost_error_bar_colors[intrahost_error_bar_colors$d_measure == 'dS', ]$this_color <- 'lightblue'
intrahost_error_bar_colors$d_measure <- factor(intrahost_error_bar_colors$d_measure, levels = c('dN', 'dS'))

# Boxes for all based on which π is higher


### PREPARE PLOT

### normalized dNdS based on difference
intrahost_results_bootstrap_LONG$dNdS_norm <- NA
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$dNdS), ]$dNdS_norm <- 
  intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$dNdS), ]$dN_m_dS / 
  max(abs(intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$dNdS), ]$dN_m_dS))

### PLOT diversity chart ###
gene_names_allLevels <- gene_names_sorted # old: c('1ab', 'S', '3a', "3c", 'E', 'M', '6', '7a', '7b', '8', 'N', '9b', '9c', '10') 
pi_corr_factor <- 10000
max_dNdS <- max(intrahost_results_bootstrap_LONG$dNdS, na.rm = TRUE)
max_d_SE_max <- max(intrahost_results_bootstrap_LONG$d_SE_max, na.rm = TRUE)


### PLOT
(intrahost_allProductsByOutbreak_PLOT <- ggplot(data = filter(intrahost_results_bootstrap_LONG, gene_name %in% gene_names_allLevels),
                                    mapping = aes(x = gene_name, y = d_value * pi_corr_factor, group = d_measure)) + #fill = d_measure)) +
    # Backdrop boxes based on which pi is higher
    geom_bar(data = filter(intrahost_results_bootstrap_LONG, d_measure == 'dN', gene_name %in% gene_names_allLevels), mapping = aes(y = Inf, fill = dNdS_norm), stat = 'identity') +
    
    geom_errorbar(data = filter(intrahost_error_bar_colors, gene_name %in% gene_names_allLevels), 
                  mapping =  aes(ymin = d_SE_min * pi_corr_factor, ymax = d_SE_max * pi_corr_factor), # , color = d_measure
                  color = rev(rep(c('pink', 'lightblue'), 78)),
                  position = position_dodge(width = 0.5), width = 0, size = 1) +
    geom_point(data = filter(intrahost_error_bar_colors, gene_name %in% gene_names_allLevels), 
               mapping =  aes(y = d_SE_min * pi_corr_factor), # , color = d_measure
               color = rev(rep(c('pink', 'lightblue'), 78)),
               position = position_dodge(width = 0.5), size = 0.29) + # size = 1
    geom_point(data = filter(intrahost_error_bar_colors, gene_name %in% gene_names_allLevels), 
               mapping =  aes(y = d_SE_max * pi_corr_factor), # color = d_measure
               color = rev(rep(c('pink', 'lightblue'), 78)),
               position = position_dodge(width = 0.5), size = 0.29) +
    geom_point(stat = 'identity', position = position_dodge(width = 0.5), pch = 21, size = 2, stroke = 0, 
               fill = rev(rep(c(brewer.pal(9, 'Set1')[1], brewer.pal(9, 'Set1')[2]), 78))) + # aes(color = d_measure), 
    
    facet_grid(outbreak ~ ., scales = 'free', space = 'free_x',) +
    
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none',
          legend.title = element_blank(),
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), # angle = 90, hjust = 1, vjust = 1
          #axis.text.x = element_text(size = 9, face = "italic"), # angle = 90, hjust = 1, vjust = 1
          axis.text.y = element_text(size = 8),
          #axis.ticks.y.right = element_line(colour = '#DCB118'),
          #axis.text.y.right = element_text(color = '#DCB118'),
          axis.title.y = element_text(size = 8),
          panel.border = element_rect(),
          strip.text = element_text(size = 8),
          strip.text.y = element_blank(),
          #legend.background = element_rect(fill = 'white', color = 'black'),
          strip.background = element_blank()) +
    xlab("") + 
    ylab(bquote('Differences per site ('*'x 10'^'-4'*')')) + 
    scale_y_continuous(breaks = scales::pretty_breaks(3), expand = expand_scale(mult = c(0, 0.1))) + # 0.1
    scale_fill_gradient2(low = brewer.pal(9, "Blues")[5], mid = 'white', high = brewer.pal(9, "Reds")[5], midpoint = 0, na.value = "white")) 

# SAVE PLOT
#png(filename = "/SARS-CoV-2-South-Africa/intrahost_allProductsByOutbreak_PLOT.png", width = 6, height = 3.5, units = 'in', res = 500) # height = 2
intrahost_allProductsByOutbreak_PLOT
#dev.off()

# SAVE SOURCE
intrahost_results_bootstrap_LONG$significance_P <- NA
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$P_value) & intrahost_results_bootstrap_LONG$P_value < 0.1, ]$significance_P <- '*'
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$P_value) & intrahost_results_bootstrap_LONG$P_value < 0.01, ]$significance_P <- '**'
intrahost_results_bootstrap_LONG[! is.na(intrahost_results_bootstrap_LONG$P_value) & intrahost_results_bootstrap_LONG$P_value < 0.001, ]$significance_P <- '***'
#write_tsv(intrahost_results_bootstrap_LONG, "/SARS-CoV-2-South-Africa/intrahost_allProductsByOutbreak_PLOT_SOURCE.tsv")

View(dplyr::select(intrahost_results_bootstrap_LONG, outbreak, gene_name, dNdS, P_value, significance_P, Q_Z_BH, significance))

# PIVOT WIDE AND SAVE
intrahost_results_bootstrap_LONG_WIDE <- pivot_wider(data = dplyr::select(intrahost_results_bootstrap_LONG, -significance, -dNdS_norm, -significance_P), 
                                                                names_from = d_measure, values_from = c(d_value, d_SE_min, d_SE_max))
#write_tsv(intrahost_results_bootstrap_LONG_WIDE, "/SARS-CoV-2-South-Africa/intrahost_results_bootstrap_LONG_WIDE.tsv")
View(intrahost_results_bootstrap_LONG_WIDE)

# Both ratios above/below 1?
(ratio_direction_by_outbreak <- filter(intrahost_results_bootstrap, outbreak != "Combined") %>%
  group_by(gene_name) %>%
  summarise(
    pos = sum(dNdS > 1.3),
    neg = sum(dNdS < 0.7)
  ))
ratio_direction_by_outbreak$diff <- FALSE
ratio_direction_by_outbreak[! is.na(ratio_direction_by_outbreak$pos) & ratio_direction_by_outbreak$pos > 0 & 
                              ! is.na(ratio_direction_by_outbreak$neg) & ratio_direction_by_outbreak$neg > 0, ]$diff <- TRUE
#View(ratio_direction_by_outbreak)

genes_diff_ratios <- c("N", "nsp1", "nsp5", "nsp10", "nsp14", "ORF7a", "ORF10")


###############################################################################
### SLIDING WINDOWS

# SAVE CODON-MEAN RESULTS for each PRODUCT; INCLUDE OL and NOL codons

###############################################################################
### Per SAMPLE per GENE: SUMMARIZE RESULTS BY CODON MEANS
(codon_results_byProductCodon <- filter(codon_results) %>%
   group_by(product_segment, codon_num) %>%
   summarise(
     N_sites = mean(N_sites, na.rm = TRUE),
     S_sites = mean(S_sites, na.rm = TRUE),
     N_diffs = mean(N_diffs, na.rm = TRUE),
     S_diffs = mean(S_diffs, na.rm = TRUE)
   ))

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
codon_results_byProductCodon$num_defined_seqs <- 6

# LOOP AND WRITE
#for (this_product in unique(codon_results_byProductCodon$product_segment)) {
#  write_tsv(filter(codon_results_byProductCodon, product_segment == this_product), 
#            paste0("/SARS-CoV-2-South-Africa/intrahost_sliding_windows/codon_results_byProductCodon_", this_product, ".tsv"))
#  
#}

# custom make nsp12 (RdRp), which has a frameshift and is named as nsp12_1 nsp12_2
#write_tsv(
#  rbind(
#    filter(codon_results_byProductCodon, product_segment == "nsp12_1"), 
#    filter(codon_results_byProductCodon, product_segment == "nsp12_2")
#  ),
#  paste0("/SARS-CoV-2-South-Africa/intrahost_sliding_windows/codon_results_byProductCodon_nsp12.tsv")
#) # already numbered perfectly well!

# Write SNPGenie sliding window script for R, akin to what's used for OLGenie: SNPGenie_sliding_windows.R

# Manually prepare command line arguments in text editor:
# IN: /SARS-CoV-2-South-Africa/intrahost_sliding_windows
# DO: ls

# FIND: codon_results_byProductCodon_(\w+)\.tsv
# REPLACE: SNPGenie_sliding_windows.R codon_results_byProductCodon_\1.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_\1.out

# IN: /SARS-CoV-2-South-Africa/intrahost_sliding_windows
# DO: [manually make sure just one for nsp12]

# size=10
#SNPGenie_sliding_windows.R codon_results_byProductCodon_E.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_E.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF3d.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3d.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_S.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_S.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp4.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp4.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_M.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_M.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF6.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF6.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_SiORF1.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF1.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp13.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp13.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp5.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp5.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_N.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_N.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF7a.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7a.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_SiORF2.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF2.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp14.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp14.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp6.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp6.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF10.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF10.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF7b.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7b.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp1.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp1.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp15.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp15.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp7.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp7.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF3a.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3a.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF8.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF8.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp10.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp10.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp16.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp16.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp8.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp8.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF3b.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3b.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF9b.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9b.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp11.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp11.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp2.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp2.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp9.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp9.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF3c.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3c.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_ORF9c.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9c.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp3.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp3.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp12.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp12.out

#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp12_1.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp12_1.out
#SNPGenie_sliding_windows.R codon_results_byProductCodon_nsp12_2.tsv N S 10 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp12_2.out


# size=20
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_E.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_E.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3d.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3d.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_S.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_S.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp4.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp4.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_M.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_M.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF6.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF6.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_SiORF1.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF1.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp13.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp13.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp5.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp5.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_N.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_N.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF7a.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7a.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_SiORF2.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF2.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp14.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp14.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp6.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp6.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF10.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF10.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF7b.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp1.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp1.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp15.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp15.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp7.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp7.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3a.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3a.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF8.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF8.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp10.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp10.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp16.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp16.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp8.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp8.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3b.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF9b.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp11.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp11.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp2.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp2.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp9.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp9.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3c.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3c.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF9c.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9c.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp3.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp3.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp12.tsv N S 20 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp12.out


# size=30
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_E.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_E.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3d.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3d.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_S.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_S.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp4.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp4.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_M.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_M.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF6.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF6.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_SiORF1.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF1.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp13.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp13.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp5.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp5.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_N.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_N.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF7a.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7a.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_SiORF2.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF2.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp14.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp14.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp6.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp6.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF10.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF10.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF7b.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp1.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp1.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp15.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp15.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp7.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp7.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3a.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3a.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF8.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF8.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp10.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp10.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp16.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp16.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp8.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp8.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3b.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF9b.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp11.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp11.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp2.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp2.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp9.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp9.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3c.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3c.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF9c.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9c.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp3.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp3.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp12.tsv N S 30 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp12.out

# size=40
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_E.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_E.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3d.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3d.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_S.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_S.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp4.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp4.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_M.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_M.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF6.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF6.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_SiORF1.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF1.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp13.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp13.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp5.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp5.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_N.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_N.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF7a.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7a.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_SiORF2.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF2.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp14.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp14.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp6.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp6.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF10.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF10.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF7b.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp1.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp1.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp15.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp15.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp7.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp7.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3a.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3a.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF8.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF8.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp10.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp10.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp16.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp16.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp8.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp8.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3b.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF9b.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp11.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp11.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp2.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp2.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp9.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp9.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3c.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3c.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF9c.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9c.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp3.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp3.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp12.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp12.out

# size=50
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_E.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_E.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3d.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3d.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_S.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_S.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp4.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp4.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_M.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_M.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF6.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF6.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_SiORF1.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF1.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp13.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp13.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp5.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp5.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_N.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_N.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF7a.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7a.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_SiORF2.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_SiORF2.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp14.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp14.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp6.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp6.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF10.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF10.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF7b.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF7b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp1.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp1.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp15.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp15.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp7.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp7.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3a.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3a.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF8.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF8.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp10.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp10.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp16.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp16.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp8.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp8.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3b.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF9b.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9b.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp11.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp11.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp2.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp2.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp9.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp9.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF3c.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF3c.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_ORF9c.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_ORF9c.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp3.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp3.out
#SNPGenie_sliding_windows.R ../codon_results_byProductCodon_nsp12.tsv N S 50 1 1000 6 NONE 6 > SNPGenie_sliding_windows_nsp12.out


###############################################################################
###############################################################################
### PLOT sliding windows
# Combine files at command line:

# IN: /SARS-CoV-2-South-Africa/intrahost_sliding_windows_size30
# DO: cat codon_results_byProductCodon_*.tsv | grep -v product_segment > codon_results_byProductCodon_WINDOWS_POOLED.tsv

# CHANGE THIS
#WINDOW_SIZE <- 10
#WINDOW_SIZE <- 20
WINDOW_SIZE <- 30
#WINDOW_SIZE <- 40
#WINDOW_SIZE <- 50

# CHANGE THIS
codon_results_byProductCodon_WINDOWS_POOLED <- 
  #read_tsv("/SARS-CoV-2-South-Africa/intrahost_sliding_windows_size10/codon_results_byProductCodon_WINDOWS_POOLED.tsv", 
  #read_tsv("/SARS-CoV-2-South-Africa/intrahost_sliding_windows_size20/codon_results_byProductCodon_WINDOWS_POOLED.tsv", 
  read_tsv("/SARS-CoV-2-South-Africa/intrahost_sliding_windows_size30/codon_results_byProductCodon_WINDOWS_POOLED.tsv", 
  #read_tsv("/SARS-CoV-2-South-Africa/intrahost_sliding_windows_size40/codon_results_byProductCodon_WINDOWS_POOLED.tsv", 
  #read_tsv("/SARS-CoV-2-South-Africa/intrahost_sliding_windows_size50/codon_results_byProductCodon_WINDOWS_POOLED.tsv", 
           col_names = c("product", "codon_num", "N_sites", "S_sites", "N_diffs", "S_diffs", "num_defined_seqs", "sw_ratio", 
                         "sw_start", "sw_center", "sw_end", "sw_num_replicates", "sw_N_diffs", "sw_S_diffs", "sw_N_sites", "sw_S_sites", 
                         "sw_dN", "sw_dS", "sw_dNdS", "sw_dN_m_dS", "sw_boot_dN_SE", "sw_boot_dS_SE", "sw_boot_dN_over_dS_SE", "sw_boot_dN_over_dS_P", 
                         "sw_boot_dN_m_dS_SE", "sw_boot_dN_m_dS_P", "sw_boot_dN_gt_dS_count", "sw_boot_dN_eq_dS_count", "sw_boot_dN_lt_dS_count", 
                         "sw_ASL_dN_gt_dS_P", "sw_ASL_dN_lt_dS_P", "sw_ASL_dNdS_P"))

#View(codon_results_byProductCodon_WINDOWS_POOLED)

# Renumber codons starting at 1 for each product
#codon_results_byProductCodon_WINDOWS_POOLED$codon_num_ORF <- codon_results_byProductCodon_WINDOWS_POOLED$codon_num
codon_results_byProductCodon_WINDOWS_POOLED$sw_center_ORF <- codon_results_byProductCodon_WINDOWS_POOLED$sw_center
for (this_product in unique(codon_results_byProductCodon_WINDOWS_POOLED$product)) {
  #this_product <- "nsp11"
  
  # codon
  codon_results_byProductCodon_WINDOWS_POOLED[codon_results_byProductCodon_WINDOWS_POOLED$product == this_product, ]$codon_num <- 
    codon_results_byProductCodon_WINDOWS_POOLED[codon_results_byProductCodon_WINDOWS_POOLED$product == this_product, ]$codon_num - 
    min(codon_results_byProductCodon_WINDOWS_POOLED[codon_results_byProductCodon_WINDOWS_POOLED$product == this_product, ]$codon_num) + 1
}

# Renumber sw_center
codon_results_byProductCodon_WINDOWS_POOLED$sw_center <- codon_results_byProductCodon_WINDOWS_POOLED$codon_num + ((WINDOW_SIZE - 1) / 2)
codon_results_byProductCodon_WINDOWS_POOLED[is.na(codon_results_byProductCodon_WINDOWS_POOLED$sw_start), ]$sw_center <- NA
#View(codon_results_byProductCodon_WINDOWS_POOLED)

# Rename nsp12
codon_results_byProductCodon_WINDOWS_POOLED[codon_results_byProductCodon_WINDOWS_POOLED$product %in% c("nsp12_1", "nsp12_2"), ]$product <- "nsp12"

# Make a long form for plotting dN and dS
(codon_results_byProductCodon_WINDOWS_POOLED_LONG <- codon_results_byProductCodon_WINDOWS_POOLED %>%
    pivot_longer(cols = c('sw_dN', 'sw_dS'), names_to = "d_measure", values_to = "d_value"))

### DEFINE products
(uniq_products_LONG <- unique(codon_results_byProductCodon_WINDOWS_POOLED_LONG$product)) # nsp not included because only 8 codons!
# "E", "M", "N", "ORF10", "ORF3a", "ORF3c", "ORF3d", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "ORF9c", "S", "SiORF1", "SiORF2", "nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", "nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", 
gene_names_sorted # recall
#"ORF1ab" "S"      "ORF3a"  "E"      "M"      "ORF6"   "ORF7a"  "ORF7b"  "ORF8"   "N"      "ORF10" 

uniq_products_LONG_SORTED <- c(
  "nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16",
  "S", "SiORF1", "SiORF2", "ORF3a", "ORF3c", "ORF3d", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF9c", "ORF10"
)

uniq_products_LONG_SORTED_KEEPERS <- c(
  "nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16",
  "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"
)

# FACTOR products
codon_results_byProductCodon_WINDOWS_POOLED_LONG$product <- factor(codon_results_byProductCodon_WINDOWS_POOLED_LONG$product, 
                                                                   levels = uniq_products_LONG_SORTED)

###############################################################################
### SUMMARIZE piN and piS for each product (NOT LONG)
(codon_results_byProductCodon_byProduct <- filter(codon_results_byProductCodon_WINDOWS_POOLED) %>%
   group_by(product) %>%
   summarise(
     N_sites = sum(N_sites, na.rm = TRUE),
     S_sites = sum(S_sites, na.rm = TRUE),
     N_diffs = sum(N_diffs, na.rm = TRUE),
     S_diffs = sum(S_diffs, na.rm = TRUE)
   ))

#View(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10"))
#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10")$S_diffs) # 0.008516171
#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10")$S_sites) # 26.33383
#0.008516171/26.33383 # 0.0003233928

#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10", codon_num >= 1, codon_num <= 30)$S_diffs) # 0.008516171
#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10", codon_num >= 1, codon_num <= 30)$S_sites) # 21.0009
#0.008516171/21.0009 # 0.0004055146

#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10", codon_num >= 10, codon_num <= 39)$S_diffs) # 0.008516171
#sum(filter(codon_results_byProductCodon_WINDOWS_POOLED, product == "ORF10", codon_num >= 10, codon_num <= 39)$S_sites) # 20.66626
#0.008516171/20.66626 # 0.0004120809

# UNDERSTAND IT NOW! EXPLANATION:
# I was confused that the mean piS for all of ORF10 could be LOWER than piS for EVERY WINDOW. However, the only S differences are observed in
# the middle of the gene. Thus, because the window size is 30 codons, every window overlaps the only S differences that exist. As a result,
# the only effect of calculating the total piS for the whole gene is to add S sites but now S differences, and the denominator but not
# numerator grows. This is why the total value of piS for the gene (horitzontal dashed line) can be lower than that of any given window.

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
codon_results_byProductCodon_byProduct$num_defined_seqs <- 6

# ADD STATS
codon_results_byProductCodon_byProduct$dN <- codon_results_byProductCodon_byProduct$N_diffs / codon_results_byProductCodon_byProduct$N_sites
codon_results_byProductCodon_byProduct$dS <- codon_results_byProductCodon_byProduct$S_diffs / codon_results_byProductCodon_byProduct$S_sites
codon_results_byProductCodon_byProduct$dNdS <- codon_results_byProductCodon_byProduct$dN / codon_results_byProductCodon_byProduct$dS
#View(codon_results_byProductCodon_byProduct)

# FACTOR MEANS
codon_results_byProductCodon_byProduct$product <- factor(codon_results_byProductCodon_byProduct$product, 
                                                         levels = uniq_products_LONG_SORTED)

# Create STANDARD ERROR COLUMN
#View(codon_results_byProductCodon_WINDOWS_POOLED_LONG)
codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_SE <- codon_results_byProductCodon_WINDOWS_POOLED_LONG$sw_boot_dN_SE
codon_results_byProductCodon_WINDOWS_POOLED_LONG[codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_measure == "sw_dS", ]$d_value_SE <- 
  codon_results_byProductCodon_WINDOWS_POOLED_LONG[codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_measure == "sw_dS", ]$sw_boot_dS_SE
codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_CI_min <- 
  codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value - codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_SE
codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_CI_max <- 
  codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value + codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_SE

# reset negative mins to 0
codon_results_byProductCodon_WINDOWS_POOLED_LONG[!is.na(codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_CI_min) & 
                                                   codon_results_byProductCodon_WINDOWS_POOLED_LONG$d_value_CI_min < 0, ]$d_value_CI_min <- 0

# Overlapping genes OLGs to mask
OLGs_to_mask <- data.frame(
  product = c("S", "ORF3a", "ORF3a", "ORF3a", "N", "N"), # "nsp12", 
  name = c("S-iORF", "ORF3c", "ORF3d", "ORF3b", "ORF9b", "ORF9c"), # "nsp11", 
  start_codon = c(61, 22+(WINDOW_SIZE-22), 44, 141, 4+(WINDOW_SIZE-4), 154), # 8, 
  end_codon = c(101, 64, 102, 164, 102, 228), # 14, 
  sw_center = rep(0, 6), # 
  d_value = rep(0, 6), # 
  d_measure = rep("sw_dN", 6) # 
)

# PLOT
pi_corr_factor2 <- 1000
(codon_results_byProductCodon_WINDOWS_POOLED_LONG_PLOT <- ggplot(data = filter(codon_results_byProductCodon_WINDOWS_POOLED_LONG, product %in% uniq_products_LONG_SORTED_KEEPERS), # product == "ORF3a"), 
                                                                 mapping = aes(x = sw_center, y = pi_corr_factor2 * d_value, color = d_measure)) + 
    
    facet_wrap(product ~ ., scales = "free") +
    geom_rect(data = OLGs_to_mask, mapping = aes(xmin = start_codon, xmax = end_codon, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "grey", color = NA, 
              size = .3) + 
    geom_line(size = .3) +
    #geom_vline(xintercept = 614) + # to see S 614
    geom_hline(data = filter(codon_results_byProductCodon_byProduct, product %in% uniq_products_LONG_SORTED_KEEPERS),
               mapping = aes(yintercept = pi_corr_factor2 * dS), linetype = "dashed", color = brewer.pal(9, 'Set1')[2], 
               size = .3) +
    geom_ribbon(mapping = aes(ymin = pi_corr_factor2 * d_value_CI_min, ymax = pi_corr_factor2 * d_value_CI_max, fill = d_measure), alpha = 0.25, linetype = 0, 
                size = .3) +
    theme_bw() +
    theme(strip.background = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          axis.ticks = element_line(size = 0.3),
          legend.title = element_blank(),
          legend.position = 'none',
          axis.text = element_text(size = 8),
          strip.text = element_text(face = c("bold")),
          panel.grid = element_blank()) +
    xlab("Codon") + 
    ylab(bquote('Nucleotide diversity ('*'x 10'^'-3'*')')) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3), expand = c(0, 0)) + 
    scale_color_manual(values = c(brewer.pal(9, 'Set1')[1], brewer.pal(9, 'Set1')[2])) 
)


# SAVE FIGURE SLIDING WINDOW
png(
  #filename = "/SARS-CoV-2-South-Africa/codon_results_byProductCodon_WINDOWS_POOLED_LONG_size10_PLOT.png", 
  #filename = "/SARS-CoV-2-South-Africa/codon_results_byProductCodon_WINDOWS_POOLED_LONG_size20_PLOT.png", 
  filename = "/SARS-CoV-2-South-Africa/codon_results_byProductCodon_WINDOWS_POOLED_LONG_size30_PLOT.png", 
  #filename = "/SARS-CoV-2-South-Africa/codon_results_byProductCodon_WINDOWS_POOLED_LONG_size40_PLOT.png", 
  #filename = "/SARS-CoV-2-South-Africa/codon_results_byProductCodon_WINDOWS_POOLED_LONG_size50_PLOT.png", 
    width = 8, height = 6, units = 'in', res = 500)
codon_results_byProductCodon_WINDOWS_POOLED_LONG_PLOT
dev.off()


### CANDIDATE WINDOWS for 30 codons

# PIVOT to SHORT again
#View(codon_results_byProductCodon_WINDOWS_POOLED_LONG)
nrow(codon_results_byProductCodon_WINDOWS_POOLED_LONG) # 20198
codon_results_byProductCodon_WINDOWS_POOLED_WIDE <- pivot_wider(data = codon_results_byProductCodon_WINDOWS_POOLED_LONG, 
                                                                 names_from = d_measure, values_from = c(d_value, d_value_SE, d_value_CI_min, d_value_CI_max))
#View(codon_results_byProductCodon_WINDOWS_POOLED_WIDE)
nrow(codon_results_byProductCodon_WINDOWS_POOLED_WIDE) # 10099. 10099*2 = 20198, CONFIRMED!

# JOIN the PRODUCT'S values
product_dS_values <- dplyr::select(codon_results_byProductCodon_byProduct, product, dS)
names(product_dS_values) <- c("product", "product_dS")
codon_results_byProductCodon_WINDOWS_POOLED_WIDE <- left_join(x = codon_results_byProductCodon_WINDOWS_POOLED_WIDE, product_dS_values, by = "product")
nrow(codon_results_byProductCodon_WINDOWS_POOLED_WIDE) # 10099

# Get the candidates
CANDIDATE_WINDOWS_SIZE30 <- filter(codon_results_byProductCodon_WINDOWS_POOLED_WIDE, 
                                   d_value_CI_min_sw_dN > d_value_CI_max_sw_dS, 
                                   d_value_CI_min_sw_dN > product_dS)
#View(CANDIDATE_WINDOWS_SIZE30)

# SAVE CANDIDATES
#write_tsv(CANDIDATE_WINDOWS_SIZE30, "/SARS-CoV-2-South-Africa/CANDIDATE_WINDOWS_SIZE30.tsv")


################################################################################
# MANUALLY DEFINED START TO END OF CANDIDATE WINDOWS IN CANDIDATE_WINDOWS_SIZE30.xlsx, saved in CANDIDATE_WINDOWS_SIZE30_table.tsv
# LOOP and define/test each window; list codons with nonsyn diffs; etc.
CANDIDATE_WINDOWS_SIZE30_table <- read_tsv("/SARS-CoV-2-South-Africa/CANDIDATE_WINDOWS_SIZE30_table.tsv")
CANDIDATE_WINDOWS_SIZE30_table$N_diffs <- NA
CANDIDATE_WINDOWS_SIZE30_table$S_diffs <- NA
CANDIDATE_WINDOWS_SIZE30_table$N_sites <- NA
CANDIDATE_WINDOWS_SIZE30_table$S_sites <- NA
CANDIDATE_WINDOWS_SIZE30_table$N_diff_codons <- NA
CANDIDATE_WINDOWS_SIZE30_table$piN_max_codon <- NA

# Additional columns
CANDIDATE_WINDOWS_SIZE30_table$num_bootstraps <- NA
CANDIDATE_WINDOWS_SIZE30_table$num_codons <- NA
#CANDIDATE_WINDOWS_SIZE30_table$N_sites_check <- NA
#CANDIDATE_WINDOWS_SIZE30_table$S_sites_check <- NA
#CANDIDATE_WINDOWS_SIZE30_table$N_diffs_check <- NA
#CANDIDATE_WINDOWS_SIZE30_table$S_diffs_check <- NA
CANDIDATE_WINDOWS_SIZE30_table$num_replicates <- NA
CANDIDATE_WINDOWS_SIZE30_table$dN_check <- NA
CANDIDATE_WINDOWS_SIZE30_table$dS_check <- NA
CANDIDATE_WINDOWS_SIZE30_table$dNdS_check <- NA
CANDIDATE_WINDOWS_SIZE30_table$dN_m_dS_check <- NA
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_SE <- NA
CANDIDATE_WINDOWS_SIZE30_table$boot_dS_SE <- NA
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_over_dS_SE <- NA
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_over_dS_P <- NA
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_m_dS_SE <- NA
CANDIDATE_WINDOWS_SIZE30_table$P_value <- NA
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_gt_dS_count <- NA
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_eq_dS_count <- NA
CANDIDATE_WINDOWS_SIZE30_table$boot_dN_lt_dS_count <- NA
CANDIDATE_WINDOWS_SIZE30_table$ASL_dN_gt_dS_P <- NA
CANDIDATE_WINDOWS_SIZE30_table$ASL_dN_lt_dS_P <- NA

# Bootstrap parameters
NBOOTSTRAPS <- 10000
NCPUS <- 6
# no min number codons because NA here

#View(CANDIDATE_WINDOWS_SIZE30_table)

# LOOP
for (i in 1:nrow(CANDIDATE_WINDOWS_SIZE30_table)) {
  #i <- 15 # 21 # 15
  cat(i, " ")
  this_row <- CANDIDATE_WINDOWS_SIZE30_table[i, ]
  this_product <- as.character(this_row$product)
  this_start <- as.integer(this_row$codon_start)
  this_end <- as.integer(this_row$codon_end)
  
  this_data_subset <- filter(codon_results_byProductCodon_WINDOWS_POOLED_WIDE, # codon_results_byProductCodon, 
                             product == this_product, # product_segment == this_product, # 
                             codon_num >= this_start,
                             codon_num <= this_end)
  
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$num_codons <- nrow(this_data_subset)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_diffs <- sum(this_data_subset$N_diffs, na.rm = TRUE)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$S_diffs <- sum(this_data_subset$S_diffs, na.rm = TRUE)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_sites <- sum(this_data_subset$N_sites, na.rm = TRUE)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$S_sites <- sum(this_data_subset$S_sites, na.rm = TRUE)
  
  # Find codons with N diffs
  #View(this_data_subset)
  this_data_subset_NDiffCodons <- filter(this_data_subset, N_diffs > 0)$codon_num
  this_data_subset_NDiffCodons_string <- paste(this_data_subset_NDiffCodons, collapse = ", ")
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_diff_codons <- this_data_subset_NDiffCodons_string
  
  # Find codon with MAX piN
  this_data_subset$piN <- this_data_subset$N_diffs / this_data_subset$N_sites
  #this_data_subset[is.nan(this_data_subset$piN), ]$piN <- NA
  this_data_subset_maxPiNCodon <- as.integer(this_data_subset[! is.na(this_data_subset$piN) & this_data_subset$piN == max(this_data_subset$piN, na.rm = TRUE), ]$codon_num)
  #View(this_data_subset)
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$piN_max_codon <- this_data_subset_maxPiNCodon
  
  #View(this_data_subset)
  
  # BOOTSTRAP THIS FULL WINDOW
  boot_dNdS <- dNdS_diff_boot_fun(this_data_subset, 'N', 'S', NBOOTSTRAPS, NCPUS)
  
  # RECORD HEADER
  boot_vector_names <- c(#'num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                         'num_replicates', 
                         'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                         'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
  
  #boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
  boot_dNdS_vector <- unlist(str_split(string = paste(boot_dNdS, sep = "\t"), pattern = "\t"))
  
  # Add names
  names(boot_dNdS_vector) <- boot_vector_names
  
  # Populate additional rows
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$num_bootstraps = NBOOTSTRAPS
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$num_codons = as.integer(boot_dNdS_vector['num_codons'])
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_sites_check = as.numeric(boot_dNdS_vector['N_sites'])
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$S_sites_check = as.numeric(boot_dNdS_vector['S_sites'])
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$N_diffs_check = as.numeric(boot_dNdS_vector['N_diffs'])
  #CANDIDATE_WINDOWS_SIZE30_table[i, ]$S_diffs_check = as.numeric(boot_dNdS_vector['S_diffs'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$num_replicates = as.integer(boot_dNdS_vector['num_replicates'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$dN_check = as.numeric(boot_dNdS_vector['dN'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$dS_check = as.numeric(boot_dNdS_vector['dS'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$dNdS_check = as.numeric(boot_dNdS_vector['dNdS'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$dN_m_dS_check = as.numeric(boot_dNdS_vector['dN_m_dS'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$P_value = as.numeric(boot_dNdS_vector['P_value'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P'])
  CANDIDATE_WINDOWS_SIZE30_table[i, ]$ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P'])
  
}

# Calculate basic estimates
CANDIDATE_WINDOWS_SIZE30_table$dN <- CANDIDATE_WINDOWS_SIZE30_table$N_diffs / CANDIDATE_WINDOWS_SIZE30_table$N_sites
CANDIDATE_WINDOWS_SIZE30_table$dS <- CANDIDATE_WINDOWS_SIZE30_table$S_diffs / CANDIDATE_WINDOWS_SIZE30_table$S_sites
CANDIDATE_WINDOWS_SIZE30_table$dNdS <- CANDIDATE_WINDOWS_SIZE30_table$dN / CANDIDATE_WINDOWS_SIZE30_table$dS

# SAVE
#write_tsv(CANDIDATE_WINDOWS_SIZE30_table, "/SARS-CoV-2-South-Africa/CANDIDATE_WINDOWS_SIZE30_table_FULL.tsv")
write_tsv(CANDIDATE_WINDOWS_SIZE30_table, "/SARS-CoV-2-South-Africa/CANDIDATE_WINDOWS_SIZE30_table_FULL2.tsv")

# RELOAD
#CANDIDATE_WINDOWS_SIZE30_table <- read_tsv("/SARS-CoV-2-South-Africa/CANDIDATE_WINDOWS_SIZE30_table_FULL.tsv")
CANDIDATE_WINDOWS_SIZE30_table <- read_tsv("/SARS-CoV-2-South-Africa/CANDIDATE_WINDOWS_SIZE30_table_FULL2.tsv")

View(CANDIDATE_WINDOWS_SIZE30_table)



