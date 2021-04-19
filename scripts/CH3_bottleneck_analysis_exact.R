library(readxl)
library(kableExtra)
library(dplyr)
library(stringr)
library(purrr)
# setwd("/Applications/code/src/R/PlayGround")
setwd("/Users/sanem/temp/Collaborations/covid/houriiyah/")
df_all_ch3_pairs <- readRDS("CH3_transmission_pairs.RDS")
df_snvs_isnvs_purged <- read_xlsx("masked/Supp_table3_ch1_ch3_raw_variants_100.xlsx")

##old code
# for (row in 1:nrow(df_all_ch3_pairs)) {
#   donorId <- df_all_ch3_pairs[row, "donor"]
#   recipientId <- df_all_ch3_pairs[row, "recipient"]
#   
#   donor <- df_snvs_isnvs_purged[df_snvs_isnvs_purged$Sample == donorId, c("Sample","POS","REF","ALT","AF", "R1","R2", "Depth")]
#   colnames(donor) <- paste0("D_", colnames(donor))
#   donor <- donor %>% mutate(D_AF = as.numeric(D_AF),
#                             D_isnv = ifelse(D_AF > 0.5, D_REF, D_ALT),
#                             D_isnv_freq = ifelse(D_isnv == D_ALT, D_AF, (1 - D_AF)),
#                             D_isnv_reads = ifelse(D_isnv == D_ALT, D_R2, D_R1)) %>% 
#     filter(D_isnv_reads >= 10)
#   
#   
#   recipient <- df_snvs_isnvs_purged[df_snvs_isnvs_purged$Sample == recipientId, c("Sample","POS","REF","ALT","AF", "R1","R2", "Depth")]
#   colnames(recipient) <- paste0("R_", colnames(recipient))
#   
#   pair_001 <- left_join(donor, recipient, by = c("D_POS"="R_POS"))
#   pair_001 <- pair_001 %>% mutate(R_AF   = as.numeric(ifelse(is.na(R_AF),0, R_AF)),
#                                   R_isnv = ifelse(D_isnv == R_REF, R_REF, 
#                                                   ifelse(D_isnv == R_ALT, R_ALT, NA)),
#                                   R_isnv_freq = ifelse(R_isnv == R_ALT, R_AF, 
#                                                        ifelse(R_isnv == R_REF, (1 - R_AF), NA)),
#                                   R_isnv_reads = ifelse(R_isnv == R_ALT, R_R2, 
#                                                         ifelse(R_isnv == R_REF, R_R1, NA))
#   )
#   
#   ## if no shared isnv, the field will not be numeric, so simply skip processing that file
#   if(!is.numeric(pair_001$R_isnv_freq)){
#     next
#   }
#   message("processing... ", paste0(donorId, "_", recipientId))
#   pair_001[is.na(pair_001$R_isnv_freq), "R_isnv_freq"] <- 0
#   pair_001 <- pair_001 %>% filter(!is.na(R_isnv)) ##Only keep rows for which we have shared isnvs
#   
#   message(head(pair_001[, "R_isnv_freq"]))
#   write.table(pair_001[, c("D_isnv_freq","R_isnv_freq", "R_Depth", "R_isnv_reads")], paste0("pair_freqs_ch3/exact/1000/",donorId, "_", recipientId,"-isnv_Freqs.txt"), 
#               quote = F, row.names = F, col.names = F, sep = "\t")
# }

##New code

for (row in 1:nrow(df_all_ch3_pairs)) {
  donorId <- df_all_ch3_pairs[row, "donor_id"]
  recipientId <- df_all_ch3_pairs[row, "recipient_id"]
  
  donor <- df_snvs_isnvs_purged[df_snvs_isnvs_purged$Sample == donorId, c("Sample","POS","REF","ALT","AF", "R1","R2", "Depth", "Mutation")]
  colnames(donor) <- paste0("D_", colnames(donor))
  donor <- donor %>% filter(D_AF <= 0.5)
  
  
  
  recipient <- df_snvs_isnvs_purged[df_snvs_isnvs_purged$Sample == recipientId, c("Sample","POS","REF","ALT","AF", "R1","R2", "Depth", "Mutation")]
  colnames(recipient) <- paste0("R_", colnames(recipient))
  
  pair_001 <- left_join(donor, recipient, by = c("D_Mutation"="R_Mutation"))
  pair_001 <- pair_001 %>% mutate(R_Sample = ifelse(is.na(R_Sample), recipientId, R_Sample),
                                  R_POS = ifelse(is.na(R_POS), D_POS, R_POS),
                                  R_REF = ifelse(is.na(R_REF), D_REF, R_REF),
                                  R_ALT = ifelse(is.na(R_ALT), D_ALT, R_ALT),
                                  R_AF = ifelse(is.na(R_AF), 0, R_AF),
                                  R_R1 = ifelse(is.na(R_R1), 0, R_R1),
                                  R_R2 = ifelse(is.na(R_R2), 0, R_R2),
                                  R_Depth = ifelse(is.na(R_Depth),0,R_Depth))
  
  
  if(dim(pair_001)[1] > 0){
    write.table(pair_001[, c("D_AF","R_AF", "R_Depth", "R_R2")], paste0("pair_freqs_ch3/exact/1000/masked/",donorId, "_", recipientId,"-isnv_Freqs.txt"),
                quote = F, row.names = F, col.names = F, sep = "\t")
  }
}

cmd_isnv_freqs <- "ls -1 pair_freqs_ch3/exact/1000/masked/*-isnv_Freqs.txt"
isnv_freqsList <- system(cmd_isnv_freqs, intern = T)

df_bottle_neck_ch3 = data.frame()

for (freqs_file in isnv_freqsList) {
  #freqs_file = "pair_freqs_ch3/exact/1000/K002188_KPCOVID-0162-isnv_Freqs.txt"# Only uncoment for testing
  fname <- sapply(str_split(freqs_file, "/", n = 5, simplify = FALSE), `[`, 5)
  fname <- gsub("-isnv_Freqs.txt","",fname)
  message("BB_bottleneck processing...", fname)
  source("scripts/Bottleneck_size_estimation_exact_sj.R")
  #source("/Users/sanem/Downloads/BB_bottleneck-master/Bottleneck_size_estimation_approx.r")
  
  if(is_empty(Max_LL_bottleneck)) next
  bottle_neck <- data.frame(transPair = fname, "BottleNeck" = Max_LL_bottleneck, "lowerCI" = lower_CI_bottleneck, 
                            "upperCI" = upper_CI_bottleneck)
  df_bottle_neck_ch3 <- rbind(df_bottle_neck_ch3, bottle_neck)
}

write.table(df_bottle_neck_ch3, "pair_freqs_ch3/exact/1000/masked/df_bottle_neck_ch3_NbMax1000_raw.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

df_bottle_neck_ch3  %>%  kbl(caption = "Transmission Bottleneck Analsysis (BB_Bottleneck) - CH1 Paired Samples") %>%
  kable_paper(bootstrap_options = "striped", full_width = F)

non_convergent <- df_bottle_neck_ch3 %>% group_by(transPair) %>% count()
table(non_convergent$n)

## 5402 pairs
## 2405 - had shared isnvs

## After runnin bb bottleneck, with max 200, 
## 32 did not converge after a thousand iterations

## 2284 had results,
#non_convergent_ids <- non_convergent[non_convergent$n >= 1000, "transPair"]
#non_convergent_ids <- non_convergent[non_convergent$n >= 200, "transPair"]
#df_bottle_neck_ch3_converged <- df_bottle_neck_ch3[!df_bottle_neck_ch3$transPair %in% non_convergent_ids$transPair, ]

library(sqldf)
df_bottle_neck_ch3_trimmed <- sqldf("select transPair, max(BottleNeck) as BottleNeck, lowerCI, upperCI from df_bottle_neck_ch3 group by transPair") ##Max possible bottleneck for those that did not converge
writexl::write_xlsx(df_bottle_neck_ch3_trimmed, "pair_freqs_ch3/exact/1000/masked/df_bottle_neck_ch3_NbMax1000.xlsx")

setwd("/Applications/code/src/R/PlayGround")
