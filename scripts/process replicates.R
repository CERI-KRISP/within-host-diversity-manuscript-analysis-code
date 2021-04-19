library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

cmd7 <- "ls -1 ~/temp/Collaborations/covid/whd/Replicates/*_sf.table"
snpSiftListR <- system(cmd7, intern = T)

df_snpsR = data.frame()

for (sfreportR in snpSiftListR) {
  psR <- read_delim(sfreportR, "\t", escape_double = FALSE, col_types = cols(.default = "c"), trim_ws = TRUE)
  #message("Processing: ", sfreport)
  
  fnameR <- sapply(str_split(sfreportR, "/", n = 9, simplify = FALSE), `[`, 9)
  fnameR <- sapply(str_split(fnameR, "_", n = 2, simplify = FALSE), `[`, 1)
  
  psR$Sample <- fnameR
  
  df_snpsR <- rbind(df_snpsR, psR)
  #message("...done!")
}


colnames(df_snpsR)[c(6,9:11)] <- c("Depth","Biological Effect", "IMPACT", "GENE")

df_snpsR$`Biological Effect` <- gsub("_variant","",df_snpsR$`Biological Effect`)

df_snpsR <- df_snpsR %>% separate(DP4, c("R1+", "R1-", "R2+", "R2-"))


df_snpsR[, c(2,5:11)] <- apply(df_snpsR[, c(2,5:11)], 2, as.numeric)
df_snpsR$R1 <- df_snpsR$`R1+` + df_snpsR$`R1-`
df_snpsR$R2 <- df_snpsR$`R2+` + df_snpsR$`R2-`

df_snvs_isnvsR <- df_snpsR[,c(15,1:6,16:17,7:14)]
df_snvs_isnvsR[, c(3,6:13)] <- apply(df_snvs_isnvsR[, c(3,6:13)],2,as.numeric)



df_snvs_isnvsR <- df_snvs_isnvsR  %>%
 mutate(`Biological Effect Cat` = case_when(`Biological Effect` == "missense" ~ "non-synonymous",
                                            `Biological Effect` == "stop_gained" ~ "nonsense",
                                            `Biological Effect` == "stop_lost" ~ "nonsense",
                                            `Biological Effect` == "start_lost" ~ "nonsense",
                                            `Biological Effect` == "splice_region&synonymous" ~ "synonymous",
                                            `Biological Effect` == "stop_lost&splice_region" ~ "nonsense",
                                            TRUE ~ `Biological Effect`))


##filter out multiallelic positions
#df_snvs_isnvs %>% count(Sample,POS,REF) %>% filter(n > 1) %>% View()
df_snvs_isnvsR <- df_snvs_isnvsR %>% mutate(alleles = paste0(Sample,"_",POS,"_",REF))
dupps <- df_snvs_isnvsR %>% group_by(Sample,POS,REF) %>% count() %>% filter(n > 1) %>% mutate(alleles = paste0(Sample,"_",POS,"_",REF))
df_snvs_isnvs_fillR <- df_snvs_isnvsR %>% 
  filter(!(alleles %in% dupps$alleles))
df_snvs_isnvs_fillR$alleles <- NULL


## protein length
c19_protein_locations <- read.table(file='c19_genome_protein_lengths.txt',header=FALSE,
                                    sep=':',col.names=c('Protein','position'))
c19_protein_locations <- c19_protein_locations %>% 
  mutate(Protein = str_trim(Protein),
         position =  str_replace(position, "\\n", "|")) %>% 
         separate(position, sep = "\\|", into = c("begin", "end"))

rownames(c19_protein_locations) <- c19_protein_locations$Protein


df_snvs_isnvs_annR <- sqldf::sqldf("select df_snvs_isnvs_fillR.*, protein from df_snvs_isnvs_fillR left join c19_protein_locations on df_snvs_isnvs_fillR.POS >= begin and df_snvs_isnvs_fillR.POS <= end")
#dim(df_snvs_isnvs_ann) #4215


## keep only one entry for each allele
#df_snvs_isnvs_ann %>% count(Sample,POS,REF,ALT,AF) %>% filter(n > 1) %>% View()
df_snvs_isnvs_fillR_ <- df_snvs_isnvs_annR %>% distinct(Sample,POS,REF,ALT,AF, .keep_all = TRUE)

df_snvs_isnvs_ann_fillR_ <- df_snvs_isnvs_fillR_ %>%  mutate( Protein = factor(Protein, levels = rownames(c19_protein_locations)),
                                                            GENE = factor(GENE, levels = c19_genes))
#dim(df_snvs_isnvs_ann_fil) #4165
df_snvs_isnvs_ann_fillR__ <- process_intrahost_variants(df_snvs_isnvs_ann_fillR_) 

### Separate minor allelels from 
df_isnvsR <- df_snvs_isnvs_ann_fillR__ %>%  filter( AF <= 0.5) ## 31
df_snvsR <- df_snvs_isnvs_ann_fillR__ %>% filter(AF > 0.5 ) ## 76

dim(df_isnvsR)
View(df_isnvsR)

## Compare no. of minor alleles in each sample to genome coverage
alleleCount_all_isnvsR <- df_isnvsR %>% 
  count(Sample)
colnames(alleleCount_all_isnvsR) <- c("Seq", "iSNVs")
alleleCount_all_isnvsR$Seq <- sapply(str_split(alleleCount_all_isnvsR$Seq, "_", n = 2, simplify = FALSE), `[`, 1)

## Compare no. of major alleles in each sample to genome coverage
alleleCount_all_snvsR <- df_snvsR %>% 
  count(Sample)
colnames(alleleCount_all_snvsR) <- c("Seq", "SNVs")
alleleCount_all_snvsR$Seq <- sapply(str_split(alleleCount_all_snvsR$Seq, "_", n = 2, simplify = FALSE), `[`, 1)

df_snvs_isnvs_ann_fillR__ <- df_snvs_isnvs_ann_fillR__ %>% mutate(repNo = as.factor(ifelse(grepl("-2", Sample), 2,1)),
                                            Sample = gsub("-2", "", Sample),
                                            Mutation = paste0(REF,POS,ALT))


## join allele counts and coverage details
#cov_ <- left_join(cov_, alleleCount_all_isnvsR, by = c("Seq" = "Seq")) #Check why some samples do not have minor alleles?
#cov_ <- left_join(cov_, alleleCount_all_snvsR, by = c("Seq" = "Seq")) #Check why some samples do not have minor alleles?
#cov_ <- cov_ %>% mutate(iSNVs = ifelse(is.na(iSNVs), 0, iSNVs))

#genomicStatsR <- cov_[, c(8:12,2:7,14:15)]

#genomicStatsR[is.na(genomicStatsR$iSNVs), "iSNVs"] <- 0
#genomicStatsR[is.na(genomicStatsR$SNVs), "SNVs"] <- 0
#genomicStatsR  <-  genomicStatsR %>% mutate(TotalAlleles = iSNVs + SNVs)
# write_xlsx(genomicStats, "Supp_table2_genomic_stats.xlsx")

## Genome stats purged
#genomicStats_purgedR <- genomicStatsR[genomicStatsR$genomecov >= 90, ] #2 samples dropped with coverage less than 90
#genomicStats_purged <- genomicStats_purged[genomicStats_purged$DepthCov10PercGT80 == 1, ]
# write.table(genomicStats_purged, "Supp_table2_genomicStats_ch1_ch3.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
#write_xlsx(genomicStats_purgedR, "Supp_table5_sequencingReplicates.xlsx")

### Exploratory Plots and Stats




### Compare frequencies of replicates:
## http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
# ggpaired(ToothGrowth, x = "supp", y = "len",
#          color = "supp", line.color = "gray", line.size = 0.4,
#          palette = "jco")+
#   stat_compare_means(paired = TRUE)
# 
# ggpaired(filter(df_snvs_isnvsR, Sample=="CC00686518T6"), x = "repNo", y = "AF",
#          color = "POS", line.color = "gray", line.size = 0.4,
#          palette = "jco")+
#   facet_grid(rows = vars(Sample))+
#   stat_compare_means(paired = TRUE)

ggplot(df_snvs_isnvs_ann_fillR__, aes(x = repNo, y = AF)) +
  geom_line(aes(group = POS, color = factor(POS)), show.legend = FALSE) +
  facet_grid(cols = vars(Sample)) +
  theme(legend.position='none') +
  xlab("Replicate No.") + 
  ylab("Allele Frequency")+
  theme_light()


df_snvs_isnvsR_sprd <-  df_snvs_isnvs_ann_fillR__ %>% 
        select(c(Sample, Mutation, repNo, AF )) %>% 
        spread(key = "repNo", value = "AF", fill = 0)

colnames(df_snvs_isnvsR_sprd)[c(3,4)] <- c("A", "B")

df_snvs_isnvsR_sprd <- df_snvs_isnvsR_sprd %>% filter(!(Sample == "CC00686518T6" & Mutation %in% c("A22206G", "C1059T", "G22813T", "C28887T")) &
                                                        !(Sample == "EG00465499T6" & Mutation %in% c( "C1059T", "C28887T")))
ggscatter(df_snvs_isnvsR_sprd, x = "A", y = "B",
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE # Add confidence interval
  )+ xlab("Replicate 2 AF") +
  ylab("Replicate 1 AF") + 
  facet_grid(cols = vars(Sample)) +
  stat_cor(method = "pearson") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_bw()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text=element_text(size=12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))


write_xlsx(df_snvs_isnvs_ann_fillR__, "df_snvs_isnvs_ann_fillR__.xlsx")

