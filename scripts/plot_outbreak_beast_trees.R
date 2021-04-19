library(ggplot2)
library(ape)
library(repr)
library("readxl")
library('gridExtra')
library(dplyr)
library(hrbrthemes)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library("ggsci")
library(ggalt)
library("Hmisc")
library(ggtree)
#library(tidyverse)
library(treeio)
library("lubridate")
library(readr)
library("wesanderson")
library(stringr)

#https://github.com/GuangchuangYu/plotting_tree_with_data/blob/master/supplemental_file.Rmd

tree <- read.beast("../whd/beast/CH1_purged_final_BeastSet_mafft.aln_strict_constant_est.tree")

##Tip labels
ch1_tips <- tree@phylo$tip.label
#get.tree(tree)$tip.label

ch1_sample_df <- data.frame(sample = ch1_tips, 
                            tip_dates_decimal = as.numeric(sapply(str_split(ch1_tips, "/", n = 6, simplify = FALSE), `[`, 6)))
ch1_sample_df <- ch1_sample_df %>% mutate(tip_dates = format(date_decimal(tip_dates_decimal), "%Y-%m-%d"),
                                          krisp_id = sapply(str_split(ch1_tips, "/", n = 3, simplify = FALSE), `[`, 2),
                                          outbreak_id = gsub("CH1-", "", sapply(str_split(ch1_tips, "/", n = 5, simplify = FALSE), `[`, 4)),
                                          new_tip_label = paste0(krisp_id, "_", outbreak_id))

ch1_sample_df[ch1_sample_df$krisp_id == "KRISP-0035", "outbreak_id"] <- "R1A"
ch1_sample_df[ch1_sample_df$krisp_id == "KRISP-0041", "outbreak_id"] <- "R1A"
ch1_sample_df[ch1_sample_df$krisp_id == "KRISP-0118", "outbreak_id"] <- "ND0118"
ch1_sample_df[ch1_sample_df$krisp_id == "KRISP-0119", "outbreak_id"] <- "ND0119"
ch1_sample_df[ch1_sample_df$krisp_id == "KRISP-0120", "outbreak_id"] <- "ND0120"
ch1_sample_df[ch1_sample_df$krisp_id == "KRISP-0121", "outbreak_id"] <- "ND0121"
ch1_sample_df[ch1_sample_df$krisp_id == "KRISP-0122", "outbreak_id"] <- "ND0122"


ch1_dept_lin <- read_delim("~/temp/Collaborations/covid/whd/ch1_ch3/ch1_dept_lin.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)

ch1_dept_lin$BeastId <- as.character(ch1_dept_lin$BeastId)
ch1_sample_df <- inner_join(ch1_sample_df, ch1_dept_lin, by =c("sample" = "BeastId"))
                                        
ch1DeptPanel <- c("ICU A"="#dfc9f0", "ICU C"="#fffe00", "Ward A"="#bdd7ee", "Other"="#79d2ee", "Ward D"="#ac0454", "Ward C"="#fcc004", "Nursing Home"="#7f7f7f", "Emerg. Dept."="#a9d08e")

p_ch1 <- ggtree(tree, mrsd="2020-04-7", as.Date=TRUE,color='gray30',size=1) + theme_tree2()  #+ geom_tiplab()
p_ch1 <- p_ch1 %<+% ch1_sample_df + 
            scale_fill_manual(name="Department", values = ch1DeptPanel)+
            geom_tippoint(aes(shape = PangolinLineage, fill = Department), size=4) +
            geom_tiplab(aes(label=outbreak_id),  linesize=.3, hjust = -0.3) +
            #scale_fill_manual(name='PangolinLineage', values = c("#9bc2e6","#548235"))+
            #scale_color_manual(name = "Department", values = c("#dfc9f0", "#fffe00", "#bdd7ee", "#79d2ee", "#ac0454", "#fcc004", "#7f7f7f", "#a9d08e")) +
            theme(axis.text=element_text(size=10),
                  legend.position = "right")+
            scale_shape_manual(name = "Pangolin Lineage", values = c(21:23)) +
            guides(fill = guide_legend(override.aes = list(shape = 21)))+
            theme(legend.position = "none")

p_ch1

CH1_tempest_data <- read_delim("CH1_tempest.data.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
CH1_tempest_data$krisp_id <- sapply(str_split(CH1_tempest_data$tip, "_", n = 3, simplify = FALSE), `[`, 2)

## add dept and lineage info
CH1_tempest_data <- inner_join(CH1_tempest_data, ch1_sample_df[, c("krisp_id", "Department", "PangolinLineage")], by = "krisp_id")

ggscatter(CH1_tempest_data, x = "distance", y = "residual",
          fill = "Department", shape = "PangolinLineage",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
)+ xlab("Distance") +
  ylab("Residual") + 
  #facet_grid(cols = vars(Sample)) +
  stat_cor(method = "pearson") + 
  theme_bw()+
  scale_fill_manual(name="Department", values = ch1DeptPanel)+
  scale_shape_manual(name = "Pangolin Lineage", values = c(21:23)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text=element_text(size=12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "left")



### CH3 Tree
treeCH3 <- read.beast("../whd/beast/CH3_new_purged_final_BeastSet_mafft.aln_strict_constant_est.tree")
##get.tree(treeCH3)$tip.label
ch3_tips <- treeCH3@phylo$tip.label
ch3_sample_df <- data.frame(sample = ch3_tips, 
                            tip_dates_decimal = as.numeric(sapply(str_split(ch3_tips, "/", n = 6, simplify = FALSE), `[`, 6)))
ch3_sample_df <- ch3_sample_df %>% mutate(tip_dates = format(date_decimal(tip_dates_decimal), "%Y-%m-%d"),
                                          krisp_id = sapply(str_split(ch3_tips, "/", n = 3, simplify = FALSE), `[`, 2),
                                          outbreak_id = gsub("CH3-", "", sapply(str_split(ch3_tips, "/", n = 5, simplify = FALSE), `[`, 4)),
                                          new_tip_label = paste0(krisp_id, "_", outbreak_id))

ch3_dept_lin <- read_delim("~/temp/Collaborations/covid/whd/ch1_ch3/ch3_dept_lin.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)
#ch3_dept_lin$BeastId <- as.character(ch3_dept_lin$BeastId)

ch3_sample_df <- inner_join(ch3_sample_df, ch3_dept_lin, by =c("sample" = "BeastId"))

ch3DeptPanel <- c("Recovery Room"="#9bc2e6", "Theatre Cleaner"="#f3ecfb", "Theatre Nurse"="#f4e9d0", 
                  "Other"="#79d2ee", "HR Dept"="#c9c9c9", "Stores"="#f2f2f3", "TSSU"="#e5eef7", 
                  "Packing Room"="#e2efda")

p_ch3 <- ggtree(treeCH3, mrsd="2020-06-8", as.Date=TRUE,color='gray30',size=1) + theme_tree2()  #+ geom_tiplab()
p_ch3 <- p_ch3 %<+% ch3_sample_df + 
  scale_fill_manual(name="Department", values = ch3DeptPanel)+
  geom_tippoint(aes(shape = PangolinLineage, fill = Department), size=4) +
  geom_tiplab(aes(label=outbreak_id),  linesize=.3, hjust = -0.3) +
  #scale_fill_manual(name='PangolinLineage', values = c("#9bc2e6","#548235"))+
  #scale_color_manual(name = "Department", values = c("#dfc9f0", "#fffe00", "#bdd7ee", "#79d2ee", "#ac0454", "#fcc004", "#7f7f7f", "#a9d08e")) +
  theme(axis.text=element_text(size=10),
        axis.title.y = element_text(size=15),
        legend.position = "right")+
  scale_shape_manual(name = "Pangolin Lineage", values = c(21:24)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  theme(legend.position = "none")

p_ch3

CH3_tempest_data <- read_delim("CH3_tempest.data.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
CH3_tempest_data$krisp_id <- sapply(str_split(CH3_tempest_data$tip, "_", n = 3, simplify = FALSE), `[`, 2)

## add dept and lineage info
CH3_tempest_data <- inner_join(CH3_tempest_data, ch3_sample_df[, c("krisp_id", "Department", "PangolinLineage")], by = "krisp_id")

ggscatter(CH3_tempest_data, x = "distance", y = "residual",
          fill = "Department", shape = "PangolinLineage",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
)+ xlab("Distance") +
  ylab("Residual") + 
  #facet_grid(cols = vars(Sample)) +
  stat_cor(method = "pearson") + 
  scale_fill_manual(name="Department", values = ch3DeptPanel)+
  scale_shape_manual(name = "Pangolin Lineage", values = c(21:24)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text=element_text(size=12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "left")

