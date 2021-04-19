samples_metadata_ <- read_xlsx("~/temp/Collaborations/covid/whd/ch1_ch3/CH1_CH3_Metadata.xlsx")

##Add dates
samples_metadata_$beastIds <- paste0(samples_metadata_$strain,"-" ,
                                     sapply(str_split(samples_metadata_$date, "-", n = 2, simplify = FALSE), `[`, 2))


##CH1
samples_metadata_beast_CH1 <- samples_metadata_ %>% filter(outbreak == "CH1")
write.table(samples_metadata_beast_CH1[, c("strain", "beastIds")], file = "../whd/beast/CH1_beastIds.txt", row.names = F, col.names = F, sep = ",", quote = F)
writexl::write_xlsx(samples_metadata_beast_CH1,  "../whd/beast/CH1_metadata_beast.xlsx")



## CH3
CH3_Epi_Ids <- read_csv("~/temp/Collaborations/covid/whd/ch1_ch3/100/CH3_Epi_Ids.txt", 
                        col_names = FALSE)
CH3_Epi_Ids_collapse <- paste0(CH3_Epi_Ids$X1, collapse = "|")
samples_metadata_beast_CH3 <- samples_metadata_ %>% filter(grepl(CH3_Epi_Ids_collapse, strain))
write.table(samples_metadata_beast_CH3[, c("strain", "beastIds")], file = "../whd/beast/CH3_beastIds.txt", row.names = F, col.names = F, sep = ",", quote = F)
write.table(samples_metadata_beast_CH3[, "strain"], file = "../whd/beast/CH3_Ids_forBeast.txt", row.names = F, col.names = F, sep = "\t", quote = F)

## Beast dates
write.table(samples_metadata_beast_CH3[, c("strain", "date")], file = "../whd/beast/CH3_datess_forBeast.txt", row.names = F, col.names = F, sep = "\t", quote = F)


writexl::write_xlsx(samples_metadata_beast_CH3, "../whd/beast/CH3_metadata_beast.xlsx")
##