get_freqs <- function(pairs, snv)
{
  snv <- subset(snv, Sample %in% pairs, select = c(Sample, Mutation, CHROM, POS, REF, ALT, AF))
  
  # just need the snv in the samples we're looking at.
  if(nrow(snv) > 0)
  { # There are mutations.
    
    mut_table <- tidyr::spread(snv, Sample, AF, fill = 0) # a data frame with mutation down the first row and then frequency in either sample in the next 2.
    
    mut_table$source <- pairs[1] # add column with first sample ID
    mut_table$recipient <- pairs[2] # add column with second sample ID
    names(mut_table)[which(names(mut_table) == as.character(pairs[1]))] <-'freq1'
    names(mut_table)[which(names(mut_table) == as.character(pairs[2]))] <-'freq2'
    
    # This function can only be run on samples that qualified for variant identification.
    # If no variants were found in the sample then the SPECID will be missing from mut_table column and so
    # freq1 or freq2 will be missing since nothing was found we set that to 0 here and add the column.
    # equal compare will replace these cases with the reference at 1.
    if(!('freq1' %in% names(mut_table)))
    {
      mut_table$freq1 <- 0
    }
    if(!('freq2' %in% names(mut_table)))
    {
      mut_table$freq2 <- 0
    }
    mut_table <- dplyr::select(mut_table, Mutation, CHROM, POS, REF, ALT, freq1,freq2, source, recipient)
    mut_table <- mut_table[order(mut_table$CHROM, mut_table$POS),]
    
    # Fill in differences based on inferred. In this case we are left with only sites that are polymorphic in one or between the 2 samples
    #mut_table %>% dplyr::group_by(CHROM,POS) %>% dplyr::do(equal_compare(.)) -> all.freq
    
    if(nrow(mut_table) > 0)
    { # Some differences exist
      return(mut_table)
    }
    else
    { # some snv were initially present (before frequency filter) but after taking into account differences in inference the samples are the same.
      x <- tibble(Mtation = NA, freq1 = NA, freq2 = NA, source = NA, recipient = NA, chr = NA, ref = NA, "pos" = NA, var = NA)
      return(x[F,])
    }
  }
  else
  { # No variants found in either sample
    x <- tibble(Mutation = NA, freq1 = NA, freq2 = NA, donor_id = NA, recipient_id = NA, chr = NA, ref = NA, "pos" = NA, var = NA)
    return(x[F,])
  }
}


## dist_tp calls get_freqs above
dist_tp <- function(pairs, snv)
{
  data.df <- get_freqs(pairs = pairs, snv = snv)
  
  if(nrow(data.df) == 0) # This is the case when there are no variants found in either sample
  {
    d = 0
  } else
  {
    y <- as.matrix(cbind(data.df$freq1, data.df$freq2))
    d = dist(t(y), method = "manhattan")
  }
  
  return(as.numeric(d))
}

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 


## QC intrahost variants and add Ref intrahost variants
snvs_ <- data.frame()
process_intrahost_variants <- function(snvs){
  message(dim(snvs))
  ## Add ref minor alleles
  for (row in 1:nrow(snvs)) {
    if(snvs[row, "AF"] > 0.5 & snvs[row, "AF"] <= 0.95){
      ref_minor <- snvs[row, ]
      ref_minor <- ref_minor %>% mutate(Sample = Sample, CHROM = CHROM, POS = POS, REF = REF, ALT = REF, AF = (1 - AF), Depth = Depth, R1 = R1, R2 = R1, SB = SB, `R1+` = `R1+`, `R1-` = `R1-`, 
                              `R2+` = `R1+`, `R2-` = `R1-`, IMPACT = "Wild Type", `Biological Effect` = "Wild Type", GENE = GENE,
                              `Biological Effect Cat` = "Wild Type", Protein = Protein)
      snvs_ = rbind(snvs_, ref_minor)
    }else if(snvs[row, "AF"] > 0.95){
        ##There is no alternative allele here, fix at 1,
      snvs[row, "AF"] <- 1
      snvs[row, "R1"] <- 0
      snvs[row, "R1+"] <- 0
      snvs[row, "R1-"] <- 0
      }
  }
  snvs <- rbind(snvs, snvs_)
  message(dim(snvs))
  ## Filter minor alleles - only those passing the criteria should be kept
  snvs <- snvs %>% filter(AF> 0.05 & `R2+` >= 5 & `R2-` >= 5 & (`R2+`/R2) >= 0.02 & (`R2-`/R2) >= 0.02)
  message(dim(snvs))
  return(snvs)
}


process_intrahost_alt_minor <- function(snvs){
  message(dim(snvs))
  ## Add ref minor alleles
  for (row in 1:nrow(snvs)) {
    if(snvs[row, "AF"] > 0.95){
      ##There is no alternative allele here, fix at 1,
      snvs[row, "AF"] <- 1
      snvs[row, "R1"] <- 0
      snvs[row, "R1+"] <- 0
      snvs[row, "R1-"] <- 0
      snvs[row, "AF"] <- (snvs[row, "AF"] + (1 - snvs[row, "AF"]))
    }
  }
  snvs <- rbind(snvs, snvs_)
  message(dim(snvs))
  ## Filter minor alleles - only those passing the criteria should be kept
  snvs <- snvs %>% filter(AF> 0.05 & `R2+` >= 5 & `R2-` >= 5 & (`R2+`/R2) >= 0.02 & (`R2-`/R2) >= 0.02)
  message(dim(snvs))
  return(snvs)
}

## Filter minor alleles - only those passing the criteria should be kept
# filter_artifacts <- function(snv){
#   snvs <- snvs %>% filter(`R2+` >= 5 & `R2-` >= 5 & (`R2+`/R2) >= 0.02 & (`R2-`/R2) >= 0.02)
#   return(snvs)
# }


###

#complete_variants <- process_intrahost_variants(df_snvs_isnvs)
