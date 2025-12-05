#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Imports ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(pacman)
p_load(GenomicRanges,IRanges,S4Vectors,rtracklayer,rio,here,readr,rtracklayer,
       tidyverse,dplyr,summarytools,janitor,stringr,fs,beepr,lubridate)
timeformat <- "%I:%M:%S %p"
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Function #1 ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
ProcessGTF <- function (gtf_path){
  ### Read GTF ###
  gtf <- rtracklayer::import(gtf_path)
  ### Process GTF ###
  ## Filtering to Autosomes 
  chr = paste0('chr',1:22)
  gtf <- gtf[seqnames(gtf) %in% chr]
  ## Filtering to Protein coding genes/Exons
  exons <- gtf[mcols(gtf)$type == "exon"]
  exons <- exons[!is.na(mcols(exons)$gene_type) & mcols(exons)$gene_type == "protein_coding"]
  ## Splitting by genes
  gr <- GenomicRanges::split(exons ,exons$gene_name)
  ### Reducing to remove overlaps and sort -> Union 
  gr <- GenomicRanges::reduce(gr)
  ### Unlist (ungrp)
  final_gr <- unlist(gr)
  ### adding the metadeta again
  gene_ids <- rep(names(gr), elementNROWS(gr))
  mcols(final_gr)$gene <- gene_ids
  ### Returns the final Genomic Range of Unionized Exons from a gene 
  return (final_gr)
}
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Function #2 ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
ProcessSEG <- function (seg_df){
  ## Tagging
  seg_df <- seg_df %>%  filter(npcr  > 10 & npaf > 5) %>% 
    mutate(
      seg_length_bp = end - start + 1,
      type = case_when(
        log2cr <= -1.0                            ~ "deep_deletion",
        log2cr >=  1.0                            ~ "amplification",
        # shallow events
        between(log2cr, -1.0, -0.2)               ~ "shallow_deletion",
        between(log2cr,  0.2,  1.0)               ~ "gain",
        # allelic imbalance categories
        abs(log2cr) < 0.2  & maf < 0.1            ~ "cnLOH",
        TRUE                                      ~ NA_character_
      ))                                      %>% 
    drop_na()                                 %>%
    mutate(seg_idx = row_number())
  # Converting to GRanges
  seg_gr <- makeGRangesFromDataFrame(
    seg_df,
    seqnames.field     = "contig",
    start.field        = "start",
    end.field          = "end",
    keep.extra.columns = TRUE               # Will keep other cols
  )
  return (list(df = seg_df,
               gr = seg_gr ))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Function #3 ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
ProcessOverlap <- function (case_id,gtf_gr,seg_gr,seg_df){
  ### Finding Overlaps ####
  hits <- findOverlaps(gtf_gr, seg_gr, ignore.strand = TRUE)
  ### Compute overlapping bases ###
  ov_width <- width(
    pintersect(
      ranges(gtf_gr)[queryHits(hits)],
      ranges(seg_gr)[subjectHits(hits)]))
  ## idx DF 
  hit_df <- data.frame(
    gene      = mcols(gtf_gr)$gene[queryHits(hits)],
    gene_idx  = queryHits(hits),
    seg_idx   = mcols(seg_gr)$seg_idx[subjectHits(hits)],
    overlap_bp = ov_width
  )
  ## Total overlap of a gene with a seg
  hit_df <- hit_df %>%
    group_by(gene, seg_idx) %>%
    summarise(overlap_bp = sum(overlap_bp), .groups = "drop")
  ## Filter to keep a gene only in one segment with max overlap
  hit_df <- hit_df %>%
    group_by(gene) %>%
    slice_max(overlap_bp, n = 1,with_ties = FALSE) %>%  ungroup()
  ## Merge ##
  final_df <- hit_df %>%
    left_join(seg_df, by = "seg_idx") %>%
    select(
      seg_idx,
      gene,
      contig, start, end,seg_length_bp,
      npcr,
      log2cr,
      npaf,
      maf,
      type
    ) %>% mutate(case_id = case_id)
  return(final_df)
}
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Function #4 ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
Readseg <- function(seg_path){
  ### Reads the Seg tsv and returns changes col names ###
  seg_df <- read_tsv(seg_path, comment = "@",show_col_types = FALSE)        %>% 
    janitor::clean_names()                                                  %>% 
    select(contig:num_points_allele_fraction,ends_with('50'))               %>% 
    rename(log2cr  = log2_copy_ratio_posterior_50,
           maf     = minor_allele_fraction_posterior_50, 
           npaf    = num_points_allele_fraction,
           npcr    = num_points_copy_ratio)                                 %>% 
    drop_na(maf)
  return(seg_df)
}
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Main ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
start_time <- now()
gtf_path  <- here("data","gencode.v36.annotation.gtf")
input_dir <- here("results","final","seg")
#### Step #1 ####
print(str_glue("[ {format(now(),timeformat)} ] Starting ....\n"))
gtf_gr <- ProcessGTF(gtf_path)
final_df <- NULL
#### Step #2 ####
for (path in fs::dir_ls(input_dir,recurse = T,glob = "*.modelFinal.seg")) {
  case_id <- fs::path_file(fs::path_dir(path))
  print(str_glue("[ {format(now(),timeformat)} ] Working on {case_id}\n"))
  seg_df <- Readseg(path)
  seg_list <- ProcessSEG(seg_df)
  ## For the first run final df is empty so it gets new_df
  new_df <- ProcessOverlap(case_id,gtf_gr,seg_list$gr,seg_list$df)
  final_df <- bind_rows(final_df,new_df)
  print(str_glue("[ {format(now(),timeformat)} ] Done!!\n"))
}
rio::export(final_df,here('results','final','summary_modelseg.csv'))
time_taken <- now() - start_time
print(str_glue(" Time Taken : {round(time_taken,2)} seconds"))
beepr::beep(11)








