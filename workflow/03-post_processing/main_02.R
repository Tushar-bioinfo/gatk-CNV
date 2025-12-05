#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Imports ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(pacman)
p_load(rio,here,tidyverse,dplyr,summarytools,janitor,stringr)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Load Data ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
seg <- rio::import(here("results",'final','summary_modelseg.csv'))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Main ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## Unique patients 
n_cases <- seg %>% pull(case_id) %>% n_distinct()

# fq_tbl <- seg %>% group_by(gene) %>% mutate(
#   gainfq      = (sum(type == "gain"            )/29) * 100,
#   ampfq       = (sum(type == "amplification"   )/29) * 100,
#   deepdelfq   = (sum(type == "deep_deletion"   )/29) * 100,
#   shlldelfq   = (sum(type == "shallow_deletion")/29) * 100,
#   cnLOHfq     = (sum(type == "cnLOH"           )/29) * 100
# ) %>% arrange(gene) %>% ungroup()

gene_freq <- seg %>% 
  group_by(gene) %>% 
  summarise(
    gain_n           = sum(type == "gain"),
    amp_n            = sum(type == "amplification"),
    deepdel_n        = sum(type == "deep_deletion"),
    shalldel_n       = sum(type == "shallow_deletion"),
    cnLOH_n          = sum(type == "cnLOH"),
    gainfq           = 100 * gain_n    / n_cases,
    ampfq            = 100 * amp_n     / n_cases,
    deepdelfq        = 100 * deepdel_n / n_cases,
    shlldelfq        = 100 * shalldel_n/ n_cases,
    cnLOHfq          = 100 * cnLOH_n   / n_cases,
    n_cases          = n_distinct(case_id),
    .groups = "drop"
  )

amp_list  <- gene_freq %>% arrange(desc(ampfq)) %>% head(20) %>% pull(gene) %>% str_c("amp_",.)
ddel_list <- gene_freq %>% arrange(desc(deepdelfq)) %>% head(10) %>% pull(gene) %>% str_c("ddel_",.)
gene_freq %>% arrange(desc(ampfq)) %>% head(20)
gene_freq %>% arrange(desc(deepdelfq)) 

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Function ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
surv_func <- function (gene,type,tag,type1,type2) {
  as.integer(any(gene == tag & type %in% c(type1,type2)))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
### Main ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# lst <- c(amp_list,ddel_list)
lst <- c(amp_list)

final_df <- NULL
for (i in lst){ 
### Creating a label
    gene_name <- str_replace(i,"^(amp_|ddel_)","")
    print(str_glue('Working on {gene_name} ...'))
    if (startsWith(i , "amp")) {
      temp_df  <- seg %>% 
        group_by(case_id) %>% 
        summarise(
          !!i:= surv_func(gene,type,gene_name,"amplification","gain"))
    } else {
      temp_df  <- seg %>% 
        group_by(case_id) %>% 
        summarise(
          !!i := surv_func(gene,type,gene_name,"deep_deletion","shallow_deletion"))
    }
  ## For the first iteration
  if (is.null(final_df)) {
    final_df <- temp_df 
  } else {
    final_df <- final_df %>% left_join(temp_df, by = "case_id")
 }
}
### Write ####
## Results to a csv 
rio::export(final_df,here('results','final','analysis','final.csv'))

  