t1 <- Sys.time()

### Load and attach required libraries ###

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(R.utils))

### Define Functions ###

# function that returns the number of "floored" values per row
# then adds to your table

get_nflr <- function(q) {
  tmp <- c()
  if (!("nflr" %in% colnames(q))) {
    for (i in 1:length(rownames(q))) {
      this_flr = q[i,]$flr
      x <- q[i,6:length(q[i,])]
      x <- x[!is.na(x)]
      tmp <- rbind(tmp, sum(x == this_flr))
    }
    q <- q %>% 
      bind_cols(as.data.frame(tmp)) %>%
      rename(nflr = V1) %>%
      relocate(nflr, .after=flr)
  }
  return(q) 
}

# function that sets floored values to NA if the number 
# of FLOORED entries is large (flr_cap default is 3)

flr.to.na <- function(q, flr_cap = 3) {
  x <- c()
  t <- q %>% select(CHR, start, end, gene_id, flr)
  for (i in 1:length(rownames(q))) {
    if (q[i,]$nflr >= flr_cap) { #if equal to or over cap, set to NA
      this_flr = q[i,]$flr
      x <- x %>% bind_rows(na_if(q[i,6:length(q[i,])],this_flr))
    } else { # if not over cap, keep as is
      x <- x %>% bind_rows(q[i,6:length(q[i,])]) 
    }
  }
  x <- t %>%
    bind_cols(x)
  return(x)
}

# function that kills the job without producing errors

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

### Assign input arguments ###

args <- commandArgs();

if (length(args) != 15) {
        print("improper number of arguments. Exiting");
        print("Usage: %> Rscript circ_eQTL.R known_chr start_interval end_interval tissue_name covariates_file time_file ge_file snp_file disttogene out_dir");
        quit(save="no");
}

if ((args[6]) == 'X') {
        known_chr <- as.character(args[6])
} else {
        known_chr <- as.numeric(args[6])
}

start_interval <- as.numeric(args[7])
end_interval <- as.numeric(args[8])
tissue_name <- as.character(args[9])
covariates_file <- as.character(args[10])
time_file <- as.character(args[11])
ge_file <- as.character(args[12])
snp_file <- as.character(args[13])
disttogene <- as.numeric(args[14])
out_dir <- as.character(args[15])

pre_si <- format(start_interval, scientific = FALSE)
pre_ei <- format(end_interval, scientific = FALSE)

### Print chromosome and region for job output ###

print(paste("Chromosome:",known_chr))
print(paste("Start position:",pre_si))
print(paste("End position:",pre_ei))

### Set up result file names ###

fn_expr <- paste(out_dir,tissue_name,"_chr",known_chr,"_",pre_si,"-",pre_ei, "_expr.txt", sep = "")
fn_base <- paste(out_dir,tissue_name,"_chr",known_chr,"_",pre_si,"-",pre_ei, "_base.txt", sep = "")
fn_baseexpr <- paste(out_dir,tissue_name,"_chr",known_chr,"_",pre_si,"-",pre_ei, "_baseexpr.txt", sep = "")
fn_full <- paste(out_dir,tissue_name,"_chr",known_chr,"_",pre_si,"-",pre_ei, "_full.txt", sep = "")
fn_models <- paste(out_dir,tissue_name,"_chr",known_chr,"_",pre_si,"-",pre_ei, "_model_comparison.txt", sep = "")
fn_nflr <- paste(out_dir,tissue_name,"_chr",known_chr,"_",pre_si,"-",pre_ei, "_nflr.txt", sep = "")

### Read in data ###

## Transcripts ##

ge <- read.table(file = ge_file, header = TRUE, comment.char = "") %>%
  rename(CHR = chr) %>%
  filter(CHR == known_chr & start >= start_interval - disttogene & end <= end_interval + disttogene)

if(nrow(ge) == 0){
  print("No transcripts in this region")
  t2<- Sys.time()
  print("Runtime:")
  print(t2 - t1)
  stop_quietly()
}

# Format Transcripts #

###############################################################
PROP_FILT_VAL <- 0.5

## take table %>% for rows rowwise, get floor value across GTEX entries
## move flr around, ungroup rowwise, count floored values (nflr),
## then get total non-NA values, proportion that are floored
## then filter out transcript that have too MANY floored values (PROP_FILT_VAL)
## excise these entries
## and finally set values that match the FLOORED value to NA if too numerous

q <- ge %>% 
  rowwise(gene_id) %>% 
  mutate(flr = min(c_across(starts_with("GTEX")), na.rm=TRUE)) %>%
  relocate(flr, .after=gene_id) %>%
  ungroup() %>%  
  get_nflr() %>%
  mutate(nmiss = rowSums(!is.na(select(.,-gene_id,-flr,-nflr, -CHR, -start, -end)))) %>%
  mutate(prop = nflr/nmiss) %>%
  filter(prop <= PROP_FILT_VAL)

if(nrow(q) == 0){
  print("No transcripts in this region")
  t2<- Sys.time()
  print("Runtime:")
  print(t2 - t1)
  stop_quietly()
}

q <- q
  select(-nmiss,-prop) %>%
  flr.to.na()
#########################################################

ge <- q
remove(q)


## SNPs w/ reference alleles ##

snps <- read.table(file = snp_file, header=TRUE) %>% filter(CHR == known_chr & POS >= start_interval & POS <= end_interval)

if(nrow(snps) == 0){
  print("No SNPs in this region")
  t2<- Sys.time()
  print("Runtime:")
  print(t2 - t1)
  stop_quietly()
}

## Covariates ##

covariates <- read.table(file = covariates_file, header = TRUE)
COV <- covariates

## Time Ordering ##

time <- read.csv(file = time_file, header = TRUE)
TD <- time %>% 
  select(ID, CIRC_TIME = Phase) %>% # Can change to other "Phase" variable
  mutate(ID = str_replace(ID, "-","."))

### Nested for loop with regressions ###

headerflag <- TRUE

for (i in 1:nrow(snps)) {
  this_snp <- snps[i,]$POS
  this_ge <- ge %>%
    filter(this_snp >= start - disttogene & this_snp <= start + disttogene)

  if (nrow(this_ge) != 0) {
    
    this_snpid <- snps[i,]$SNPID
    MAPxpr <- snps[i,]
    SNPxpr <- data.frame(ID = names(snps[i,-(1:5)]), ADD = as.numeric(snps[i,-(1:5)])) 
    
    for (j in 1:nrow(this_ge)) {
      
      GExpr <- this_ge[j,] %>%
        select(-CHR, -start, -end, -gene_id, -nflr, -flr) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var = "ID") %>%
        rename(EXPR = V1)
      
      transcript2 <- this_ge[j,]$gene_id
      this_flr <- this_ge[j,]$flr
      this_nflr <- this_ge[j,]$nflr
      
      # to join multiple dataframes, use nested left_join
      joinxpr <- left_join(GExpr, SNPxpr, by = "ID") %>%
        left_join(., COV, by = "ID") %>%
        left_join(., TD, by = "ID")
      
      
      ### Regressions
      ## basic regression Does it work? null
      null_glm <- glm(EXPR ~ . -ID -CIRC_TIME -ADD, data = joinxpr)
      
      # Q: Does genotype (mean shift) predict gene expression? ADD only
      # join genotypes and transcripts values
      # ADD == dosage of genotypes
      expr_glm <- glm(EXPR ~ . -ID -CIRC_TIME, data=joinxpr)
      
      # Does circadian cycle predict patterns of gene expression over time (theta)? circtime + covariates
      base_glm <- glm(EXPR ~ . -ID -CIRC_TIME -ADD + sin(CIRC_TIME) + cos(CIRC_TIME), data=joinxpr)
      
      
      # Does circadian cycle plus genotype (mean shift) predict gene expression over time?
      base_expr_glm <- glm(EXPR ~ . -ID -CIRC_TIME + sin(CIRC_TIME) + cos(CIRC_TIME), data=joinxpr)
      
      
      # Does circadian cycle plus genotype (mean shift) plus cycle x genotype interactionb predict gene expression over time?
      full_glm <- glm(EXPR ~ . -ID -CIRC_TIME + sin(CIRC_TIME) + cos(CIRC_TIME)  + ADD:sin(CIRC_TIME) + ADD:cos(CIRC_TIME), data=joinxpr)
      
      # Model comparisons
      m1 <- anova(null_glm,expr_glm, test="F")
      m2 <- anova(null_glm,base_glm, test="F")
      m3 <- anova(expr_glm,base_expr_glm, test="F")
      m4 <- anova(base_expr_glm,full_glm, test="F")
      
      model_results <- t(c(as.character(this_snpid),as.character(MAPxpr$REF), as.character(MAPxpr$ALT),as.character(transcript2),m1$F[2],m1$`Pr(>F)`[2],m2$F[2],m2$`Pr(>F)`[2],m3$F[2],m3$`Pr(>F)`[2],m4$F[2],m4$`Pr(>F)`[2]))
      colnames(model_results) = c("SNPID", "REF", "ALT", "GENEID", "m1_Fstat", "m1_P", "m2_Fstat", "m2_P","m3_Fstat", "m3_P","m4_Fstat", "m4_P")
      
      #floored number results
      nflr_results <- t(c(as.character(transcript2),as.character(this_flr), as.character(this_nflr)))
      colnames(nflr_results) = c("GENEID", "FLR", "NFLR")
      
      
      ## Print output of regressions to files
      # print results, but only summary stats for parameters of interest, ignoring covariates.
      
      # expr_results
      expr_results_new_add <- summary(expr_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "ADD") %>%
        select(-params)
      expr_results_new2 <- t(c(as.character(this_snpid), as.character(MAPxpr$REF), as.character(MAPxpr$ALT), as.character(transcript2))) %>%
        as.data.frame()
      expr_results_new <- merge(expr_results_new2, expr_results_new_add) %>%
        as.data.frame()
      colnames(expr_results_new) <-(c("SNPID", "REF", "ALT", "GENEID", "ADD", "ADD_SE", "ADD_tstat", "ADD_P"))
      
      
      # base_results
      base_results_new_gamma <- summary(base_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "sin(CIRC_TIME)") %>%
        select(-params)
      colnames(base_results_new_gamma) <- (c("c_GAMMA", "c_GAMMA_SE", "c_GAMMA_tstat", "c_GAMMA_P"))
      base_results_new_beta <- summary(base_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "cos(CIRC_TIME)") %>%
        select(-params)
      colnames(base_results_new_beta) <- (c("c_BETA", "c_BETA_SE", "c_BETA_tstat", "c_BETA_P"))
      base_results_new_merge <- merge(base_results_new_gamma, base_results_new_beta) %>%
        as.data.frame() 
      base_results_new2 <- t(c(as.character(this_snpid), as.character(MAPxpr$REF), as.character(MAPxpr$ALT), as.character(transcript2))) %>%
        as.data.frame()
      base_results_new <- merge(base_results_new2, base_results_new_merge)
      colnames(base_results_new) <- (c("SNPID", "REF", "ALT", "GENEID", "c_GAMMA", "c_GAMMA_SE", "c_GAMMA_tstat", "c_GAMMA_P", "c_BETA", "c_BETA_SE", "c_BETA_tstat", "c_BETA_P"))     
      
      # base_expr_results
      base_expr_results_new_add <- summary(base_expr_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "ADD") %>%
        select(-params)
      colnames(base_expr_results_new_add) <- (c("ADD", "ADD_SE", "ADD_tstat", "ADD_P"))
      base_expr_results_new_gamma <- summary(base_expr_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "sin(CIRC_TIME)") %>%
        select(-params)
      colnames(base_expr_results_new_gamma) <- (c("c_GAMMA", "c_GAMMA_SE", "c_GAMMA_tstat", "c_GAMMA_P"))
      base_expr_results_new_beta <- summary(base_expr_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "cos(CIRC_TIME)") %>%
        select(-params)
      colnames(base_expr_results_new_beta) <- (c("c_BETA", "c_BETA_SE", "c_BETA_tstat", "c_BETA_P"))
      base_expr_results_new_merge <- merge(base_expr_results_new_add, base_expr_results_new_gamma) %>%
        merge(., base_expr_results_new_beta) %>%
        as.data.frame() 
      base_expr_results_new2 <- t(c(as.character(this_snpid), as.character(MAPxpr$REF), as.character(MAPxpr$ALT), as.character(transcript2))) %>%
        as.data.frame()
      base_expr_results_new <- merge(base_expr_results_new2, base_expr_results_new_merge)
      colnames(base_expr_results_new) <- (c("SNPID", "REF", "ALT", "GENEID", "ADD", "ADD_SE", "ADD_tstat", "ADD_P", "c_GAMMA", "c_GAMMA_SE", "c_GAMMA_tstat", "c_GAMMA_P", 
                                            "c_BETA", "c_BETA_SE", "c_BETA_tstat", "c_BETA_P"))      
      
      # full_results
      full_results_new_add <- summary(full_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "ADD") %>%
        select(-params)
      colnames(full_results_new_add) <- (c("ADD", "ADD_SE", "ADD_tstat", "ADD_P"))
      full_results_new_gamma <- summary(full_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "sin(CIRC_TIME)") %>%
        select(-params)
      colnames(full_results_new_gamma) <- (c("c_GAMMA", "c_GAMMA_SE", "c_GAMMA_tstat", "c_GAMMA_P"))
      full_results_new_beta <- summary(full_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "cos(CIRC_TIME)") %>%
        select(-params)
      colnames(full_results_new_beta) <- (c("c_BETA", "c_BETA_SE", "c_BETA_tstat", "c_BETA_P"))
      full_results_new_addxgamma <- summary(full_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "ADD:sin(CIRC_TIME)") %>%
        select(-params)
      colnames(full_results_new_addxgamma) <- (c("ADDxGAMMA", "ADDxGAMMA_SE", "ADDxGAMMA_tstat", "ADDxGAMMA_P"))
      full_results_new_addxbeta <- summary(full_glm)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "params") %>%
        filter(params == "ADD:cos(CIRC_TIME)") %>%
        select(-params)
      colnames(full_results_new_addxbeta) <- (c("ADDxBETA", "ADDxBETA_SE", "ADDxBETA_tstat", "ADDxBETA_P"))
      full_results_new_merge <- merge(full_results_new_add, full_results_new_gamma) %>%
        merge(., full_results_new_beta) %>%
        merge(., full_results_new_addxgamma) %>%
        merge(., full_results_new_addxbeta) %>%
        as.data.frame() 
      full_results_new2 <- t(c(as.character(this_snpid), as.character(MAPxpr$REF), as.character(MAPxpr$ALT), as.character(transcript2))) %>%
        as.data.frame()
      full_results_new <- merge(full_results_new2, full_results_new_merge)
      colnames(full_results_new) <- (c("SNPID", "REF", "ALT", "GENEID", "ADD", "ADD_SE", "ADD_tstat", "ADD_P", "c_GAMMA", "c_GAMMA_SE", "c_GAMMA_tstat", "c_GAMMA_P", "c_BETA", 
                                       "c_BETA_SE", "c_BETA_tstat", "c_BETA_P", "ADDxGAMMA", "ADDxGAMMA_SE", "ADDxGAMMA_tstat", "ADDxGAMMA_P", "ADDxBETA",
                                       "ADDxBETA_SE", "ADDxBETA_tstat", "ADDxBETA_P")) 
      
      write.table(file = fn_expr, expr_results_new, append = T, row.names = F, col.names = headerflag, quote = F)
      write.table(file = fn_base, base_results_new, append = T, row.names = F, col.names = headerflag, quote = F)
      write.table(file = fn_baseexpr, base_expr_results_new, append = T, row.names = F, col.names = headerflag, quote = F)
      write.table(file = fn_full, full_results_new, append = T, row.names = F, col.names = headerflag, quote = F)
      write.table(file = fn_models, model_results, append = T, row.names = F, col.names = headerflag, quote = F)
      write.table(file = fn_nflr, nflr_results, append = T, row.names = F, col.names = headerflag, quote = F)
      
      headerflag <- FALSE
    }
  } else {
    # print("There were no genes associated with this SNP after filtering for this_ge.")
  }
  
}

### Print runtime to output file ###

t2<- Sys.time()
print("Runtime:")
print(t2 - t1)
