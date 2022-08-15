
library(dplyr)
library(DescTools)
library(rstan)
library(brms)
library(ggplot2)
library(GGally)
library(tidybayes)
library(tidyr)


load_sequencing_summary <- function(data_folder)
{
  rb_fastq_summary_file<-Sys.glob(file.path(data_folder, "fastq", "sequencing_summary.txt"))
  rb_fastq_summary<-read.csv(rb_fastq_summary_file, header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  
  rb_summary_file<-Sys.glob(file.path(data_folder,"sequencing_summary_*.txt"))
  rb_fast5_summary<-read.csv(rb_summary_file, header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  
  rb_end_reason<-rb_fast5_summary[c("read_id","end_reason")]
  rb_sequencing_summary<-merge(rb_fastq_summary,rb_end_reason,by.x = "read_id", by.y = "read_id")
  return(rb_sequencing_summary)
}

load_coverage_report <- function(dirlist, min_rd_len)
{
  mydirlist <- unlist(strsplit(dirlist, ","))
  coverage_report <- data.frame()
  for (dir in mydirlist)
  {
    #dir <- "D:/Plasmide/Data/20220321_1207_MN24598_FAR91003_ff83ee47"
    cov <- read.csv(file.path(dir, "results", "coverage_report.tsv"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    sequencing_summary <- load_sequencing_summary(dir)
    sequencing_summary <- sequencing_summary[sequencing_summary$end_reason == 'signal_positive',]
    sequencing_summary <- sequencing_summary[sequencing_summary$passes_filtering == 'TRUE',]
    
    dir_name <- SplitPath(dir)
    fc_name <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    exp_date <- unlist(strsplit(dir_name$filename,split = "_"))[1]
    
    sample_sheet <- read.csv(file.path(dir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    sample_sheet <- sample_sheet %>%
      mutate(Reference = case_when(Reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ Reference  ))
    
    as_channel_vec <- integer()
    norm_channel_vec <- integer()
    
    for (minutes in seq.int(from = 30, to = 24*60, by = 30))
    {
      start <- (minutes - 30) * 60
      end <- minutes * 60
      tmp_summary <- sequencing_summary[sequencing_summary$template_start < end,]
      tmp_summary <- tmp_summary[tmp_summary$template_start >= start,]
      tmp_summary <- tmp_summary[tmp_summary$sequence_length_template >= min_rd_len,]
      
      first_half <- tmp_summary[tmp_summary$channel < 257,]
      second_half <- tmp_summary[tmp_summary$channel > 256,]
      
      as_channel_vec <- append(as_channel_vec, length(unique(first_half$channel)))
      norm_channel_vec <- append(norm_channel_vec, length(unique(second_half$channel)))
      
    }
    
    as_active_channels <- 1 #median(as_channel_vec)
    normal_active_channels <- 1 #median(norm_channel_vec)
    
    
    
    depth_summary <- read.csv(file.path(dir, "plasmid_depth_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    depth_summary <- depth_summary %>%
      rename(timepoint = timesplit) %>%
      group_by(reference, timepoint) %>%
      summarise(ef_depth = (meandepth[channel == "1-256"]/as_active_channels)/(meandepth[channel == "257-512"]/normal_active_channels)) %>%
      mutate(experiment = fc_name) %>%
      mutate(date = exp_date) %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      inner_join(sample_sheet, by=c("reference" = "Reference"))
    
    
    
    
    ef_data <- data.frame()
    for (minutes in seq.int(from = 30, to = 24*60, by = 30))
    {
      
      tmp_cov <- cov[cov$timepoint == minutes,]
      
      as <- tmp_cov[tmp_cov$as == "YES",]  %>%
        mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
        mutate(bases_per_channel_as = bases / as_active_channels) %>%
        mutate(reads_per_channel_as = reads / as_active_channels)
      as_plasmid <- as[as$contig == "plasmid",]
      
      normal <- tmp_cov[tmp_cov$as == "NO",]  %>%
        mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
        mutate(bases_per_channel_norm = bases / normal_active_channels) %>%
        mutate(reads_per_channel_norm = reads / normal_active_channels)
      normal_plasmid <- normal[normal$contig == "plasmid",]
      
      tmp_ef_data <- merge(as_plasmid[c("reference","timepoint","bases_per_channel_as","reads_per_channel_as")], 
                           normal_plasmid[c("reference","timepoint","bases_per_channel_norm","reads_per_channel_norm")], 
                           by=c("reference", "timepoint")) %>%
        mutate(efactor_bases = bases_per_channel_as / bases_per_channel_norm) %>%
        mutate(efactor_reads = reads_per_channel_as / reads_per_channel_norm)
      
      ef_data <- bind_rows(ef_data, tmp_ef_data)
      
    }
    
    ef_data <- ef_data %>%
      inner_join(depth_summary) %>%
      mutate(timepoint = timepoint / 60)
    
    coverage_report <- bind_rows(coverage_report, ef_data)
  }
  return(coverage_report)
}



collect_sample_info <- function(rb_dirs, mk_dirs)
{
  
  sample_infos <- data.frame()
  
  rb_cov <- load_coverage_report(rb_dirs, opt$min_read_length)
  mk_cov <- load_coverage_report(mk_dirs, opt$min_read_length)
  
  rb_cov <- rb_cov %>%
    mutate(tool = "ReadBouncer")
  
  mk_cov <- mk_cov %>%
    mutate(tool = "MinKNOW")
  
  ef <- bind_rows(rb_cov, mk_cov)
  
  # calculate the mean enrichment factor => response variable
  rb_ef_mean <- rb_cov %>%
    select(reference, efactor_bases) %>%
    group_by(reference) %>%
    summarize(ef_bases_mean = mean(efactor_bases))
  
  mydirlist <- unlist(strsplit(rb_dirs, ","))
  
  data <- data.frame()
  for (dir in mydirlist)
  {
    #dir <- "D:/Plasmide/Data/20220324_1028_MN24598_FAR92672_9a263022"
    tmp_cov <- read.csv(file.path(dir, "results", "coverage_report.tsv"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    mapping_summary <- read.csv(file.path(dir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    sequencing_summary <- load_sequencing_summary(dir)
    
    # rejected mean read lengths
    merged1 <- sequencing_summary %>%
      filter(end_reason == "data_service_unblock_mux_change" & channel < 257) %>%
      inner_join(mapping_summary, by="read_id") %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      filter(start_time < 24*60*60)
    
    rej_mean_len <- merged1 %>%
      select(reference, sequence_length_template) %>%
      group_by(reference) %>%
      summarize(rejected_mean_length = mean(sequence_length_template))
    
    # mean library size
    merged2 <- sequencing_summary %>%
      filter(end_reason == "signal_positive" & passes_filtering == "TRUE" & channel > 256) %>%
      inner_join(mapping_summary, by="read_id") %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      filter(start_time < 24*60*60)
    
    mean_lib_size <- merged2 %>%
      select(reference, sequence_length_template) %>%
      group_by(reference) %>%
      summarize(mean_lib_size = mean(sequence_length_template))
    
    lens <- rej_mean_len %>% inner_join(mean_lib_size, by = "reference")
    
    # species abundance calculation
    yield <- tmp_cov %>%
      filter(as == "NO") %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      group_by(reference) %>%
      summarise(sum_bases = sum(bases)) %>%
      mutate(abundance = sum_bases / sum(sum_bases)) %>%
      select(reference, abundance)
    
    # difference in active channels between depletion and control region
    channel_report <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(channel_report) <- c("timepoint","act_channels_depl","act_channels_ctrl")
    channel_report$timepoint <- as.numeric(channel_report$timepoint)
    channel_report$act_channels_depl <- as.integer(channel_report$act_channels_depl)
    channel_report$act_channels_ctrl <- as.integer(channel_report$act_channels_ctrl)
   
    for (minutes in seq.int(from = 30, to = 24*60, by = 30))
    {
      start <- (minutes - 30) * 60
      end <- minutes * 60
      tmp_summary <- sequencing_summary %>%
        filter(template_start < end & template_start > start)
      
      
      first_half <- tmp_summary[tmp_summary$channel < 257,]
      second_half <- tmp_summary[tmp_summary$channel > 256,]
      
      channel_report <- channel_report %>%
        add_row(timepoint = minutes/60, act_channels_depl = n_distinct(first_half$channel),act_channels_ctrl = n_distinct(second_half$channel))
     }
    
    channel_report <- channel_report %>%
      mutate(active_diff = act_channels_ctrl - act_channels_depl)
    
    
    # sensitivity and specificity 
    
    overall_metrics <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(overall_metrics) <- c("recall","specificity", "tool", "reference")
    
    overall_metrics$recall <- as.numeric(overall_metrics$recall)
    overall_metrics$specificity <- as.numeric(overall_metrics$specificity)
    overall_metrics$tool <- as.character(overall_metrics$tool)
    overall_metrics$reference <- as.character(overall_metrics$reference)
    
    if (file.exists(file.path(dir, "other_reports","read_until_decision_stats.csv")))
    {
      ru_log <- read.csv(file.path(dir, "other_reports","read_until_decision_stats.csv"), header=TRUE, sep=";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
      ru_log_unblock <- ru_log[ru_log$channel < 257,] %>%
        group_by(read_id) %>%
        filter(decision == "unblock") %>%
        filter(row_number() == n())
      
      ru_log_other <- ru_log[ru_log$channel < 257,] %>%
        group_by(read_id) %>%
        filter(decision != "unblock") %>%
        filter(row_number() == n()) %>%
        anti_join(ru_log_unblock, by="read_id")
      
      
      ru_log <- bind_rows(ru_log_unblock, ru_log_other)
      
      ru_mapping_summary <- read.csv(file.path(dir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
      
      ru_mapping_summary <- ru_mapping_summary[ru_mapping_summary$channel < 257,]
      
      mapped_fastq <- ru_log[c("read_id","decision")] %>%
        inner_join(ru_mapping_summary, by="read_id")
      
      refs <- mapped_fastq %>%
        select(reference) %>%
        group_by(reference) %>%
        distinct(reference)
      
      
      
      for (r in unlist(refs)) 
      {
        if (r == "-")
          next
        
        matching_stats <- mapped_fastq %>%
          group_by(decision, contig, reference) %>%
          count() %>%
          filter(reference == r) %>%
          mutate(unblocked = case_when(decision == "no_decision" ~ "NO", decision == "stop_receiving" ~ "NO", decision == "unblock" ~ "YES")) %>%
          group_by(contig, unblocked) %>%
          summarise(n = sum(n))
        
        tp <- matching_stats %>%
          filter(unblocked == "YES" & contig == "chromosome") %>%
          pull(n)
        
        fp <- matching_stats %>%
          filter(unblocked == "YES" & contig != "chromosome") %>%
          group_by(unblocked) %>%
          summarise(n = sum(n)) %>%
          pull(n)
        
        tn <- matching_stats %>%
          filter(unblocked == "NO" & contig != "chromosome") %>%
          group_by(unblocked) %>%
          summarise(n = sum(n)) %>%
          pull(n)
        
        fn <- matching_stats %>%
          filter(unblocked == "NO" & contig == "chromosome") %>%
          pull(n)
        
        overall_metrics <- overall_metrics %>%
          add_row(
            recall = tp/(tp+fn),
            specificity = tn/(tn+fp),
            tool = "ReadBouncer",
            reference = r)
      }
      
    }
    
    
    
    data <- bind_rows(data, lens %>% 
                        inner_join(yield, by = "reference") %>% 
                        mutate(mean_channel_diff = mean(channel_report$active_diff)) %>%
                        left_join(overall_metrics, by = "reference")
                      
                      )
    
  }
  
  rb_ef_mean <- rb_ef_mean %>%
    left_join(data, by = "reference")
  
  
  
  # now minknow data
  
  mk_ef_mean <- mk_cov %>%
    select(reference, efactor_bases) %>%
    group_by(reference) %>%
    summarize(ef_bases_mean = mean(efactor_bases))
  
  
  mydirlist <- unlist(strsplit(mk_dirs, ","))
  
  data <- data.frame()
  for (dir in mydirlist)
  {
   # dir <- "D:/Plasmide/Data/20220419_1137_MN24598_FAR92750_43bb5c8c"
    tmp_cov <- read.csv(file.path(dir, "results", "coverage_report.tsv"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    mapping_summary <- read.csv(file.path(dir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    sequencing_summary <- load_sequencing_summary(dir)
    
    # rejected mean read lengths
    merged1 <- sequencing_summary %>%
      filter(end_reason == "data_service_unblock_mux_change" & channel < 257) %>%
      inner_join(mapping_summary, by="read_id") %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      filter(start_time < 24*60*60)
    
    rej_mean_len <- merged1 %>%
      select(reference, sequence_length_template) %>%
      group_by(reference) %>%
      summarize(rejected_mean_length = mean(sequence_length_template))
    
    # mean library size
    merged2 <- sequencing_summary %>%
      filter(end_reason == "signal_positive" & passes_filtering == "TRUE" & channel > 256) %>%
      inner_join(mapping_summary, by="read_id") %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      filter(start_time < 24*60*60)
    
    mean_lib_size <- merged2 %>%
      select(reference, sequence_length_template) %>%
      group_by(reference) %>%
      summarize(mean_lib_size = mean(sequence_length_template))
    
    lens <- rej_mean_len %>% inner_join(mean_lib_size, by = "reference")
    
    # species abundance calculation
    yield <- tmp_cov %>%
      filter(as == "NO") %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      group_by(reference) %>%
      summarise(sum_bases = sum(bases)) %>%
      mutate(abundance = sum_bases / sum(sum_bases)) %>%
      select(reference, abundance)
    
    # difference in active channels between depletion and control region
    channel_report <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(channel_report) <- c("timepoint","act_channels_depl","act_channels_ctrl")
    channel_report$timepoint <- as.numeric(channel_report$timepoint)
    channel_report$act_channels_depl <- as.integer(channel_report$act_channels_depl)
    channel_report$act_channels_ctrl <- as.integer(channel_report$act_channels_ctrl)
    
    for (minutes in seq.int(from = 30, to = 24*60, by = 30))
    {
      start <- (minutes - 30) * 60
      end <- minutes * 60
      tmp_summary <- sequencing_summary %>%
        filter(template_start < end & template_start > start)
      
      
      first_half <- tmp_summary[tmp_summary$channel < 257,]
      second_half <- tmp_summary[tmp_summary$channel > 256,]
      
      channel_report <- channel_report %>%
        add_row(timepoint = minutes/60, act_channels_depl = n_distinct(first_half$channel),act_channels_ctrl = n_distinct(second_half$channel))
    }
    
    channel_report <- channel_report %>%
      mutate(active_diff = act_channels_ctrl - act_channels_depl)
    
    
    # sensitivity and specificity 
    
    overall_metrics <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(overall_metrics) <- c("recall","specificity", "tool", "reference")
    
    overall_metrics$recall <- as.numeric(overall_metrics$recall)
    overall_metrics$specificity <- as.numeric(overall_metrics$specificity)
    overall_metrics$tool <- as.character(overall_metrics$tool)
    overall_metrics$reference <- as.character(overall_metrics$reference)
    
    if (file.exists(file.path(dir, "other_reports","read_until_decision_stats.csv")))
    {
      ru_log <- read.csv(file.path(dir, "other_reports","read_until_decision_stats.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
      ru_log_unblock <- ru_log[ru_log$channel < 257,] %>%
        group_by(read_id) %>%
        filter(decision == "unblock") %>%
        filter(row_number() == n())
      
      ru_log_other <- ru_log[ru_log$channel < 257,] %>%
        group_by(read_id) %>%
        filter(decision != "unblock") %>%
        filter(row_number() == n()) %>%
        anti_join(ru_log_unblock, by="read_id")
      
      
      ru_log <- bind_rows(ru_log_unblock, ru_log_other)
      
      ru_mapping_summary <- read.csv(file.path(dir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
      
      ru_mapping_summary <- ru_mapping_summary[ru_mapping_summary$channel < 257,]
      
      mapped_fastq <- ru_log[c("read_id","decision")] %>%
        inner_join(ru_mapping_summary, by="read_id")
      
      refs <- mapped_fastq %>%
        select(reference) %>%
        group_by(reference) %>%
        distinct(reference) 
      
      for (r in unlist(refs)) 
      {
        if (r == "-")
          next
        
        matching_stats <- mapped_fastq %>%
          group_by(decision, contig, reference) %>%
          count() %>%
          filter(reference == r) %>%
          mutate(unblocked = case_when(decision == "no_decision" ~ "NO", decision == "stop_receiving" ~ "NO", decision == "unblock" ~ "YES")) %>%
          group_by(contig, unblocked) %>%
          summarise(n = sum(n))
        
        tp <- matching_stats %>%
          filter(unblocked == "YES" & contig == "chromosome") %>%
          pull(n)
        
        fp <- matching_stats %>%
          filter(unblocked == "YES" & contig != "chromosome") %>%
          group_by(unblocked) %>%
          summarise(n = sum(n)) %>%
          pull(n)
        
        tn <- matching_stats %>%
          filter(unblocked == "NO" & contig != "chromosome") %>%
          group_by(unblocked) %>%
          summarise(n = sum(n)) %>%
          pull(n)
        
        fn <- matching_stats %>%
          filter(unblocked == "NO" & contig == "chromosome") %>%
          pull(n)
        
        overall_metrics <- overall_metrics %>%
          add_row(
            recall = tp/(tp+fn),
            specificity = tn/(tn+fp),
            tool = "MinKNOW",
            reference = r)
      }
      
      overall_metrics <- overall_metrics %>%
        mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  ))
        
      
    }
    
    
    data <- bind_rows(data, lens %>% 
                        inner_join(yield, by = "reference") %>% 
                        mutate(mean_channel_diff = mean(channel_report$active_diff)) %>%
                        left_join(overall_metrics, by = "reference")
                      )
  }
    
  mk_ef_mean <- mk_ef_mean %>%
    left_join(data, by = "reference")
    
  sample_infos <- bind_rows(rb_ef_mean, mk_ef_mean)
  return(sample_infos)
}

my_diag <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_density(fill = "firebrick4", size = 0)
}

my_lower <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) + 
    geom_smooth(method = "lm", color = "firebrick4", size = 1/3, 
                fill = "firebrick", alpha = 1/5) +
    geom_point(color = "firebrick", alpha = .8, size = 1/4)
}


median_standardize <- function(x) {
  (x- median(x) )/mad(x)
}

mean_standardize <- function(x)
{
  (x - mean(x)) / sd(x)
}

fit_model <- function(sample_infos)
{
  
  sample_infos <- infos
  
  sample_infos <- sample_infos %>% 
    drop_na() %>%
    mutate(log_lib_size = log(mean_lib_size)) %>%
    mutate(log_rej_len = log(rejected_mean_length)) %>%
    mutate(ef_mean_z = mean_standardize(ef_bases_mean)) %>%
    mutate(rej_mean_z = mean_standardize(rejected_mean_length)) %>%
    mutate(lib_mean_z = mean_standardize(mean_lib_size)) %>%
    mutate(abund_mean_z = mean_standardize(abundance)) %>%
    mutate(channel_mean_z = mean_standardize(mean_channel_diff)) %>%
    mutate(recall_mean_z = mean_standardize(recall)) %>%
    mutate(spec_mean_z = mean_standardize(specificity)) %>%
    mutate(ef_median_z = median_standardize(ef_bases_mean)) %>%
    mutate(rej_median_z = median_standardize(rejected_mean_length)) %>%
    mutate(lib_median_z = median_standardize(mean_lib_size)) %>%
    mutate(abund_median_z = median_standardize(abundance)) %>%
    mutate(channel_median_z = median_standardize(mean_channel_diff)) %>%
    mutate(recall_median_z = median_standardize(recall)) %>%
    mutate(spec_median_z = median_standardize(specificity))
  
  mean_standardized_fit_ef =  brm(data = sample_infos,
                             family = gaussian(),
                             formula = ef_mean_z ~ 0 + rej_mean_z + lib_mean_z + abund_mean_z + channel_mean_z + recall_mean_z + spec_mean_z + rej_mean_z:recall_mean_z
  )
  
  mean_standardized_fit_ef <- add_criterion(mean_standardized_fit_ef, "waic", "loo")
 
  priors = c(prior(lognormal(0, 0.25), class = Intercept, lb = 0),
            prior(normal(600, 150), class = b, coef = rejected_mean_length),
            prior(normal(0.15, 0.07), class = b, coef = abundance),
            prior(exponential(1), class = sigma),
            prior(normal(20, 9), class = b, coef = mean_channel_diff),
            prior(normal(95,3), class = b, coef = recall),
            prior(normal(95,2), class = b, coef = specificity)
            )
  
  model <- bf(ef_bases_mean ~ 0 + rejected_mean_length + mean_lib_size + 
                abundance + mean_channel_diff + recall + specificity)
  
  
  fitted_ef <- brm(data = sample_infos,
                   family = gaussian(),
                   formula = model,
                   control = list(adapt_delta = 0.9)
                  )
  
  get_prior(fitted_ef)
  
  fitted_ef <- add_criterion(fitted_ef, "waic", "loo")
  
  median_standardized_fit_ef =  brm(data = sample_infos,
                                   family = gaussian(),
                                   formula = ef_median_z ~ rej_median_z + lib_median_z + abund_median_z + channel_median_z + recall_median_z + spec_median_z
  )
  
  median_standardized_fit_ef <- add_criterion(median_standardized_fit_ef, "waic", "loo")
  
  
  mean_interaction <- brm(data = sample_infos,
                          family = gaussian(),
                          formula = ef_mean_z ~ 0 + rej_mean_z + lib_mean_z + abund_mean_z + channel_mean_z + recall_mean_z + spec_mean_z + rej_mean_z:recall_mean_z + spec_mean_z:rej_mean_z,
                          control = list(adapt_delta = 0.9)
  )
  
  
  mean_interaction <- add_criterion(mean_interaction, "waic", "loo")
  
  model_weights(mean_standardized_fit_ef, mean_interaction, fitted_ef,
                weights = "loo") %>% 
    round(digits = 2)
  
  
  rbind(bayes_R2(mean_standardized_fit_ef), 
        bayes_R2(mean_interaction), 
        bayes_R2(fitted_ef)) %>%
    as_tibble() %>%
    mutate(model = c("mean", "interaction", "normal"),
           r_square_posterior_mean = round(Estimate, digits = 2)) %>%
    select(model, r_square_posterior_mean)
  
  print(mean_standardized_fit_ef)
  print(fitted_ef)
  print(mean_interaction)
  
 # print(fitted_ef)
  posterior_summary()
  
  post <- posterior_samples(mean_standardized_fit_ef)
  
  post %>% 
    select(-lp__) %>% 
    gather() %>% 
    
    ggplot(aes(x = value, y = reorder(key, value))) +  # note how we used `reorder()` to arrange the coefficients
    geom_vline(xintercept = 0, color = "firebrick4", alpha = 1/10) +
    stat_pointintervalh(point_interval = mode_hdi, .width = .95, 
                        size = 3/4, color = "firebrick4") +
    labs(title = "My tidybayes-based coefficient plot",
         x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.grid   = element_blank(),
          panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
          axis.text.y  = element_text(hjust = 0),
          axis.ticks.y = element_blank())
  
  
  
  ggpairs(data  = sample_infos, columns = c(2:8),
          diag  = list(continuous = my_diag),
          lower = list(continuous = my_lower)) + 
    theme_bw() +
    theme(strip.background = element_rect(fill = "white", color = "white"),
          axis.text        = element_blank(),
          axis.ticks       = element_blank(),
          panel.grid       = element_blank())

  }



mkdir1 <- file.path(opt$data_dir,"20220419_1137_MN24598_FAR92750_43bb5c8c")
mkdir2 <- file.path(opt$data_dir,"20220420_1351_MN24598_FAP84921_7d466390")
rbdir1 <- file.path(opt$data_dir,"20220321_1207_MN24598_FAR91003_ff83ee47")
rbdir2 <- file.path(opt$data_dir,"20220324_1028_MN24598_FAR92672_9a263022")
rbdir3 <- file.path(opt$data_dir,"20220608_1124_MN24598_FAT03513_da85ca38")

rb_dirs <- paste(rbdir1, rbdir2, rbdir3, sep = ",")
mk_dirs <- paste(mkdir1, mkdir2, sep = ",")

dirlist <- paste(rb_dirs, mk_dirs,sep = ",")


infos <- collect_sample_info(rb_dirs, mk_dirs)



