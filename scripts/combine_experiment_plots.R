if(!require(optparse)) install.packages("optparse")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(scales)) install.packages("scales")
if(!require(ggsci)) install.packages("ggsci")
if(!require(gridExtra)) install.packages("gridExtra")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(DescTools)) install.packages("DescTools")
if(!require(ggh4x)) install.packages("ggh4x")

library(optparse)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(grid)
library(gridExtra)
library(DescTools)
library(ggpubr)
library(ggh4x)

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

create_ef_plots <- function(rb_cov, mk_cov, output_file)
{
  
  rb_cov <- rb_cov %>%
    mutate(tool = "ReadBouncer")
  
  mk_cov <- mk_cov %>%
    mutate(tool = "MinKNOW")
  
  ef <- bind_rows(rb_cov, mk_cov)
  
  
  ef <- ef %>%
    mutate(experiment = 
             case_when(experiment == "FAR91003" ~ "ReadBouncer1",
                       experiment == "FAR92672" ~ "ReadBouncer2",
                       experiment == "FAR92750" ~ "MinKNOW1",
                       experiment == "FAP84921" ~ "MinKNOW2",
                       TRUE ~ experiment)
          ) %>%
    mutate(Species = 
             case_when(Species == "Campylobacter coli" ~ "C. coli",
                       Species == "Campylobacter jejuni" ~ "C. jejuni",
                       Species == "Salmonella enterica" ~ "S. enterica",
                       Species == "Klebsiella pneumoniae" ~ "K. pneumoniae",
                       Species == "Enterobacter hormaechei" ~ "E. hormaechei",
                       TRUE ~ Species)
    )
  
  
  sort <- ef %>%
    select(experiment, date) %>%
    filter(!duplicated(.)) %>%
    arrange(date)
  
  
  
  ef_yield <- ggplot(ef, aes(x=timepoint, y=efactor_bases, group=Species)) +
    theme_minimal()+
    geom_line(aes(color=Species), size=1.0)+
    xlim(c(0, 25))+
    ylim(c(0, 2.5))+
    #geom_point(aes(color=design))+
    scale_fill_jco()+
    scale_color_brewer(palette="Set2")+
    labs(title ="(a)",x="Time (hours)", y = "Enrichment by yield")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    theme(legend.position = "none")+
    #theme(legend.title = element_blank())+
    facet_nested(~ factor(experiment, levels = sort$experiment))+
    theme(strip.text.x = element_text(size = 8, face="bold"))+
    theme(strip.text.y = element_text(size = 8, face="bold"))+
    theme(strip.background=element_rect(color="grey30"))
    
    #facet_wrap(~factor(experiment, levels = sort$experiment), ncol = 4)+
  
  ef_reads <- ggplot(ef, aes(x=timepoint, y=efactor_reads, group=Species)) +
    theme_minimal()+
    geom_line(aes(color=Species), size=1.0)+
    xlim(c(0, 25))+
    ylim(c(0, 2.5))+
    #geom_point(aes(color=design))+
    scale_fill_jco()+
    scale_color_brewer(palette="Set2")+
    labs(title ="(b)",x="Time (hours)", y = "Enrichment by reads")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    theme(legend.position = "none")+
    #theme(legend.title = element_blank())+
    facet_nested(~ factor(experiment, levels = sort$experiment))+
    #facet_wrap(~factor(experiment, levels = sort$experiment), ncol = 4)+
    #theme(strip.text.x = element_text(size = 8, face="bold"))
    theme(strip.text.x = element_blank())
  
  ef_depth <- ggplot(ef, aes(x=timepoint, y=ef_depth, group=Species)) +
    theme_minimal()+
    geom_line(aes(color=Species), size=1.0)+
    xlim(c(0, 25))+
    ylim(c(0, 2.5))+
    #geom_point(aes(color=design))+
    scale_fill_jco()+
    scale_color_brewer(palette="Set2")+
    labs(title ="(c)",x="Time (hours)", y = "Enrichment by depth")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    theme(legend.position = "bottom")+
    theme(legend.title = element_blank())+
    facet_nested(~ factor(experiment, levels = sort$experiment))+
    #facet_wrap(~factor(experiment, levels = sort$experiment), ncol = 4)+
    #theme(strip.text.x = element_text(size = 8, face="bold"))
    theme(strip.text.x = element_blank())
  
  # creates Figure 4
  combined <- grid.arrange(ef_yield, ef_reads, ef_depth, nrow = 3)
  ggsave(filename=output_file, 
         plot = combined, 
         device = cairo_pdf, 
         width = 210, 
         height = 297, 
         units = "mm")
  
  dev.off()
}

create_matching_stats <- function(rbdir, mkdir, rb2_dir, output_file)
{
 
  overall_metrics <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(overall_metrics) <- c("accuracy","precision","recall","specificity", "f1_score","balanced_acc","tool", "reference")
  overall_metrics$accuracy <- as.numeric(overall_metrics$accuracy)
  overall_metrics$precision <- as.numeric(overall_metrics$precision)
  overall_metrics$recall <- as.numeric(overall_metrics$recall)
  overall_metrics$specificity <- as.numeric(overall_metrics$specificity)
  overall_metrics$f1_score <- as.numeric(overall_metrics$f1_score)
  overall_metrics$balanced_acc <- as.numeric(overall_metrics$balanced_acc)
  #overall_metrics$mcc <- as.numeric(overall_metrics$mcc)
  overall_metrics$tool <- as.character(overall_metrics$tool)
  overall_metrics$reference <- as.character(overall_metrics$reference)

  # first mkdir  
  ru_log <- read.csv(file.path(mkdir1, "other_reports","read_until_decision_stats.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
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
    
  ru_mapping_summary <- read.csv(file.path(mkdir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    
  ru_mapping_summary <- ru_mapping_summary[ru_mapping_summary$channel < 257,]
    
  sequencing_summary <- load_sequencing_summary(mkdir)
  
  mapped_fastq <- ru_log[c("read_id","decision")] %>%
    inner_join(ru_mapping_summary, by="read_id")

  matching_stats <- mapped_fastq %>%
    group_by(decision, contig) %>%
    count() %>%
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
      accuracy = (tp+tn)/(tp+tn+fp+fn),
      precision = tp/(tp+fp),
      recall = tp/(tp+fn),
      specificity = tn/(tn+fp),
      f1_score = 2*tp/(2*tp + fp + fn),
      balanced_acc = (recall + specificity)/2,
     # mcc = (as.double(tp)*as.double(tn) - fp*fn)/sqrt(as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)),
      tool = "MinKNOW",
      reference = "all")
  
  # now for each species
  
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
        accuracy = (tp+tn)/(tp+tn+fp+fn),
        precision = tp/(tp+fp),
        recall = tp/(tp+fn),
        specificity = tn/(tn+fp),
        f1_score = 2*tp/(2*tp + fp + fn),
        balanced_acc = (recall + specificity)/2,
        #mcc = (as.double(tp)*as.double(tn) - fp*fn)/sqrt(as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)),
        tool = "MinKNOW",
        reference = r)
    
  }
  
  sample_sheet <- read.csv(file.path(mkdir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  
  # now add readbouncer metrics
  
  ru_log <- read.csv(file.path(rbdir, "other_reports","read_until_decision_stats.csv"), header=TRUE, sep=";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
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
  
  ru_mapping_summary <- read.csv(file.path(rbdir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  
  ru_mapping_summary <- ru_mapping_summary[ru_mapping_summary$channel < 257,]
  
  sequencing_summary <- load_sequencing_summary(rbdir)
  
  mapped_fastq <- ru_log[c("read_id","decision")] %>%
    inner_join(ru_mapping_summary, by="read_id")
  
  matching_stats <- mapped_fastq %>%
    group_by(decision, contig) %>%
    count() %>%
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
      accuracy = (tp+tn)/(tp+tn+fp+fn),
      precision = tp/(tp+fp),
      recall = tp/(tp+fn),
      specificity = tn/(tn+fp),
      f1_score = 2*tp/(2*tp + fp + fn),
      balanced_acc = (recall + specificity)/2,
      #mcc = (as.double(tp)*as.double(tn) - fp*fn)/sqrt(as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)),
      tool = "ReadBouncer",
      reference = "all")
     
  
  # now for each species again
  
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
        accuracy = (tp+tn)/(tp+tn+fp+fn),
        precision = tp/(tp+fp),
        recall = tp/(tp+fn),
        specificity = tn/(tn+fp),
        f1_score = 2*tp/(2*tp + fp + fn),
        balanced_acc = (recall + specificity)/2,
        #mcc = (as.double(tp)*as.double(tn) - fp*fn)/sqrt(as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)),
        tool = "ReadBouncer",
        reference = r)
    
  }
  
  
  # now add second readbouncer directory metrics
  
  ru_log <- read.csv(file.path(rbdir2, "other_reports","read_until_decision_stats.csv"), header=TRUE, sep=";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
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
  
  ru_mapping_summary <- read.csv(file.path(rbdir2, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  
  ru_mapping_summary <- ru_mapping_summary[ru_mapping_summary$channel < 257,]
  
  sequencing_summary <- load_sequencing_summary(rbdir2)
  
  mapped_fastq <- ru_log[c("read_id","decision")] %>%
    inner_join(ru_mapping_summary, by="read_id")
  
  matching_stats <- mapped_fastq %>%
    group_by(decision, contig) %>%
    count() %>%
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
      accuracy = (tp+tn)/(tp+tn+fp+fn),
      precision = tp/(tp+fp),
      recall = tp/(tp+fn),
      specificity = tn/(tn+fp),
      f1_score = 2*tp/(2*tp + fp + fn),
      balanced_acc = (recall + specificity)/2,
      #mcc = (as.double(tp)*as.double(tn) - fp*fn)/sqrt(as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)),
      tool = "ReadBouncer2",
      reference = "all")
  
  
  # now for each species again
  
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
        accuracy = (tp+tn)/(tp+tn+fp+fn),
        precision = tp/(tp+fp),
        recall = tp/(tp+fn),
        specificity = tn/(tn+fp),
        f1_score = 2*tp/(2*tp + fp + fn),
        balanced_acc = (recall + specificity)/2,
        #mcc = (as.double(tp)*as.double(tn) - fp*fn)/sqrt(as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)),
        tool = "ReadBouncer2",
        reference = r)
    
  }
  
  
  sample_sheet <- read.csv(file.path(mkdir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  
  overall_metrics <- sample_sheet %>%
    select(Reference, Species) %>%
    right_join(overall_metrics, by=c("Reference" = "reference")) %>%
    mutate(Species = case_when(Reference == "all" ~ "all", TRUE ~ Species  )) %>%
    select(-"Reference")
 
  pdf(output_file, paper = "a4r", width = 20)
  grid.table(overall_metrics, rows = NULL)#, rows = NULL, cols = c("Species", "Flow Cell ID", "Mean Read Length", "Median Read Length", "Standard Deviation"))
  #grid.table(as, rows = NULL, cols = c("Species", "Flow Cell ID", "Mean Read Length", "Median Read Length", "Standard Deviation"))
  dev.off()
 
}

create_species_abundance_plot <- function(dirlist)
{
  #output_dir <- "D:/Plasmide/results"
  dirlist <- paste(rb_dirs, mk_dirs,sep = ",")
  mydirlist <- unlist(strsplit(dirlist, ","))
  abundance_report <- data.frame()
  for (dir in mydirlist)
  {
    #dir <- "D:/Plasmide/Data/20220321_1207_MN24598_FAR91003_ff83ee47"
    tmp_cov <- read.csv(file.path(dir, "results", "coverage_report.tsv"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    sample_sheet <- read.csv(file.path(dir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    
    dir_name <- SplitPath(dir)
    fc_name <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    exp_date <- unlist(strsplit(dir_name$filename,split = "_"))[1]
    yield <- tmp_cov %>%
      filter(as == "NO") %>%
      group_by(reference) %>%
      summarise(sum_bases = sum(bases)) %>%
      mutate(experiment = fc_name) %>%
      mutate(date = exp_date) %>%
      mutate(perc = sum_bases / sum(sum_bases)) %>%
      inner_join(sample_sheet, by=c("reference" = "Reference"))
    
    abundance_report <- bind_rows(abundance_report, yield)
  }
  
  abundance_report <- abundance_report %>%
    mutate(experiment = 
             case_when(experiment == "FAR91003" ~ "ReadBouncer1",
                       experiment == "FAR92672" ~ "ReadBouncer2",
                       experiment == "FAR92750" ~ "MinKNOW1",
                       experiment == "FAP84921" ~ "MinKNOW2",
                       TRUE ~ experiment)
          )
  
  sort <- abundance_report %>%
    select(experiment, date) %>%
    filter(!duplicated(.)) %>%
    arrange(date)
  
  sp_ab <- ggplot(abundance_report, aes(fill=Species, y=perc, x=factor(experiment, levels = sort$experiment))) + 
    geom_bar(position="fill", stat="identity")+
    scale_y_continuous(labels=scales::percent) +
    theme_minimal()+
    scale_fill_brewer(palette="YlGnBu")+
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    #theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    labs(title = "(b)", x="", y = "abundance")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    geom_text(aes(label = paste0(round(100 * perc, 2),"%")), 
              position = position_stack(vjust = 0.5), size = 4)
  
  return(sp_ab)
  
}

create_plasmid_abundance_plot <- function(rb_dirs, mk_dirs, output_file)
{
  
  mydirlist <- unlist(strsplit(rb_dirs, ","))
  abundance_report<- data.frame()
  for (dir in mydirlist)
  {
    #dir <- "D:/Plasmide/Data/20220324_1028_MN24598_FAR92672_9a263022"
    #tmp_cov <- read.csv(file.path(dir, "results", "coverage_report.tsv"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    mapping_summary <- read.csv(file.path(dir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "") %>%
      filter(contig != "none")
    sequencing_summary <- load_sequencing_summary(dir)
    
    merged <- sequencing_summary %>%
      filter(sequence_length_template >= 1000 & end_reason == "signal_positive" & passes_filtering == "TRUE") %>%
      select("read_id","channel","start_time","sequence_length_template") %>%
      inner_join(mapping_summary[c("read_id","reference","contig")], by="read_id") %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      filter(start_time < 24*60*60)
    
    dir_name <- SplitPath(dir)
    fc_name <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    exp_date <- unlist(strsplit(dir_name$filename,split = "_"))[1]
    
    sample_sheet <- read.csv(file.path(dir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    sample_sheet <- sample_sheet %>%
      mutate(Reference = case_when(Reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ Reference  ))
    
    yield <- merged %>%
      filter(channel > 256) %>%
      group_by(reference,contig) %>%
      summarise(sum_bases = sum(sequence_length_template)) %>%
      mutate(perc = sum_bases / sum(sum_bases)) %>%
      mutate(experiment = fc_name) %>%
      mutate(date = exp_date) %>%
      mutate(region = "Control") %>%
      mutate(tool = "ReadBouncer") %>%
      filter(contig == "plasmid" | contig == "chromosome") %>%
      #mutate(contig = case_when(plasmid == "YES" ~ "Plasmid", plasmid == "NO" ~ "Chromosome")) %>%
      inner_join(sample_sheet, by=c("reference" = "Reference"))
    
    abundance_report <- bind_rows(abundance_report, yield)
    
    yield <- merged %>%
      filter(channel < 257) %>%
      group_by(reference,contig) %>%
      summarise(sum_bases = sum(sequence_length_template)) %>%
      mutate(perc = sum_bases / sum(sum_bases)) %>%
      mutate(experiment = fc_name) %>%
      mutate(date = exp_date) %>%
      mutate(region = "Depletion") %>%
      mutate(tool = "ReadBouncer") %>%
      filter(contig == "plasmid" | contig == "chromosome") %>%
      #mutate(contig = case_when(plasmid == "YES" ~ "Plasmid", plasmid == "NO" ~ "Chromosome")) %>%
      inner_join(sample_sheet, by=c("reference" = "Reference"))
    
    abundance_report <- bind_rows(abundance_report, yield)
  }
  
  mydirlist <- unlist(strsplit(mk_dirs, ","))
  for (dir in mydirlist)
  {
    #dir <- "D:/Plasmide/Data/20220324_1028_MN24598_FAR92672_9a263022"
    #tmp_cov <- read.csv(file.path(dir, "results", "coverage_report.tsv"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    mapping_summary <- read.csv(file.path(dir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "") %>%
    filter(contig != "none")
    
    sequencing_summary <- load_sequencing_summary(dir)
    
    merged <- sequencing_summary %>%
      filter(sequence_length_template >= 1000 & end_reason == "signal_positive" & passes_filtering == "TRUE") %>%
      select("read_id","channel","start_time","sequence_length_template") %>%
      inner_join(mapping_summary[c("read_id","reference","contig")], by="read_id") %>%
      mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
      filter(start_time < 24*60*60)
    
    dir_name <- SplitPath(dir)
    fc_name <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    exp_date <- unlist(strsplit(dir_name$filename,split = "_"))[1]
    
    sample_sheet <- read.csv(file.path(dir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    sample_sheet <- sample_sheet %>%
      mutate(Reference = case_when(Reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ Reference  ))
    
    yield <- merged %>%
      filter(channel > 256) %>%
      group_by(reference,contig) %>%
      summarise(sum_bases = sum(sequence_length_template)) %>%
      mutate(perc = sum_bases / sum(sum_bases)) %>%
      mutate(experiment = fc_name) %>%
      mutate(date = exp_date) %>%
      mutate(region = "Control") %>%
      mutate(tool = "MinKNOW") %>%
      filter(contig == "plasmid" | contig == "chromosome") %>%
      #mutate(contig = case_when(plasmid == "YES" ~ "Plasmid", plasmid == "NO" ~ "Chromosome")) %>%
      inner_join(sample_sheet, by=c("reference" = "Reference"))
    
    abundance_report <- bind_rows(abundance_report, yield)
    
    yield <- merged %>%
      filter(channel < 257) %>%
      group_by(reference,contig) %>%
      summarise(sum_bases = sum(sequence_length_template)) %>%
      mutate(perc = sum_bases / sum(sum_bases)) %>%
      mutate(experiment = fc_name) %>%
      mutate(date = exp_date) %>%
      mutate(region = "Depletion") %>%
      mutate(tool = "MinKNOW") %>%
      filter(contig == "plasmid" | contig == "chromosome") %>%
      #mutate(contig = case_when(plasmid == "YES" ~ "Plasmid", plasmid == "NO" ~ "Chromosome")) %>%
      inner_join(sample_sheet, by=c("reference" = "Reference"))
    
    abundance_report <- bind_rows(abundance_report, yield)
  }
  
  abundance_report <- abundance_report %>%
    mutate(experiment = 
             case_when(experiment == "FAR91003" ~ "ReadBouncer1",
                       experiment == "FAR92672" ~ "ReadBouncer2",
                       experiment == "FAR92750" ~ "MinKNOW1",
                       experiment == "FAP84921" ~ "MinKNOW2",
                       TRUE ~ experiment)
    )
  
  sort <- abundance_report %>%
    select(experiment, date) %>%
    filter(!duplicated(.)) %>%
    arrange(date)
  
  pl_ab <- ggplot(abundance_report, aes(fill=contig, y=perc, x=Species)) + 
    geom_bar(position="fill", stat="identity")+
    scale_y_continuous(labels=scales::percent) +
    labs(x="", y = "")+
    scale_fill_brewer(palette="Blues")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.title = element_blank())+
    theme(legend.position = "top")+
    geom_text(aes(label = paste0(round(100 * perc, 2),"%")), 
              position = position_stack(vjust = 0.5), size = 4)+
    facet_nested(region ~ factor(experiment, levels = sort$experiment),scales = "free_x")+
    theme(strip.text.x = element_text(size = 10, face="bold"))+
    theme(strip.text.y = element_text(size = 10, face="bold"))+
    theme(strip.background=element_rect(color="grey30"))
  
  pdf(output_file, paper = "a4r", width = 20)
  print(pl_ab)
  dev.off()
}


print_active_channel_stats <- function(rb_dirs, mk_dirs, output_file)
{
  channel_report <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(channel_report) <- c("experiment","tool","timepoint","channels","design", "yield")
  channel_report$experiment <- as.character(channel_report$experiment)
  channel_report$tool <- as.character(channel_report$tool)
  channel_report$date <- as.character(channel_report$date)
  channel_report$timepoint <- as.numeric(channel_report$timepoint)
  channel_report$channels <- as.integer(channel_report$channels)
  channel_report$design <- as.character(channel_report$design)
  channel_report$yield <- as.integer(channel_report$yield)
  mydirlist <- unlist(strsplit(rb_dirs, ","))
  for (dir in mydirlist)
  {
    
    sequencing_summary <- load_sequencing_summary(dir)
    dir_name <- SplitPath(dir)
    fc_name <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    exp_date <- unlist(strsplit(dir_name$filename,split = "_"))[1]
    
    for (minutes in seq.int(from = 30, to = 24*60, by = 30))
    {
      
      start <- (minutes - 30) * 60
      end <- minutes * 60
      tmp_summary <- sequencing_summary %>%
        filter(template_start < end & template_start > start)
      
      
      first_half <- tmp_summary[tmp_summary$channel < 257,]
      second_half <- tmp_summary[tmp_summary$channel > 256,]
      
     channel_report <- channel_report %>%
       add_row(experiment = fc_name, tool = "ReadBouncer", date=exp_date, timepoint = minutes/60, channels = n_distinct(first_half$channel),design="Depletion", yield=sum(first_half$sequence_length_template)*2/1000000)%>%
       add_row(experiment = fc_name, tool = "ReadBouncer", date=exp_date, timepoint = minutes/60, channels = n_distinct(second_half$channel), design="Control", yield=sum(second_half$sequence_length_template)*2/1000000)
    }
  }
  
  mydirlist <- unlist(strsplit(mk_dirs, ","))
  for (dir in mydirlist)
  {
    
    sequencing_summary <- load_sequencing_summary(dir)
    dir_name <- SplitPath(dir)
    fc_name <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    exp_date <- unlist(strsplit(dir_name$filename,split = "_"))[1]
    
    for (minutes in seq.int(from = 30, to = 24*60, by = 30))
    {
     
      start <- (minutes - 30) * 60
      end <- minutes * 60
      tmp_summary <- sequencing_summary %>%
        filter(template_start < end & template_start > start)
      
      
      first_half <- tmp_summary[tmp_summary$channel < 257,]
      second_half <- tmp_summary[tmp_summary$channel > 256,]
      
      channel_report <- channel_report %>%
        add_row(experiment = fc_name, tool = "MinKNOW", date=exp_date, timepoint = minutes/60, channels = n_distinct(first_half$channel),design="Depletion", yield=sum(first_half$sequence_length_template)*2/1000000)%>%
        add_row(experiment = fc_name, tool = "MinKNOW", date=exp_date, timepoint = minutes/60, channels = n_distinct(second_half$channel), design="Control", yield=sum(second_half$sequence_length_template)*2/1000000)
    }
  }
  
  channel_report <- channel_report %>%
    mutate(experiment = 
             case_when(experiment == "FAR91003" ~ "ReadBouncer1",
                       experiment == "FAR92672" ~ "ReadBouncer2",
                       experiment == "FAR92750" ~ "MinKNOW1",
                       experiment == "FAP84921" ~ "MinKNOW2",
                       TRUE ~ experiment)
           )
   
  sort <- channel_report %>%
    select(experiment, date) %>%
    filter(!duplicated(.)) %>%
    arrange(date) 
  
  act <- ggplot(channel_report, aes(x=timepoint, y=channels, color=design)) +
    theme_minimal()+
    geom_line(aes(color=design), size=0.5)+
    xlim(c(0, 25))+
    ylim(c(0, 256))+
    #geom_point(aes(color=design))+
    scale_fill_jco()+
    #  scale_color_brewer(palette="Dark2")+
    labs(title ="(a)",x="Time (hours)", y = "Active Channels")+
    theme(legend.position = "top")+
    theme(legend.title = element_blank())+
    #theme(plot.margin=unit(c(0,0,4,0),"cm"))+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    theme(legend.title = element_blank())+
    facet_nested(~factor(experiment, levels = sort$experiment))+
    theme(strip.text.x = element_text(size = 10, face="bold"))+
    theme(strip.text.y = element_text(size = 10, face="bold"))+
    theme(strip.background=element_rect(color="grey30"))
  
  yph <- ggplot(channel_report, aes(x=timepoint, y=yield, color=design)) +
    theme_minimal()+
    geom_line(aes(color=design), size=0.5)+
    xlim(c(0, 25))+
    #ylim(c(0, 256))+
    #geom_point(aes(color=design))+
    scale_fill_jco()+
    #  scale_color_brewer(palette="Dark2")+
    labs(title ="(b)",x="Time (hours)", y = "Mbases per hour")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    theme(legend.position = "top")+
    theme(legend.title = element_blank())+
    facet_nested(~ factor(experiment, levels = sort$experiment))+
    theme(strip.text.x = element_text(size = 10, face="bold"))+
    theme(strip.text.y = element_text(size = 10, face="bold"))+
    theme(strip.background=element_rect(color="grey30"))
  
  combined <- grid.arrange(act, yph, nrow = 2)
  
  #combined
  
  ggsave(filename=output_file, 
         plot = combined, 
         device = cairo_pdf, 
         width = 297, 
         height = 210, 
         units = "mm")
  
  dev.off()
}

combine_unblocked_histograms <- function(rb_dirs, mk_dirs, output_file)
{
  rbdirlist <- unlist(strsplit(rb_dirs, ","))
  all_unblocked_reads <- data.frame()
  for (dir in rbdirlist)
  {
    tmp_summary <- load_sequencing_summary(dir)
    unblocked_reads <- tmp_summary %>%
      filter(channel < 257 & end_reason == "data_service_unblock_mux_change") %>%
      select(read_id, sequence_length_template) %>%
      mutate(tool = "ReadBouncer")
    
    all_unblocked_reads <- bind_rows(all_unblocked_reads, unblocked_reads)
  }
  
  mkdirlist <- unlist(strsplit(mk_dirs, ","))
  for (dir in mkdirlist)
  {
    tmp_summary <- load_sequencing_summary(dir)
    unblocked_reads <- tmp_summary %>%
      filter(channel < 257 & end_reason == "data_service_unblock_mux_change") %>%
      select(read_id, sequence_length_template) %>%
      mutate(tool = "MinKNOW")
    
    all_unblocked_reads <- bind_rows(all_unblocked_reads, unblocked_reads)
  }
  
  
  ubp <- ggplot(all_unblocked_reads, aes(x=sequence_length_template, fill=tool)) +
    geom_histogram(position="identity", binwidth = 5, alpha=0.5) +
    xlim(c(0, 3000))+
    ylim(c(0, 6500))+
    theme_minimal()+
    labs(x="read length", y = "count")+
    theme(legend.title = element_blank(), legend.position = "top")
  
  pdf(output_file, paper = "a4r", width = 20)
  print(ubp)
  dev.off()
}

print_read_length_info <- function(dirlist, output_file)
{
  
  mydirlist <- unlist(strsplit(dirlist, ","))
  control <- data.frame()
  as <- data.frame()
  
  for (dir in mydirlist)
  {
    #dir <- "D:/Plasmide/Data/20220324_1028_MN24598_FAR92672_9a263022"
    sequencing_summary <- load_sequencing_summary(dir)
    sample_sheet <- read.csv(file.path(dir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    
    dir_name <- SplitPath(dir)
    exp_data <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    
    df_control <- select(sample_sheet,-Plasmids,-Chromosome) %>%
      mutate(flow_cell = exp_data) %>%
      mutate(mean_rd_length = 0) %>%
      mutate(median_rd_length = 0) %>%
      mutate(sd_rd_length = 0)
    
    df_as <- select(sample_sheet,-Plasmids,-Chromosome) %>%
      mutate(flow_cell = exp_data) %>%
      mutate(mean_rd_length = 0) %>%
      mutate(median_rd_length = 0) %>%
      mutate(sd_rd_length = 0)
  
    reduced_second <- sequencing_summary %>%
      #filter(channel > 256) %>%
      filter(end_reason == 'signal_positive' & passes_filtering == 'TRUE' & channel > 256) %>%
      select("read_id","barcode_full_arrangement", "sequence_length_template")
    
    reduced_first <- sequencing_summary %>%
      #filter(channel < 257) %>%
      filter(end_reason == 'signal_positive' & passes_filtering == 'TRUE' & channel < 257) %>%
      select("read_id","barcode_full_arrangement", "sequence_length_template") %>%
      mutate(reference = NA)
  
    for (row in 1:nrow(sample_sheet))
    {
      #row <- 1
      barcode <- toString(sample_sheet[row, "Barcode"])
      rbk_second_half <- reduced_second[reduced_second$barcode_full_arrangement == barcode,]
      rbk_first_half <- reduced_first[reduced_first$barcode_full_arrangement == barcode,]
      
      df_control <- df_control %>%
        mutate(mean_rd_length = replace(mean_rd_length, Barcode == barcode, mean(rbk_second_half$sequence_length_template, na.rm = TRUE))) %>%
        mutate(median_rd_length = replace(median_rd_length, Barcode == barcode, median(rbk_second_half$sequence_length_template, na.rm = TRUE))) %>%
        mutate(sd_rd_length = replace(sd_rd_length, Barcode == barcode, sd(rbk_second_half$sequence_length_template, na.rm = TRUE)))
      
      df_as <- df_as %>%
        mutate(mean_rd_length = replace(mean_rd_length, Barcode == barcode, mean(rbk_first_half$sequence_length_template, na.rm = TRUE))) %>%
        mutate(median_rd_length = replace(median_rd_length, Barcode == barcode, median(rbk_first_half$sequence_length_template, na.rm = TRUE))) %>%
        mutate(sd_rd_length = replace(sd_rd_length, Barcode == barcode, sd(rbk_first_half$sequence_length_template, na.rm = TRUE)))
    }
    
    control <- bind_rows(control, df_control)
    as <- bind_rows(as, df_as)
  }
  
  control <- control %>%
    select(-"Barcode", -"Reference")
  
  as <- as %>%
    select(-"Barcode", -"Reference")
  
 
  pdf(output_file,paper = "a4r", width = 20)
  grid.table(control)
  dev.off()
}

create_read_length_histogram_plots <- function(dirlist)
{
  mydirlist <- unlist(strsplit(dirlist, ","))
  control <- data.frame()
 
  for (dir in mydirlist)
  {
    #dir <- "D:/Plasmide/Data/20220324_1028_MN24598_FAR92672_9a263022"
    sequencing_summary <- load_sequencing_summary(dir)
    sample_sheet <- read.csv(file.path(dir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
    
    dir_name <- SplitPath(dir)
    fc_name <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    exp_date <- unlist(strsplit(dir_name$filename,split = "_"))[1]
    
    reduced_control <- sequencing_summary %>%
      inner_join(sample_sheet,by=c("barcode_full_arrangement" = "Barcode")) %>%
      filter(end_reason == 'signal_positive' & passes_filtering == 'TRUE' & channel > 256) %>%
      select("read_id","Species", "sequence_length_template") %>%
      mutate(experiment = fc_name) %>%
      mutate(date = exp_date) 
    
    control <- bind_rows(control, reduced_control)
  }
  
  control <- control %>%
    mutate(experiment = 
             case_when(experiment == "FAR91003" ~ "ReadBouncer1",
                       experiment == "FAR92672" ~ "ReadBouncer2",
                       experiment == "FAR92750" ~ "MinKNOW1",
                       experiment == "FAP84921" ~ "MinKNOW2",
                       TRUE ~ experiment)
          )
  
  sort <- control %>%
    select(experiment, date) %>%
    filter(!duplicated(.)) %>%
    arrange(date)
  
  histo <- ggplot(control, aes(x=sequence_length_template, fill=Species)) +
    #geom_histogram(position="identity", binwidth = 100, alpha=0.25)+
    geom_density(aes(y = ..count..),position="identity", alpha=0.5)+
    scale_fill_brewer(palette="YlGnBu")+
    #geom_freqpoly(position="identity", binwidth = 100)+
    xlim(c(0, 50000))+
    #ylim(c(0, 3500))+
    theme_minimal()+
    theme(legend.position="bottom")+
    theme(legend.direction='vertical')+
    guides(fill = guide_legend(ncol = 2, bycol = TRUE))+
    theme(legend.title = element_blank())+
    #theme(plot.margin=unit(c(0,0,4,0),"cm"))+
    labs(title="(c)", x="Read Length", y = "count")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    facet_wrap(~factor(experiment, levels = sort$experiment), nrow = 4,scales = "free_y")
  
  return(histo)
}

create_read_length_violin_plot <- function(dirlist)
{
  mydirlist <- unlist(strsplit(dirlist, ","))
  control <- data.frame()
  
  for (dir in mydirlist)
  {
    
    sequencing_summary <- load_sequencing_summary(dir)
    
    dir_name <- SplitPath(dir)
    fc_name <- unlist(strsplit(dir_name$filename,split = "_"))[4]
    exp_date <- unlist(strsplit(dir_name$filename,split = "_"))[1]
    
    reduced_control <- sequencing_summary %>%
      filter(end_reason == 'signal_positive' & passes_filtering == 'TRUE' & channel > 256) %>%
      select("read_id", "sequence_length_template") %>%
      mutate(experiment = fc_name) %>%
      mutate(date = exp_date) 
    
    control <- bind_rows(control, reduced_control)
  }
  
  control <- control %>%
    mutate(experiment = 
             case_when(experiment == "FAR91003" ~ "ReadBouncer1",
                       experiment == "FAR92672" ~ "ReadBouncer2",
                       experiment == "FAR92750" ~ "MinKNOW1",
                       experiment == "FAP84921" ~ "MinKNOW2",
                       TRUE ~ experiment)
    )
  
  sort <- control %>%
    select(experiment, date) %>%
    filter(!duplicated(.)) %>%
    arrange(date)
  
  violin <- ggplot(control, aes(y=sequence_length_template, x=factor(experiment, levels = sort$experiment), fill=experiment)) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white")+
    scale_fill_brewer(palette="Accent")+
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
    theme_minimal()+
    theme(legend.position="none")+
    labs(title="(a)", y="Read Length", x = "")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot") #NEW parameter. Apply for subtitle too.
     
  
  return(violin)
}


combine_rb2_exp_plots <- function(rb2_dir, output_file)
{
  
  # create active channels plot
  
  channel_report <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(channel_report) <- c("timepoint","channels","design", "yield")
  channel_report$timepoint <- as.numeric(channel_report$timepoint)
  channel_report$channels <- as.integer(channel_report$channels)
  channel_report$design <- as.character(channel_report$design)
  channel_report$yield <- as.integer(channel_report$yield)
 
  sequencing_summary <- load_sequencing_summary(rb2_dir)

  for (minutes in seq.int(from = 30, to = 24*60, by = 30))
  {
    start <- (minutes - 30) * 60
    end <- minutes * 60
    tmp_summary <- sequencing_summary %>%
      filter(template_start < end & template_start > start)
    
      
    first_half <- tmp_summary[tmp_summary$channel < 257,]
    second_half <- tmp_summary[tmp_summary$channel > 256,]
      
    channel_report <- channel_report %>%
      add_row(timepoint = minutes/60, channels = n_distinct(first_half$channel),design="Depletion", yield=sum(first_half$sequence_length_template)*2/1000000)%>%
      add_row(timepoint = minutes/60, channels = n_distinct(second_half$channel), design="Control", yield=sum(second_half$sequence_length_template)*2/1000000)
  }
  
  act <- ggplot(channel_report, aes(x=timepoint, y=channels, color=design)) +
    theme_minimal()+
    geom_line(aes(color=design), size=0.5)+
    xlim(c(0, 25))+
    ylim(c(0, 256))+
    scale_fill_jco()+
    labs(title ="(a)",x="Time (hours)", y = "Active Channels")+
    theme(legend.position = "top")+
    theme(legend.title = element_blank())+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    theme(legend.title = element_blank())+
    theme(strip.text.x = element_text(size = 10, face="bold"))+
    theme(strip.text.y = element_text(size = 10, face="bold"))+
    theme(strip.background=element_rect(color="grey30"))
  
  yph <- ggplot(channel_report, aes(x=timepoint, y=yield, color=design)) +
    theme_minimal()+
    geom_line(aes(color=design), size=0.5)+
    xlim(c(0, 25))+
    ylim(c(0, 180))+
    scale_fill_jco()+
    labs(title ="(b)",x="Time (hours)", y = "Mbases per hour")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    theme(legend.position = "top")+
    theme(legend.title = element_blank())+
    theme(strip.text.x = element_text(size = 10, face="bold"))+
    theme(strip.text.y = element_text(size = 10, face="bold"))+
    theme(strip.background=element_rect(color="grey30"))
  
  
  # create plasmid abundance 
  
  abundance_report<- data.frame()
  mapping_summary <- read.csv(file.path(rb2_dir, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "") %>%
    filter(contig != "none")
    
  merged <- sequencing_summary %>%
    filter(sequence_length_template >= 1000 & end_reason == "signal_positive" & passes_filtering == "TRUE") %>%
    select("read_id","channel","start_time","sequence_length_template") %>%
    inner_join(mapping_summary[c("read_id","reference","contig")], by="read_id") %>%
    mutate(reference = case_when(reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ reference  )) %>%
    filter(start_time < 24*60*60)
    
  sample_sheet <- read.csv(file.path(rb2_dir, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  sample_sheet <- sample_sheet %>%
    mutate(Reference = case_when(Reference == "10-02932_Mallorca" ~ "10-02932", TRUE ~ Reference  ))
    
  yield <- merged %>%
    filter(channel > 256) %>%
    group_by(reference,contig) %>%
    summarise(sum_bases = sum(sequence_length_template)) %>%
    mutate(perc = sum_bases / sum(sum_bases)) %>%
    mutate(region = "Control") %>%
    filter(contig == "plasmid" | contig == "chromosome") %>%
    inner_join(sample_sheet, by=c("reference" = "Reference"))
    
  abundance_report <- bind_rows(abundance_report, yield)
  
 
  yield <- merged %>%
    filter(channel < 257) %>%
    group_by(reference,contig) %>%
    summarise(sum_bases = sum(sequence_length_template)) %>%
    mutate(perc = sum_bases / sum(sum_bases)) %>%
    mutate(region = "Depletion") %>%
    filter(contig == "plasmid" | contig == "chromosome") %>%
    inner_join(sample_sheet, by=c("reference" = "Reference"))
    
  abundance_report <- bind_rows(abundance_report, yield)
  
  pl_ab <- ggplot(abundance_report, aes(fill=contig, y=perc, x=Species)) + 
    geom_bar(position="fill", stat="identity")+
    scale_y_continuous(labels=scales::percent) +
    labs(title="(d)", x="", y = "")+
    scale_fill_brewer(palette="Blues")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.title = element_blank())+
    theme(legend.position = "top")+
    geom_text(aes(label = paste0(round(100 * perc, 2),"%")), 
              position = position_stack(vjust = 0.5), size = 4)+
    facet_grid(region ~ .,scales = "free_x")+
    theme(strip.text.x = element_text(size = 10, face="bold"))+
    theme(strip.text.y = element_text(size = 10, face="bold"))+
    theme(strip.background=element_rect(color="grey30"))
  

  # create enrichment plot
  
  ef <- load_coverage_report(rb2_dir, 1000)
  
  ef_yield <- 
    ggplot(ef, aes(x=timepoint, y=efactor_bases, group=Species)) +
    theme_minimal()+
    geom_line(aes(color=Species), size=1.0)+
    xlim(c(0, 25))+
    ylim(c(0, 2.5))+
    scale_fill_jco()+
    scale_color_brewer(palette="Set2")+
    labs(title ="(c)",x="Time (hours)", y = "Enrichment by yield")+
    theme(plot.title = element_text(size=16))+
    theme(plot.title.position = "plot")+
    theme(legend.position = "bottom")+
    theme(legend.title = element_blank())+
    theme(legend.direction='horizontal')+
    guides(color = guide_legend(ncol = 2, bycol=TRUE))+
    theme(strip.text.x = element_text(size = 8, face="bold"))+
    theme(strip.text.y = element_text(size = 8, face="bold"))+
    theme(strip.background=element_rect(color="grey30"))
  
  #create unblocked read length histo
  
  unblocked_reads <- sequencing_summary %>%
    filter(channel < 257 & end_reason == "data_service_unblock_mux_change") %>%
    select(read_id, sequence_length_template)
  
  ggplot(unblocked_reads, aes(x=sequence_length_template, fill="sequence_length_template")) +
    geom_histogram(aes(y=..density..), binwidth = 5, alpha=0.5) +
    xlim(c(0, 3000))+
    theme_minimal()+
    theme(axis.text=element_text(size=14))+
    labs(title ="", x="rejected read length", y = "")+
    theme(axis.title=element_text(size=16))+
    theme(legend.position ="none")+
    stat_function(fun=dnbinom, geom="line", size=1, color="black",
                  args = list(mu=650, size=8.6))
  
  ubp <- ggplot(unblocked_reads, aes(x=sequence_length_template, fill="sequence_length_template")) +
    geom_histogram(position="identity", binwidth = 5, alpha=0.5) +
    xlim(c(0, 3000))+
    ylim(c(0, 2000))+
    theme_minimal()+
    labs(title ="(e)", x="rejected read lengths", y = "count")+
    theme(legend.position ="none")
  
  
  
  rb2_plot <- grid.arrange(
    grobs = list(act, yph, pl_ab, ef_yield, ubp),
    widths = c(1, 1),
    layout_matrix = rbind(c(1, 3),c(2, 3), c(4,5))
  )
  
  ggsave(filename=output_file, 
         plot = rb2_plot, 
         device = cairo_pdf, 
         width = 210, 
         height = 297, 
         units = "mm")
  
  dev.off()
  
}


create_assembly_stats <- function(rbdir1, output_file)
{
  
  rb_stats <- read.csv(file.path(rbdir1, "assembly_quality_stats.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  sample_sheet <- read.csv(file.path(rbdir1, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  
  
  rb_stats <- rb_stats %>%
    filter(timesplit == 60 | timesplit == 120 | timesplit == 180) %>%
    mutate(time = timesplit/60) %>%
    inner_join(sample_sheet, by=c("reference" = "Reference"))%>%
    rename(depth = avg_coverage_depth) %>%
    rename(reads = reads_mapped)%>%
    select(-timesplit, -ref_cov1, -reference, -Barcode, -Plasmids, -Chromosome)%>%
    relocate(time)%>%
    relocate(Species)%>%
    arrange(Species, channel, time)
    
  
  pdf(output_file,paper = "a4r", width = 20)
  grid.table(rb_stats)
  dev.off()
}


option_list = list(
  make_option(c("", "--data_dir"), type="character", default=NULL, 
              help="directory containing the five nanopore run directories"),
  make_option(c("", "--min_read_length"), type = "integer", default=1000, 
              help="minimum read length of successfully sequenced reads"),
  make_option(c("", "--output_dir"), type="character", default=NULL, 
              help="output directory for plots")
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



if (is.null(opt$rb_dirs) || is.null(opt$mk_dirs)){
  print_help(opt_parser)
  stop("Coverage report files missing.", call.=FALSE)
}

opt$output_dir <- "D:/Plasmide/results"
opt$min_read_length <- 1000
opt$data_dir <- "D:/Plasmide/Data"


mkdir1 <- file.path(opt$data_dir,"20220419_1137_MN24598_FAR92750_43bb5c8c")
mkdir2 <- file.path(opt$data_dir,"20220420_1351_MN24598_FAP84921_7d466390")
rbdir1 <- file.path(opt$data_dir,"20220321_1207_MN24598_FAR91003_ff83ee47")
rbdir2 <- file.path(opt$data_dir,"20220324_1028_MN24598_FAR92672_9a263022")
rbdir3 <- file.path(opt$data_dir,"20220608_1124_MN24598_FAT03513_da85ca38")

rb_dirs <- paste(rbdir1, rbdir2, sep = ",")
mk_dirs <- paste(mkdir1, mkdir2, sep = ",")

dirlist <- paste(rb_dirs, mk_dirs,sep = ",")

rb_cov <- load_coverage_report(rb_dirs, opt$min_read_length)
mk_cov <- load_coverage_report(mk_dirs, opt$min_read_length)

# Figure 4
create_ef_plots(rb_cov, mk_cov, file.path(opt$output_dir, "EF_experiment_comparison.pdf"))

# Figure 2
print_active_channel_stats(rb_dirs, mk_dirs, file.path(opt$output_dir, "active_channels_comparison.pdf"))

# creates data for Table 2 & Suppl. Tabl 4
create_matching_stats(rbdir2, mkdir2, rbdir3, file.path(opt$output_dir, "ru_matching_stats.pdf") )

# creates Suppl. Table 3
print_read_length_info(paste(rb_dirs, mk_dirs,sep = ","), file.path(opt$output_dir, "read_length_info.pdf"))

# the following lines create Figure 1
sp_ab <- create_species_abundance_plot(paste(rb_dirs, mk_dirs,sep = ","))
histos <- create_read_length_histogram_plots(paste(rb_dirs, mk_dirs,sep = ","))
violins <- create_read_length_violin_plot(paste(rb_dirs, mk_dirs,sep = ","))

rl_spab <- grid.arrange(
  grobs = list(violins, sp_ab, histos),
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 3),c(2, 3))
)

ggsave(filename=file.path(opt$output_dir, "Figure1.pdf"), 
       plot = rl_spab, 
       device = cairo_pdf, 
       width = 297, 
       height = 210, 
       units = "mm")

dev.off()

# creates Figure 3
create_plasmid_abundance_plot(rb_dirs, mk_dirs, file.path(opt$output_dir, "plasmid_abundance_plots.pdf"))

# creates Suppl. Figure 5
combine_unblocked_histograms(rb_dirs, mk_dirs, file.path(opt$output_dir, "unblocked_histogram.pdf"))

# creates Suppl. Figure 6
combine_rb2_exp_plots(rbdir3, file.path(opt$output_dir, "Suppl_Figure6.pdf"))

#creates Table 2
create_assembly_stats(rbdir2,file.path(opt$output_dir, "assembly_stats.pdf"))
