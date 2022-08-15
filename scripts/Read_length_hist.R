if(!require(optparse)) install.packages("optparse")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(scales)) install.packages("scales")
if(!require(ggsci)) install.packages("ggsci")
if(!require(gridExtra)) install.packages("gridExtra")

library(optparse)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(gridExtra)


# function definitions

create_timesplit_tsv <- function(sample_sheet, sequencing_summary, data_folder, min_rd_len)
{

  csv_folder <- file.path(data_folder, "csv")
  if (!dir.exists(csv_folder))
  {
      dir.create(csv_folder)
  }
  
  first_half <- sequencing_summary[sequencing_summary$channel < 257,]
  second_half <- sequencing_summary[sequencing_summary$channel > 256,]
  
    
  first_half <- first_half[first_half$sequence_length_template >= min_rd_len, ]
  second_half <- second_half[second_half$sequence_length_template >= min_rd_len,]
    
    
    
  for (row in 1:nrow(sample_sheet))
  {
    barcode <- sample_sheet[row, "Barcode"]
    rbk_first_half <- first_half[first_half$barcode_full_arrangement == barcode,]
    rbk_second_half <- second_half[second_half$barcode_full_arrangement == barcode,]
      
    write.table(rbk_first_half, file = file.path(csv_folder, paste(barcode,"_pass_1-256.tsv", sep="")), sep = "\t", 
                  col.names = TRUE, row.names = FALSE)
    write.table(rbk_second_half, file = file.path(csv_folder,paste(barcode,"_pass_257-512.tsv", sep="")), sep = "\t", 
                  col.names = TRUE, row.names = FALSE)
  }

  csv_folder <- file.path(data_folder, "csv","timesplit")
    
    if (!dir.exists(csv_folder))
    {
      dir.create(csv_folder)
    }
    
  for (minutes in seq.int(from = 16*60 + 30, to = 24*60, by = 30))
  {
    
    csv_folder <- file.path(data_folder, "csv","timesplit",minutes)
    
    if (!dir.exists(csv_folder))
    {
      dir.create(csv_folder)
    }
    
    seconds <- minutes * 60
    tmp_summary <- sequencing_summary[sequencing_summary$template_start <= seconds,]
    first_half <- tmp_summary[tmp_summary$channel < 257,]
    second_half <- tmp_summary[tmp_summary$channel > 256,]
    
    
    first_half <- first_half[first_half$sequence_length_template >= min_rd_len, ]
    second_half <- second_half[second_half$sequence_length_template >= min_rd_len,]
    
    for (row in 1:nrow(sample_sheet))
    {
      barcode <- sample_sheet[row, "Barcode"]
      rbk_first_half <- first_half[first_half$barcode_full_arrangement == barcode,]
      rbk_second_half <- second_half[second_half$barcode_full_arrangement == barcode,]
      
      write.table(rbk_first_half, file = file.path(csv_folder, paste(barcode,"_pass_1-256.tsv", sep="")), sep = "\t", 
                  col.names = TRUE, row.names = FALSE)
      write.table(rbk_second_half, file = file.path(csv_folder,paste(barcode,"_pass_257-512.tsv", sep="")), sep = "\t", 
                  col.names = TRUE, row.names = FALSE)
    }
    
  }
}

plot_unblocked_histogram <- function(sequencing_summary, output_file)
{
  rb_unblocked <- sequencing_summary[sequencing_summary$end_reason == 'data_service_unblock_mux_change',]
  
  rb<-ggplot(rb_unblocked, aes(x=sequence_length_template)) +
    geom_histogram(binwidth = 5, fill="#0072B2") +
    xlim(c(0, 3000))+
    ylim(c(0, 3500))+
    theme_minimal()+
    ggtitle("Unblocked Read Lengths")+
    labs(x="read length", y = "count")
  
  #rb
 # unblocked_plot <- file.path(opt$output_dir, "Unblocked_histogram.pdf")
  pdf(output_file)
  print(rb)
  dev.off()
  
}

plot_channel_information <- function(sequencing_summary, output_file)
{
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  x <- c("name", "timepoint", "channels","total_yield", "yield_per_channel")
  colnames(df) <- x
    
  for (minutes in seq.int(from = 60, to = 16*60, by = 60))
  {
    start <- (minutes - 60) * 60
    end <- minutes * 60
    tmp_summary <- rb_sequencing_summary[rb_sequencing_summary$template_start < end,]
    tmp_summary <- tmp_summary[tmp_summary$template_start >= start,]
    
    first_half <- tmp_summary[tmp_summary$channel < 257,]
    second_half <- tmp_summary[tmp_summary$channel > 256,]
    active_first_half <- as.numeric(n_distinct(first_half$channel))
    active_second_half <- as.numeric(n_distinct(second_half$channel))
    
    first_half <- first_half[first_half$end_reason == "signal_positive",]
    #first_half <- first_half[first_half$sequence_length_template >= 1000,]
    first_half <- first_half[first_half$passes_filtering == 'TRUE',]
    first_half_sum_bases <- sum(first_half$sequence_length_template)
    
    second_half <- second_half[second_half$end_reason == "signal_positive",]
    #second_half <- second_half[second_half$sequence_length_template >= 1000,]
    second_half <- second_half[second_half$passes_filtering == 'TRUE',]
    second_half_sum_bases <- sum(second_half$sequence_length_template)
    
    new_row = c("1-256", as.numeric(minutes), active_first_half,first_half_sum_bases, first_half_sum_bases/active_first_half)
    df[nrow(df) + 1,] <- new_row
    new_row = c("257-512", as.numeric(minutes), active_second_half,second_half_sum_bases, second_half_sum_bases/active_second_half)
    df[nrow(df) + 1,] <- new_row
  }
  
  df$timepoint <- as.numeric(as.character(df$timepoint))
  df$channels <- as.numeric(as.character(df$channels))
  df$total_yield <- as.numeric(as.character(df$total_yield))
  df$yield_per_channel <- as.numeric(as.character(df$yield_per_channel))
  
  
  cp <- ggplot(df, aes(x=timepoint, y=channels, group=name)) +
    theme_minimal()+
    geom_line(aes(color=name))+
    xlim(c(0, 16*60 + 60))+
    ylim(c(0, 256))+
    geom_point(aes(color=name))+
    scale_fill_jco()+
    #  scale_color_brewer(palette="Dark2")+
    labs(title="Active Sequencing Channels", x="Time in minutes", y = "Number of Channels")+
    theme(legend.position = "top")+
    theme(legend.title = element_blank())
  
  #cp
  #cp_plot <- file.path(opt$output_dir, "active_channels.pdf")
  pdf(output_file)
  print(cp)
  dev.off()
  
}


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


create_first_half_csv <- function(sample_sheet, sequencing_summary, data_folder)
{
  csv_folder <- file.path(data_folder, "csv")
  if (!dir.exists(csv_folder))
  {
      dir.create(csv_folder)
  }
  
  first_half <- sequencing_summary[sequencing_summary$channel < 257,]
    
  for (row in 1:nrow(sample_sheet))
  {
    barcode <- sample_sheet[row, "Barcode"]
    rbk_first_half <- first_half[first_half$barcode_full_arrangement == barcode,]
      
    write.table(rbk_first_half, file = file.path(csv_folder, paste(barcode,"_all_1-256.tsv", sep="")), sep = "\t", 
                  col.names = TRUE, row.names = FALSE)
  }
}

# option parsing

option_list = list(
  make_option(c("", "--data_dir"), type="character", default=NULL, 
              help="Nanopore run directory"),
  make_option(c("", "--min_read_length"), type = "integer", default=1000, 
              help="minimum read length of successfully sequenced reads")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$data_dir)){
  print_help(opt_parser)
  stop("Data directory has to be provided.", call.=FALSE)
}


data_folder <- opt$data_dir
# <- "D:/Plasmide/Data/20220321_1207_MN24598_FAR91003_ff83ee47"

sample_sheet <- read.csv(file.path(data_folder, "sample_sheet.csv"), header=TRUE, sep=",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

# load and combine data from sequencing summary files
rb_sequencing_summary <- load_sequencing_summary(data_folder)


create_first_half_csv(sample_sheet, rb_sequencing_summary, data_folder)

signal_positive <- rb_sequencing_summary[rb_sequencing_summary$end_reason == 'signal_positive',]
signal_positive <- signal_positive[signal_positive$passes_filtering == 'TRUE',]
# create tsv files for 30 minute intervals
create_timesplit_tsv(sample_sheet, signal_positive, data_folder, opt$min_read_length)
# plot unblocked read length histogram of adaptive sequencing channels

