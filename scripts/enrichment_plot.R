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

create_coverage_report <- function(read_stats)
{
  #read_stats <- merged
  hours <- 24
  interval <- 30
  cov_report <- data.frame()
  for (minutes in seq.int(from = interval, to = hours*60, by = interval))
  {
    #minutes <- 30
    #start <- (minutes - interval) * 60
    end <- minutes * 60
    tmp_summary <- read_stats[read_stats$start_time < end,]
    #tmp_summary <- tmp_summary[tmp_summary$start_time >= start,]
    
    first_half <- tmp_summary[tmp_summary$channel < 257,]
    second_half <- tmp_summary[tmp_summary$channel > 256,]
    
    readnum_first <- first_half %>%
      group_by(reference, contig) %>%
      count() %>%
      mutate(as = "YES", timepoint = minutes) %>%
      rename(reads = n)
    
    bases_first <- first_half %>%
      group_by(reference, contig) %>%
      summarise(bases = sum(sequence_length_template)) %>%
      mutate(as = "YES", timepoint = minutes)
    
    first <- readnum_first %>%
      inner_join(bases_first)
    
    readnum_second <- second_half %>%
      group_by(reference, contig) %>%
      count() %>%
      mutate(as = "NO", timepoint = minutes) %>%
      rename(reads = n)
    
    bases_second <- second_half %>%
      group_by(reference, contig) %>%
      summarise(bases = sum(sequence_length_template)) %>%
      mutate(as = "NO", timepoint = minutes)
    
    second <- readnum_second %>%
      inner_join(bases_second)
    
    cov_report <- bind_rows(cov_report, first, second)
  }
  return(cov_report)
}

plot_mean_depth <- function(data_folder, output_file)
{
  #data_folder <- "D:/Plasmide/Data/20220419_1137_MN24598_FAR92750_43bb5c8c"
  
  depth_summary <- read.csv(file.path(data_folder, "plasmid_depth_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

  changed <- depth_summary %>%
    mutate(reference = case_when(channel == '1-256' ~ paste(reference,"(AS)"), channel == '257-512' ~ paste(reference,"(Control)")))
  
  rp <- ggplot(changed, aes(x=timesplit, y=meandepth, group=reference)) +
    theme_minimal()+
    geom_line(aes(color=reference))+
    xlim(c(0, 16*60))+
    ylim(c(0, 5000))+
    geom_point(aes(color=reference))+
    scale_fill_jco()+
    #  scale_color_brewer(palette="Dark2")+
    labs(title="Mean Depth of Plasmid References over time", x="Time in minutes", y = "Depth")+
    theme(legend.position = "top")+
    theme(legend.title = element_blank())
  
  pdf(output_file)
  print(rp)
  dev.off()
  
}

option_list = list(
  make_option(c("", "--data_dir"), type="character", default=NULL, 
              help="Nanopore run directory"),
  make_option(c("", "--min_read_length"), type = "integer", default=1000, 
              help="minimum read length of successfully sequenced reads"),
  make_option(c("", "--output_dir"), type="character", default=NULL, 
              help="output directory")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$data_dir) || is.null(opt$output_dir))
{
  print_help(opt_parser)
  stop("Data directory and output directory have to be provided.n", call.=FALSE)
}

data_folder <- opt$data_dir
#data_folder <- "D:/Plasmide/Data/20220419_1137_MN24598_FAR92750_43bb5c8c"
#opt$min_read_length <- 1500

sequencing_summary <- load_sequencing_summary(data_folder)

sequencing_summary <- sequencing_summary[sequencing_summary$end_reason == 'signal_positive',]
sequencing_summary <- sequencing_summary[sequencing_summary$passes_filtering == 'TRUE',]

mapping_summary <- read.csv(file.path(data_folder, "mapping_summary.txt"), header=TRUE, sep="\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

merged <- sequencing_summary[c("read_id","channel","start_time","sequence_length_template")] %>%
  inner_join(mapping_summary[c("read_id","reference","contig")], by="read_id")

coverage_report <- create_coverage_report(merged)

write.table(coverage_report, file = file.path(data_folder, "results", "coverage_report.tsv" ), sep = "\t", 
            col.names = TRUE, row.names = FALSE)






