# An R-script to process the output from Kraken

For one of my projects I analyzed the output from the Kraken classification for a set of metagenomic sequences.

 The kraken output files were generated using the command kraken_report. I then wanted to combine the files in R, so that I could start comparing the samples. Below you find the script that I made.
 
 You find the script file here: [combine kraken output files.R](R-script/combine_kraken_output_files.R)
 
 This is the script:
 
 ```
 # a script to analyze the kraken output.
# This script should be stored upstream of a data folder called: "data"
# It will make a folder called "data_edit", where it will write the reformated Kraken_output
# It then reads the reformated tables into one object (tibble) called data,
# this table is in long format
# which then is used for plotting the total number of reads of each dataset.

# required libraries
library(tidyverse) #loads dplyr, ggplot, read, tidyr, do not load plyr!


### FUNCTIONS ####

##function_1 ## create a directory when it is not present
check_create_dir <- function(the_dir) {
  if (!dir.exists(the_dir)) {
    dir.create(the_dir, recursive = TRUE) }
}


##function_2
# create a function that reads a file, adds headers to it and save the file in a specific folder
summarize_data <- function(a_tab, the_dir) {
  # open the data, fix the date and add a new column
  the_data <- read.delim(a_tab, header = FALSE, na.strings = NA) %>%
    select(V1, V2,V3, V4,V5,V6) %>%
    rename(Percentage.reads = V1,
           Assigned.reads = V2,
           reads.taxon = V3,
           Taxon.rank = V4,
           Taxon.ID = V5,
           Taxon.name = V6)
  
  # write the csv to a new file
  write_csv(the_data, file.path(the_dir, basename(a_tab)),
            na = "na")
}

#######################################
#Processing data
#######################################

### What are the in and output directories
the_dir <- "data" # input data folder
the_dir_ex <- "data_edit/"   # output folder for the processed files
check_create_dir(the_dir_ex)


# get a list of all files that you want to process.
# here a pattern is used to select a subset of the files. It can be removed
# you use the list  of files
all_precip_files <- list.files(the_dir, pattern = "*-MOCK*",  # filtering on samples MOCK*
                               full.names = TRUE)

# use lapply to reformat the kraken results files, from the list created above.
lapply(all_precip_files,
       FUN = summarize_data,
       the_dir = the_dir_ex)

# reading the modified files into one dataframe
data_path <- "data_edit"   # path to the data
files <- list.files(data_path, pattern = "*-MOCK*") # get subset of file names, from all file names.

# creating a object called: data
data <- tibble(filename = files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename, ~ read.csv(file.path(data_path, .))) # a new data column
  ) %>%
  unnest(1:2)

## reformating the first column by chopping the last 14 characters
## needed since we do not want complete file names to work with later on.
data$filename <- fct_explicit_na(as.factor(gsub('.{14}$', '', data$filename)))

## reformating the taxon.name column by removing the spaces present in the original kraken output column
data$Taxon.name <- as.factor(gsub('  ', '', data$Taxon.name, fixed = TRUE))

## saving the final table to your working directory with
write_delim(data, file.path("mock_samples_long.tab"),
            na = "na", quote=FALSE )


################################################
### summarizing counts for the datasets     ####
################################################

#Summarizing read counts with dplyr
sample_counts <- data %>%
  group_by(filename) %>%
  summarise(total_reads=sum(reads.taxon)) #reads.taxon is the reads assigned to each taxonomic level

str(sample_counts)
plot(sample_counts$filename,sample_counts$total_reads) # shows the total fragments for each dataset.

# plotting total counts
ggplot(sample_counts, aes(y = total_reads, x=filename)) + 
  geom_col( colour="black", fill="indianred")+
  ggtitle("Distribution of sample sequencing depth") + 
  ylab("Total reads (n)") +
  xlab("Samples") +
  theme(axis.title.x = element_blank())

##Save as png in working directory.
ggsave("Sequencing_depth.png",
       width = 15,
       height = 11.25,
       units = "in",
       dpi = 300)
```

