###############################################################
# Project : Differential Gene Expression Analysis using DESeq2
# Dataset : GSE220113 (Human RNA-seq)
# Study   : Smoking vs Non-Smoking Peripheral Blood B Cells
###############################################################

#==============================================================
# 1. Clear Workspace
#==============================================================

# Remove all existing objects from R memory
rm(list = ls())

#==============================================================
# 2. Load Required Libraries
#==============================================================

# DESeq2  : Differential expression analysis
# readr   : Reading data files
# dplyr   : Data manipulation

library(DESeq2)
library(readr)
library(dplyr)

#==============================================================
# 3. Set Working Directory
#==============================================================

# Change this path according to your computer

setwd("C:/Users/patel/OneDrive/Desktop/human_deseq2")

#==============================================================
# 4. Read Raw Count Matrix
#==============================================================

# Read tab-separated RNA-seq count matrix
#
# header = TRUE
# -> First row contains column names
#
# sep = "\t"
# -> Columns are separated by TAB characters

counts <- read.delim(
  "GSE220113_raw_counts_GRCh38.p13_NCBI.tsv",
  header = TRUE,
  sep = "\t"
)

#==============================================================
# 5. Explore Count Matrix
#==============================================================

# Number of genes and samples

dim(counts)

# Display first six rows

head(counts)

# Display all column names

colnames(counts)

# Number of sample columns
# (-1 because first column is GeneID)

ncol(counts) - 1

#==============================================================
# 6. Read GEO Series Matrix
#==============================================================

# Read metadata file line by line

metadata_lines <- readLines("GSE220113_series_matrix.txt")

# Total number of lines

length(metadata_lines)

#==============================================================
# 7. Extract Smoking Status
#==============================================================

# Find the line containing
# "variable-condition"

condition_line <- metadata_lines[
  grep("variable-condition", metadata_lines)
]

# Split the line into separate fields

fields <- strsplit(condition_line, "\t")[[1]]

# Remove first element
# (!Sample_characteristics_ch1)

fields <- fields[-1]

# Remove quotation marks

fields <- gsub('"', "", fields)

# Remove text:
# "variable-condition: "

condition <- sub(
  "variable-condition: ",
  "",
  fields
)

# Check group sizes

table(condition)

#==============================================================
# 8. Extract GEO Sample IDs
#==============================================================

# Find GEO accession line

sample_line <- metadata_lines[
  grep("^!Sample_geo_accession", metadata_lines)
]

# Split into separate sample IDs

sample_ids <- strsplit(sample_line, "\t")[[1]]

# Remove first element

sample_ids <- sample_ids[-1]

# Remove quotation marks

sample_ids <- gsub('"', "", sample_ids)

# Verify number of GEO samples

length(sample_ids)

#==============================================================
# 9. Create Complete Metadata
#==============================================================

metadata_all <- data.frame(
  Sample = sample_ids,
  Smoking = factor(condition),
  stringsAsFactors = FALSE
)

# View first rows

head(metadata_all)

#==============================================================
# 10. Extract Sample Names from Count Matrix
#==============================================================

# First column is GeneID
# Remaining columns are sample names

sample_names <- colnames(counts)[-1]

length(sample_names)

head(sample_names)

#==============================================================
# 11. Match Metadata with Count Matrix
#==============================================================

# Keep only samples present in count matrix

metadata <- metadata_all[
  metadata_all$Sample %in% sample_names,
]

# Reorder metadata to match count matrix

metadata <- metadata[
  match(sample_names, metadata$Sample),
]

# Set row names

rownames(metadata) <- metadata$Sample

#==============================================================
# 12. Verify Metadata
#==============================================================

# Metadata dimensions

dim(metadata)

# First rows

head(metadata)

# Check sample order

all(sample_names == metadata$Sample)

# Number of smokers/non-smokers

table(metadata$Smoking)

#==============================================================
# 13. Prepare Count Matrix for DESeq2
#==============================================================

# Move GeneID column to row names

rownames(counts) <- counts$GeneID

# Remove GeneID column

counts <- counts[, -1]

# Check dimensions

dim(counts)

# Display first rows

head(counts)