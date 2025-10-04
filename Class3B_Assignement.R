

# Bioconductor provides R packages for analyzing omics data (genomics, transcriptomics, proteomics etc).

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))

# Install CRAN packages for data manipulation
install.packages("dplyr")

# Load Required Libraries
library(GEOquery)             # Download GEO datasets (series matrix, raw CEL files)
library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(dplyr)                # Data manipulation

# -------------------------------------
#### Download Series Matrix Files ####
# -------------------------------------


# GSEMatrix = TRUE means return the series matrix objects which contain process expression data along with phentorypic and feature data
gse_data <- getGEO("GSE110223", GSEMatrix = TRUE) 
gse_data
# Extract expression data matrix (genes/probes × samples)
# Rows corresponds to probes and columns corresponds to samples
expression_data <- exprs(gse_data$GSE110223_series_matrix.txt.gz)


# Extract feature (probe annotation) data
# Rows corresponds to probes and columns corresponds to samples
feature_data <- fData(gse_data$GSE110223_series_matrix.txt.gz)

# Extract phenotype (sample metadata) data
# Rows corresponds to samples and columns corresponds to probes
phenotype_data <-  pData(gse_data$GSE110223_series_matrix.txt.gz)

# Check missing values in sample annotation
sum(is.na(phenotype_data$source_name_ch1))
# --------------------------------------
#### Download Raw Data (CEL files) ####
# --------------------------------------

# Fetch GEO supplementry files
getGEOSuppFiles("GSE110223", baseDir = "Raw_Data", makeDirectory = TRUE)
# Important Note: 
# For Affymetrix microarray data, the preprocessing pipeline is the same 
# whether raw CEL files are downloaded from NCBI GEO or ArrayExpress. 

# (See tutorial for detailed explanation of this step: https://youtu.be/DZMxkHxwWag?feature=shared) 

# Untar CEL files if compressed as .tar
untar("C:/Users/ASUS/Desktop/AI & Biotech training/AI_Omics_internship_2025/Microarray_Analysis/Assignement/Raw_Data/GSE110223_RAW.tar", exdir = "C:/Users/ASUS/Desktop/AI & Biotech training/AI_Omics_internship_2025/Microarray_Analysis/Assignement/Raw_Data/CEL_Files")


# Read CEL files into R as an AffyBatch object
raw_data <- ReadAffy(celfile.path = "C:/Users/ASUS/Desktop/AI & Biotech training/AI_Omics_internship_2025/Microarray_Analysis/Assignement/Raw_Data/CEL_Files")

raw_data   # Displays basic information about the dataset
.

# ---------------------------------------------------
#### Quality Control (QC) Before Pre-processing ####
# ---------------------------------------------------

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "C:/Users/ASUS/Desktop/AI & Biotech training/AI_Omics_internship_2025/Microarray_Analysis/Assignement/Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)
#To know the number of probes 
probes_data <- as.data.frame(exprs(raw_data))

dim(probes_data)

# -------------------------------------------------------
#### RMA (Robust Multi-array Average) Normalization ####
# -------------------------------------------------------
# This method reduces experimental variation across multiple arrays, 
# producing more symmetrical and reliable normalized expression data 
# compared to other approaches
library(affy)
normalized_data <- affy::rma(raw_data)

# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "C:/Users/ASUS/Desktop/AI & Biotech training/AI_Omics_internship_2025/Microarray_Analysis/Assignement/Results/QC_Normalized_Data",
                    force = TRUE)
# Extract normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)   # Dimensions: number of probes × number of samples

# ---------------------------------------------------------------------------
#### Filter Low-Variance Transcripts (“soft” intensity based filtering) ####
# ---------------------------------------------------------------------------

# Filtering removes probes with low or uninformative expression signals.

# Calculate median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))
row_median

# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

# Set a threshold to remove low variance probes (dataset-specific, adjust accordingly)
threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2) 

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

# Rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

# Overwrite processed data with filtered dataset
processed_data <- filtered_data 
dim(processed_data)

# -----------------------------------
#### Phenotype Data Preparation ####
# -----------------------------------

# Phenotype data contains sample-level metadata such as condition, 
# tissue type, or disease status.
# Required to define experimental groups for statistical analysis.

class(phenotype_data$source_name_ch1) 

# Define experimental groups (normal vs cancer)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("normal adjacent", "primary colorectal adenocarcinoma"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)
