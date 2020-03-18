###################################################################
#' @title Drug Sensitivity and Variant Gathering
#'
#' @author Nicole Kauer
#'
#' @description
#' Subsets the patient and drug dataset. Requires output directory
#' as an argument and assumes that the following files exist
#' in the directory: PatientsVsVariants.csv, PatientsVsDrugs.csv,
#' Patient_Subset.csv, Compound_Subset.csv.#'
#'
#' Execute from terminal:
#'     `Rscript DrugSensitivity_Variants_SubsetPatients.R <output dir>``
#'     
#' Outputs:
#'   PatientsVsVariants.csv: Overwrites original file with new
#'     subset of patients.
#'   PatientsVsDrugs.csv: Overwrites original file with new
#'     subset of patients and compounds.#'
###################################################################
library("fs", quietly=TRUE)
# Get the directory folder for the output data
# Note: Not paying accounting for input errors. Assuming input is exactly correct.
output.dir <- commandArgs(trailingOnly=TRUE)
output.dir <- fs::path(output.dir)

data.var <- read.csv(fs::path(output.dir, "PatientsVsVariants.csv"), header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
data.drug <- read.csv(fs::path(output.dir, "PatientsVsDrugs.csv"), header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

##### Subsetting data by patient list #####

# Read in file
patients.include <- read.csv(fs::path(output.dir, "Patient_Subset.csv"), check.names = FALSE)
# Get only patient IDs for which Include == 1
patient.subset <- which(patients.include$Include == 1)
# Keep only the patient subset desired in each dataset, assuming they chose at least one patient & not all patients
if (length(patient.subset) < dim(data.var)[1] && length(patient.subset) > 0) {
  data.drug <- data.drug[patient.subset,]
  data.var <- data.var[patient.subset,]
}

##### Subset data by drug compound list #####

# Read in file
compounds.include <- read.csv(fs::path(output.dir, "Compound_Subset.csv"), check.names = FALSE)
# Get only compounds for which Include == 1
compound.subset <- which(compounds.include$Include == 1)
# Keep only the drug subset desired in each dataset
# Make sure the subset is less than the current number of drugs & > 0
if (length(compound.subset) < dim(data.drug)[2] && length(compound.subset) > 0) {
  data.drug <- data.drug[,compound.subset]
}

##### Need to filter out compounds that have no useful data (all NA, all inactive drug values, etc)

# Remove columns that are all NA
#     since this subset of patients may not have data for certain drugs
drugs.NA.col <- apply(data.drug, 2, function (x) {sum(is.na(x))})
# If there is a drug that has no data, remove it; else move on
if (length(which(drugs.NA.col == dim(data.drug)[1])) > 0) {
  data.drug <- data.drug[, -which(drugs.NA.col == dim(data.drug)[1])]
}

# Filter out genes with no variants
# Count the number of 0s in each column
var.0.count <- apply(data.var, 2, function(x) {length(which(x == 0))})
# Find which columns have the same number of 0s as patients in the dataset
var.all.0 <- which(var.0.count == dim(data.var)[1])
# Remove genes that have no variants, if any
if (length(var.all.0) > 0) {
  data.var <- data.var[,-var.all.0]
}

##### Write to files #####

# Output filtered, subset, etc tables to files
write.csv(data.drug, file=fs::path(output.dir, "PatientsVsDrugs.csv"))
write.csv(data.var, file=fs::path(output.dir, "PatientsVsVariants.csv"))
