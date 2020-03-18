###################################################################
#' @title Validate clean variant data
#'
#' @author Nicole Kauer
#'
#' @description
#' Validates that the variant data is clean. Outputs text file
#' with errors if the data is not clean. Also outputs list of
#' genes.
#' Checks the following:
#'   1. Two excel sheets exist in the file, one with the word
#'      `missense` and the other with the word `indel` in the
#'      name.
#'   2. Numeric patient ID found in the filename, and the ID
#'      matches the numeric patient ID found in a column
#'      with the world `sample` in the name.
#' 
#' This script requires a `directory`` argument that has the
#' location of variant files. Note that there should
#' be no other files, other than variants, in the
#' directory.
#'
#' Outputs:
#'   Genes_Remove.csv: Two column spreadsheet. Column `Gene`
#'     has all the unique genes found across all files.
#'     Column `Include` has a 1 if the gene is to stay in
#'     the dataset; 0 otherwise.
#'   Variant_File_Errors.txt: A text file with the errors
#'   found in the data. Only generated if data was unclean.
#'
#' To execute from terminal:
#'   `Rscript Variants_CleanVerification.R <directory>`
#'
###################################################################

library("readxl", quietly = TRUE)
library("stringr", quietly = TRUE)
library("fs", quietly = TRUE)

# Get the directory folder for the variant data
# Note: Not paying accounting for input errors. Assuming input is exactly correct.
var.dir <- commandArgs(trailingOnly=TRUE)
var.dir <- fs::path(var.dir)

# The output directory is in the same input directory
output.dir <- fs::path(var.dir, "Output")

# Check if output directory already exist; if no, make it
if (!fs::dir_exists(output.dir)) {
  fs::dir_create(output.dir)
}

# Check to see if this code has been run before and "Clean_Var.txt" was created. If yes, remove file.
# This file should only exist if this code has been run and no issues found.
if (fs::file_exists(fs::path(output.dir, "Variant_File_Errors.txt"))) {
  # Invisible so it doesn't output a value
  invisible(file.remove(fs::path(output.dir, "Variant_File_Errors.txt")))
}

# Flag for checking if data had any errors
# Set to true at start and only changes to false if an error was found
clean <- TRUE

# Error codes
error.miss <- 1
error.ind <- 2
error.id <- 3
error.id.verify <- 4

# Table for holding errors and information 
errors <- data.frame(error = numeric(), filename = character(), column = character(), 
                     id1 = character(), id2 = character(), stringsAsFactors = FALSE)

# Get a list of the files in the directory
# Note: assuming only xlsx files
filenames <- list.files(var.dir, pattern = "*.xlsx")

# Vector for storing patient list
patients <- vector()
# Vector for storing gene list
genes <- vector()

# Go through all files, gathering patient IDs and genes present
for (data.file in filenames) {

  # Get names of sheets in Excel workbook
  sheets <- excel_sheets(fs::path(var.dir,data.file))
  
  ### Should have message/notice that says using only missense and indel variants
  # Find missense sheet in Excel file; stop and request fix if can't find
  if (length(which(str_detect(sheets, coll("missense", ignore_case = TRUE)))) > 0) {
    sheet.miss <- which(str_detect(sheets, coll("missense", ignore_case = TRUE)))  
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.miss, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  # Find the indel sheet in Excel file; stop and request fix if can't find
  if (length(which(str_detect(sheets, coll("indel", ignore_case = TRUE)))) > 0) {
    sheet.indel <- which(str_detect(sheets, coll("indel", ignore_case = TRUE)))   
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.ind, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  
  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))

  # Add patient to list if not NA
  if (!is.na(patient)) {
    patients <- append(patients, patient)  
  }
  
  # Open file with missense info
  var.data <- read_excel(fs::path(var.dir, data.file), sheet=sheets[sheet.miss])
  # Append genes to list
  genes <- append(genes, var.data$Gene)
  # Open file with indel info
  var.data <- read_excel(fs::path(var.dir, data.file), sheet=sheets[sheet.indel])
  # Append genes to list
  genes <- append(genes, var.data$Gene)
  
  #### Verify that user id from filename matches the user id from the file 
  # Check that a Sample column exists; generally called Sample_ID, but reducing to Sample just in case
  if (length(which(str_detect(colnames(var.data), coll("Sample", ignore_case = TRUE)))) > 0) {
    # Get index of column that has the word 'Sample' in it.
    col.index <- which(str_detect(colnames(var.data), coll("Sample", ignore_case = TRUE)))
    # Get the first string in the column to extract patient number
    sample.string <- var.data[1, col.index]
    # Check if sample ID is already just a number or if it's a string
    if (!is.na(suppressWarnings(as.numeric(sample.string)))) {
      sample.patient <- as.numeric(sample.string)
    } else {
      # Check to see if AML is in the string
      if (str_detect(sample.string, "AML")) {
        # Grab the part of the string that has AML ###, including any space, underscore, etc between AML and ###
        sample.string <- str_extract(sample.string, "AML.[0-9]+")
      } # else, just get the first number in the string
      # Get the patient number from the string
      sample.patient <- as.numeric(str_extract(sample.string, "[0-9]+"))
    }
    # Check if the patient ID number from the real filename is the same as the patient ID number listed in the 
    #     filename column in spreadsheet; also check that isn't NA
    if (is.na(sample.patient) || sample.patient != patient) {
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(var.data)[col.index], id1=patient, id2=sample.patient))
      # Set clean flag to false
      clean <<- FALSE
    } # else they match and the patient number should be fine
  } else {
    # Was unable to find 'sample' column and cannot verify patient ID number.
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.id.verify, filename=data.file, column=colnames(var.data)[col.index], id1=patient, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
}

# There may be duplicate patient numbers (duplicate files, etc); get unique set
patients.all <- unique(patients)
# Sort
patients.all <- sort(patients.all)

write.table(patients.all, file = fs::path(output.dir, "Patients_Variants.txt"), row.names = FALSE, col.names = FALSE, sep = ",")

# Get unique set of genes
genes.all <- unique(genes)

# Sort
genes.all <- sort(genes.all)

#### Initial gene removal. This is for genes that do not belong in the dataset at all.
#   On first run, this only creates the file. On the second run, it removes any that the user
#   already selected to remove.
# Check if file exists; if not, make it
if (fs::file_exists(fs::path(output.dir, "Genes_Remove.csv"))) {
  remove.genes.file <- read.csv(fs::path(output.dir,"Genes_Remove.csv"), header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  # Assumption: file is not empty
  # Get the genes that should be removed (i.e. Include = 0).
  genes.remove <- remove.genes.file[remove.genes.file$Include == 0, ]
  # Make dataframe with all genes included
  all.genes.included <- data.frame(Gene = genes.all, Include = 1)
  if (dim(genes.remove)[1] > 0) {
    # Set genes to remove that were previously set
    for (gene in genes.remove$Gene) {
      gene.ind <- which(all.genes.included$Gene == gene)
      all.genes.included$Include[gene.ind] <- 0
    }
  }
  
  write.csv(all.genes.included, file = fs::path(output.dir, "Genes_Remove.csv"), row.names = FALSE)
} else {
  # Make dataframe with all genes included
  all.genes.included <- data.frame(Gene = genes.all, Include = 1)
  # Write file to output directory
  write.csv(all.genes.included, file = fs::path(output.dir, "Genes_Remove.csv"), row.names = FALSE)
}

# Save clean variant files to new directory and save done file as trigger iff there were no problems
if (clean) {
  # Clean directory resides within output directory
  clean.dir <- fs::path(output.dir, "_Clean")
  # Check if clean output directory already exist; if no, make it
  if (!fs::dir_exists(clean.dir)) {
    dir.create(clean.dir)
  }
  # Copy files to clean folder; invisible so doesn't output values
  invisible(file.copy(fs::path(var.dir, filenames), clean.dir))
} else { # At least one error must exist
  # Print out errors to file
  # Note that the dataframe is being accessed by number instead of column name. This is due to the fact that it could
  # collapse into an atomic vector if only one error exists. Atomic vectors cannot be accessed by column name. This could
  # potentially be remedied by switching to a tibble.
  # Column - description
  #   1 - error code
  #   2 - filename
  #   3 - column name
  #   4 - patient ID from filename
  #   5 - pateint ID from inside file
  
  # Start of string message
  string <- "The following errors should be reconciled before continuing on with the analysis. Please run the cleaning code after resolving these errors."
  
  # Error message for not finding a missense tab
  err.ind <- which(errors[,1] == error.miss)
  if (length(err.ind) > 0) {
    miss.errors.message <- "\nNo missense variant data was found in the following Excel files. Please ensure a sheet/tab with the data has the word 'missense' in the name."
    
    miss.errors.list <- ""
    for (i in 1:length(err.ind)) {
      miss.errors.list <- paste(miss.errors.list, errors[err.ind[i],2], sep = "\n")  
    }
    string <- paste(string, miss.errors.message, miss.errors.list, sep = "\n")
  }
  
  # Error message for not finding an indel tab
  err.ind <- which(errors[,1] == error.ind)
  if (length(err.ind) > 0) {
    ind.errors.message <- "\nNo indel variant data was found in the following Excel files. Please ensure a sheet/tab with the data has the word 'indel' in the name."
    
    ind.errors.list <- ""
    for (i in 1:length(err.ind)) {
      ind.errors.list <- paste(ind.errors.list, errors[err.ind[i], 2], sep = "\n")  
    }
    string <- paste(string, ind.errors.message, ind.errors.list, sep = "\n")
  }
  
  # Error message for patient IDs not matching
  err.ind <- which(errors[,1] == error.id) 
  if (length(err.ind) > 0) {
    id.errors.message <- "\nThe patient number extracted from the filename was not the same as the patient number extracted from inside the file. Two columns were checked for, a column with the word File in it and a column with the word Sample in it. For the following files, the patient IDs from these locations did not match. "
    
    id.errors.list = ""
    for (i in 1:length(err.ind)) {
      temp = paste0(errors[err.ind[i], 2], " - ID from filename: ", errors[err.ind[i], 4], "; ID from column ", errors[err.ind[i], 3], ": ", errors[err.ind[i], 5])
      id.errors.list <- paste(id.errors.list, temp, sep = "\n")
    }
    string <- paste(string, id.errors.message, id.errors.list, sep = "\n")
  }
  
  # Error message for not being able to verify patient ID
  err.ind <- which(errors[,1] == error.id.verify)
  if (length(err.ind) > 0) {
    id.verify.errors.message <- "\nThere was no column with the word sample in the name found in the following files. Due to this, the patient ID could not be verified."
    
    id.verify.errors.list <- ""
    for (i in 1:err.ind) {
      id.verify.errors.list <- paste(id.verify.errors.list, errors[err.ind, 2], sep = "\n")
    }
    string <- paste(string, id.verify.errors.message, id.verify.errors.list, sep = "\n")
  }

  cat(string, file = paste0(output.dir, "Variant_File_Errors.txt"), append = FALSE)
}
