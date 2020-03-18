##### 
# Variants: Verification of Clean Data
# 
# Author: Nicole Kauer
# Date: July 16, 2019
#
# Purpose: Read in the variant files from the folder. Verify the following:
#     1) Excel sheet/tab exists with the word 'missense' in the name.
#     2) Excel sheet/tab exists with the word 'indel' in the name.
#     3) Patient ID number in file name matches patient ID number found in the 
#        file under the column with the word 'sample' in the name.
#
#     Program also allows user to update a remove file, which lists the genes 
#     that do not belong in the dataset. 
# 
# Execute from terminal:
#     'Rscript Variant_Source_to_CSV.r <variant file directory in quotes>'
# Input(s): var.dir - directory for variant files
# Output(s): Currently prints error messages if something is wrong
#
# In progress: Interaction portion for printing of gene names and updating removal list.
#####

# Package(s) needed for reading in files
# For reading Excel files
if (!require("readxl", quietly=TRUE)) {
  install.packages("readxl", quiet=TRUE, repos='http://cran.us.r-project.org')
}
library("readxl", quietly=TRUE)
# For string matching
if(!require("stringr", quietly=TRUE)) {
  install.packages("stringr", quiet=TRUE, repos='http://cran.us.r-project.org')
}
library("stringr", quietly=TRUE)

# Get the directory folder for the variant data
# Note: Not paying accounting for input errors. Assuming input is exactly correct.
var.dir <- commandArgs(trailingOnly=TRUE)
var.dir <- paste0(var.dir, "/")

# The output directory is in the same input directory
output.dir <- paste0(var.dir, "Output/")

# Check if output directory already exist; if no, make it
if (!dir.exists(output.dir)) {
  dir.create(output.dir)
}

### FOR TESTING ONLY ###
#var.dir <- "C:\\Users\\nmeka\\Documents\\AML\\PamBeckerData\\Variants\\"
#output.dir <- "C:\\Users\\nmeka\\Documents\\AML\\PamBeckerData\\Variants\\Output\\"

# Check to see if this code has been run before and "Clean_Var.txt" was created. If yes, remove file.
# This file should only exist if this code has been run and no issues found.
if (file.exists(paste0(output.dir, "Clean_Var.txt"))) {
  file.remove(paste0(output.dir, "Clean_Var.txt"))
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
  sheets <- excel_sheets(paste0(var.dir,data.file))
  
  ### Should have message/notice that says using only missense and indel variants
  # Find missense sheet in Excel file; stop and request fix if can't find
  if (length(which(str_detect(sheets, coll("missense", ignore_case = TRUE)))) > 0) {
    sheet.miss <- which(str_detect(sheets, coll("missense", ignore_case = TRUE)))  
  } else {
    #print(paste0("Unable to find missense data in ", data.file, ". Please ensure that the sheet/tab has the word 'missense' in the name."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.miss, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }
  # Find the indel sheet in Excel file; stop and request fix if can't find
  if (length(which(str_detect(sheets, coll("indel", ignore_case = TRUE)))) > 0) {
    sheet.indel <- which(str_detect(sheets, coll("indel", ignore_case = TRUE)))   
  } else {
    #print(paste0("Unable to find indel data in ", data.file, ". Please ensure that the sheet/tab has the word 'indel' in the name."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.ind, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }
  
  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))

  # Add patient to list if not NA
  if (!is.na(patient)) {
    patients <- append(patients, patient)  
  }
  
  # Open file with missense info
  var.data <- read_excel(paste0(var.dir, data.file), sheet=sheets[sheet.miss])
  # Append genes to list
  genes <- append(genes, var.data$Gene)
  # Open file with indel info
  var.data <- read_excel(paste0(var.dir, data.file), sheet=sheets[sheet.indel])
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
      #print(paste0("The patient ID number in the actual name of the file is not the same as the patient ID number listed in the ", colnames(var.data)[col.index], " column for file ", data.file, ". The file has patient ID ", patient, " and the file has patient ID ", sample.patient, " Please reconcile the issue."))
      
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(var.data)[col.index], id1=patient, id2=sample.patient))
      # Set clean flag to false
      clean <- FALSE
    } # else they match and the patient number should be fine
  } else {
    # Was unable to find 'sample' column and cannot verify patient ID number.
    #print(paste0("Unable to identify a patient number inside file ", data.file, ". Please check that there is a Sample column with the patient ID."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.id.verify, filename=data.file, column=colnames(var.data)[col.index], id1=patient, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }
}

# There may be duplicate patient numbers (duplicate files, etc); get unique set
patients.all <- unique(patients)
# Sort
patients.all <- sort(patients.all)

write.table(patients.all, file = paste0(output.dir, "Patients_Variants.txt"), row.names = FALSE, col.names = FALSE, sep = ",")

# Get unique set of genes
genes.all <- unique(genes)

#### Initial gene removal. This is for genes that do not belong in the dataset at all.
# Check if file exists; if not, make it
if (file.exists(paste0(output.dir, "Genes_Remove.txt"))) {
  remove.genes <- read.csv(paste0(output.dir,"Genes_Remove.txt"), header = FALSE, stringsAsFactors = FALSE)
  # Check if file is empty 
  if (length(remove.genes) > 0) {
    # Make into a list instead of a table; technically not necessary, but seems less...odd
    remove.genes <- remove.genes[,1]
    # Remove each compound column in list 
    for (gene in remove.genes) {
      # Get index for gene, if it exists
      gene.index <- which(tolower(genes.all) == tolower(gene))
      # Check that the drug existed in dataset
      if (length(gene.index > 0)) {
        # Remove the gene from the full set
        genes.all <- genes.all[-gene.index]
      } # else gene wasn't found for some reason and nothing should be done
    }
  } # else file was empty and move on
} else {
  # Make blank file
  file.create(paste0("Genes_Remove.txt"))
}

# Sort
genes.all <- sort(genes.all)

write(genes.all, file = paste0(output.dir, "Genes.txt"))

#### Give list with check boxes; check to remove. Autocheck OBSCN

#### Second stage gene removal. This is for genes that do not belong in the dataset at all.
remove.genes <- read.csv(paste0(output.dir,"Genes_Remove.txt"), header = FALSE, stringsAsFactors = FALSE)
# Check if file is empty 
if (length(remove.genes) > 0) {
  # Make into a list instead of a table; technically not necessary, but seems less...odd
  remove.genes <- remove.genes[,1]
  # Remove each compound column in list 
  for (gene in remove.genes) {
    # Get index for gene, if it exists
    gene.index <- which(tolower(genes.all) == tolower(gene))
    # Check that the drug existed in dataset
    if (length(gene.index > 0)) {
      # Remove the gene from the full set
      genes.all <- genes.all[-gene.index]
    } # else gene wasn't found for some reason and nothing should be done
  }
} # else file was empty and move on

write(genes.all, file = paste0(output.dir, "Genes.txt"))

# Save a done file to use as a trigger iff there were no problems
if (clean) {
  file.create(paste0(output.dir, "Clean_DS.txt"))
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
  string <- "The following errors should be reconciled before continuing on with the analysis. Please run the cleaning code, again, to ensure there are no more errors."
  
  # First set of errors
  err.ind <- which(errors[,1] == error.miss)
  if (length(err.ind) > 0) {
    miss.errors.message <- "\nNo missense variant data was found in the following Excel files. Please ensure a sheet/tab with the data has the word missense in the name."
    
    miss.errors.list <- ""
    for (i in 1:length(err.ind)) {
      miss.errors.list <- paste(miss.errors.list, errors[err.ind[i],2], sep = "\n")  
    }
    #miss.errors.list <- paste(errors[err.ind,2], sep = "\n")  
    string <- paste(string, miss.errors.message, miss.errors.list, sep = "\n")
  }
  
  # Second set of errors
  err.ind <- which(errors[,1] == error.ind)
  if (length(err.ind) > 0) {
    ind.errors.message <- "\nNo indel variant data was found in the following Excel files. Please ensure a sheet/tab with the data has the word indel in the name."
    
    ind.errors.list <- ""
    for (i in 1:length(err.ind)) {
      ind.errors.list <- paste(ind.errors.list, errors[err.ind[i], 2], sep = "\n")  
    }
    #ind.errors.list <- paste(errors[err.ind, 2], sep = "\n")  
    string <- paste(string, ind.errors.message, ind.errors.list, sep = "\n")
  }
  
  # Third set of errors
  err.ind <- which(errors[,1] == error.id) 
  if (length(err.ind) > 0) {
    id.errors.message <- "\nThe patient number extracted from the filename was not the same as the patient number extracted from inside the file. Two columns were checked for, a column with the word File in it and a column with the word Sample in it. For the following files, the patient IDs from these locations did not match. "
    
    id.errors.list = ""
    for (i in 1:length(err.ind)) {
      temp = paste0(errors[err.ind[i], 2], " - ID from filename: ", errors[err.ind[i], 4], "; ID from column ", errors[err.ind[i], 3], ": ", errors[err.ind[i], 5])
      id.errors.list <- paste(id.errors.list, temp, sep = "\n")
    }
    #id.errors.list <- paste(errors[err.ind, 2], "- ID from filename: ", errors[err.ind, 4], "; ID from column ", errors[err.ind, 3], ": ", errors[err.ind, 5], sep = "")  
    string <- paste(string, id.errors.message, id.errors.list, sep = "\n")
  }
  
  # Fourth set of errors
  err.ind <- which(errors[,1] == error.id.verify)
  if (length(err.ind) > 0) {
    id.verify.errors.message <- "\nThere was no column with the word sample in the name found in the following files. Due to this, the patient ID could not be verified."
    
    id.verify.errors.list <- ""
    for (i in 1:err.ind) {
      id.verify.errors.list <- paste(id.verify.errors.list, errors[err.ind, 2], sep = "\n")
    }
    #id.verify.errors.list <- paste(errors[err.ind, 2], sep = "\n")
    string <- paste(string, id.verify.errors.message, id.verify.errors.list, sep = "\n")
  }

  cat(string, file = paste0(output.dir, "Variant File Errors.txt"), append = FALSE)
}
