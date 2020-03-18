###################################################################
#' @title Validate clean drug sensitivity data
#'
#' @author Nicole Kauer
#'
#' @description
#' Validates that the drug sensitivity data is clean.
#' Outputs text file with errors if the data is not clean.
#' Also outputs a list of drug compounds and a synonym file.
#' Checks the following:
#'   1. Excel sheet exists with one of the following phrases
#'      in the name: `fitted parameters static`, `fitparameters``,
#'      `combined`, `fitted parameters`, `param statsort`,
#'      `parameters static`, `4Database`.
#'   2. Column exists in file with either `EC50` or `IC50`
#'      in the name, and the column has numeric values.
#'   3. Numeric patient ID is found in the filename and it
#'      matches the numeric patient ID found in a column with
#'      either the world `file` or `sample` in the name.
#'
#' This script requires a `directory`` argument that has the
#' location of drug senstivity files. Note that there should
#' be no other files, other than drug sensitivity, in the
#' directory.
#'
#' Outputs:
#'   Drugs_Remove.csv: Two column spreadsheet. Column `Drug`
#'     has all the unique drugs found across all files.
#'     Column `Include` has a 1 if the drug is to stay in
#'     the dataset; 0 otherwise.
#'   Drugs_Synonyms.csv: Spreadsheet where first column is
#'     assumed to be the preferred drug name and the
#'     following columns are the synonyms of the drug. If
#'     the file does not exist, this spreadsheet will be
#'     empty on initial creation.
#'   Drug_Sensitivity_File_Errors.txt: A text file with
#'     the errors found in the data. Only generated if
#'     data was unclean.
#'
#' To execute from terminal:
#'   `Rscript DrugSensitivity_CleanVerification.R <directory>`
###################################################################

# Line below is simply for RStudio. Turning off diagnostics because 
# keeps giving silly warnings; apparently, this is a known issue.
# !diagnostics off

library("readxl", quietly=TRUE)
library("stringr", quietly=TRUE)
library("fs", quietly=TRUE)

# Get the directory folder for the drug sensitivity data
# Note: Not paying accounting for input errors. Assuming input is exactly correct.
drug.dir <- commandArgs(trailingOnly=TRUE)
drug.dir <- as.character(fs::path(drug.dir))

# The output directory is in the same input directory
output.dir <- as.character(fs::path(drug.dir, "Output"))

# Check if output directory already exist; if no, make it
if (!fs::dir_exists(output.dir)) {
  fs::dir_create(output.dir)
}

# Check to see if this code has been run before and "Clean_DS.txt" was created. If yes, remove file.
# This file should only exist if this code has been run and no issues found.
if (fs::file_exists(as.character(fs::path(output.dir, "Drug_Sensitivity_File_Errors.txt")))) {
  # Invisible so it doesn't output a value
  invisible(
    fs::file_delete(
      as.character(fs::path(output.dir, "Drug_Sensitivity_File_Errors.txt"))
    )
  )
}

# Flag for checking if data had any errors
# Set to true at start and only changes to false if an error was found
clean <- TRUE

# Error codes
error.sheet <- 1
error.compound <- 2
error.EC50 <- 3
error.static <- 4
error.id <- 5
error.id.verify <- 6
error.auc <- 7

# Table for holding errors and information 
errors <- data.frame(error = numeric(), filename = character(),
                     column = character(), id1 = character(),
                     id2 = character(), stringsAsFactors = FALSE)

# Get a list of the files in the directory; divide up by file type
# Note: assuming there are only two types of files, xlsx and csv
filenames.excel <- list.files(drug.dir, pattern="*.xlsx")
filenames.csv <- list.files(drug.dir, pattern="*.csv")

# Vector for storing patient list
patients <- vector()
# Vector for storing compound list
compounds <- vector()

# Go through all Excel files, gathering patient IDs and compounds present
# Check whether there are issues and print message if there are
for (data.file in filenames.excel) {
  
  #print(data.file)
  
  # Get names of sheets in Excel workbook
  sheets <- excel_sheets(as.character(fs::path(drug.dir,data.file)))
  
  # Find the sheet that (hopefully) has static data based on previous sheets that (sometimes) has it
  # Assuming the first column with the name is the right one
  if (length(which(str_detect(sheets, coll("Fitted Parameters Static", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("Fitted Parameters Static", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("FitParameters", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("FitParameters", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("combined", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("combined", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("Fitted Parameters", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("Fitted Parameters", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("blast param statsort", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("blast param statsort", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("Parameters static", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("Parameters static", ignore_case = TRUE)))[1]
  }  else if (length(which(str_detect(sheets, coll("4Database", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("4Database", ignore_case = TRUE)))[1]
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.sheet, filename=data.file, column=NA, id1=NA, id2=NA))
    
    # Set clean flag to false
    clean <<- FALSE
  }
  
  # Open file
  drug.data <- suppressMessages(
    read_excel(
      as.character(fs::path(drug.dir, data.file)),
      sheet=sheets[sensitivity.sheet]
    )
  )
  
  # Get column that has drug compounds, if exists; else give error
  if (length(which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))) > 0) {
    col.compound <- which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))[1]
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.compound, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  # Get column that has either EC50 or IC50. Note, grabbing first instance 
  # since some have two sets, such as blood vs marrow samples. Also assuming
  # that the EC50 value comes before the standard deviation 
  # column that sometimes has EC50 in the name.
  if (length(which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))[1]
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.EC50, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  
  # Check if the number of NAs is the same as number of drugs
  if (sum(is.na(suppressWarnings(as.numeric(drug.data[[col.EC50]])))) == dim(drug.data)[1]) {
    # Add filename and error code to errors table
    errors <- rbind(
      errors,
      data.frame(
        error=error.static,
        filename=data.file,
        column=colnames(drug.data)[col.EC50],
        id1=NA,
        id2=NA
      )
    )
    # Set clean flag to false
    clean <<- FALSE
  } # else file has static data and should be fine
  
  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  
  #### Verify that user id from filename matches the user id from the file.
  #### Checking for column with word 'file'. If doesn't exist, check for column with word 'sample'.
  ####     If 'file' column exists, check that the value has 'AML' in it. If yes, extract patient ID. 
  ####     If no, look for 'sample' column and extract patient ID from that.
  # Check that a column exists in the sheet with the word 'File' in it.
  if (length(which(str_detect(colnames(drug.data), coll("File", ignore_case = TRUE)))) > 0) {
    # Get index of column that has the word 'file' in it.
    col.index <- which(str_detect(colnames(drug.data), coll("File", ignore_case = TRUE)))[1]
    # Get the first string in the column to extract patient number
    filename.string <- drug.data[1, col.index]
    # Check to see if AML is in the string, which generally means the ID is present in the string.
    #     If AML is not in the string, then check for the 'sample' column.
    if (str_detect(filename.string, "AML")) {
      # Grab the part of the string that has AML ###, including any space, underscore, etc between AML and ###
      filename.string <- str_extract(filename.string, "AML.[0-9]+")
    } else if (length(which(str_detect(colnames(drug.data), coll("Sample", ignore_case = TRUE)))) > 0) {
      col.index <- which(str_detect(colnames(drug.data), coll("Sample", ignore_case = TRUE)))[1]
      # Get the first string in the column to extract patient number
      filename.string <- drug.data[1, col.index]
      if (str_detect(filename.string, "AML")) {
        filename.string <- str_extract(filename.string, "AML.[0-9]+")
      } # else assume that the first numbers that appear in the file column string are the patient ID numbers
    }
    # Get the patient number from the string
    filename.patient <- as.numeric(str_extract(filename.string, "[0-9]+"))
    # Check if the patient ID number from the real filename is the same as the patient ID number listed in the 
    #     filename column in spreadsheet; check that it's not NA, either.
    if (is.na(filename.patient) || filename.patient != patient) {
      #print(paste0("The patient ID number in the actual name of the file is ", patient, ", which is not the same as the patient ID number listed (", filename.patient, ") in the ", colnames(drug.data)[col.index], " column for file ", data.file, ". Please reconcile the issue."))
      
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(drug.data)[col.index], id1=patient, id2=filename.patient))
      # Set clean flag to false
      clean <<- FALSE
    } # else they match and the patient number should be fine
  } else if (length(which(str_detect(colnames(drug.data), coll("Sample", ignore_case = TRUE)))) > 0) {
    col.index <- which(str_detect(colnames(drug.data), coll("Sample", ignore_case = TRUE)))[1]
    # Get the first string in the column to extract patient number
    filename.string <- drug.data[1, col.index]
    if (str_detect(filename.string, "AML")) {
      filename.string <- str_extract(filename.string, "AML.[0-9]+")
    }
    # Get the patient number from the string
    filename.patient <- as.numeric(str_extract(filename.string, "[0-9]+"))
    # Check if the patient ID number from the real filename is the same as the patient ID number listed in the 
    #     filename column in spreadsheet; check that it's not NA, either.
    if (is.na(filename.patient) || filename.patient != patient) {
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(drug.data)[col.index], id1=patient, id2=filename.patient))
      # Set clean flag to false
      clean <<- FALSE
    } # else they match and the patient number should be fine
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.id.verify, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  
  if (length(which(str_detect(colnames(drug.data), coll("AUC", ignore_case = TRUE)))) > 0) {
    col.AUC <- which(str_detect(colnames(drug.data), coll("AUC", ignore_case = TRUE)))[1]
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.auc, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }

  # Add patient to list
  patients <- append(patients, patient)
    
  # Append drug compounds to list
  for (i in 1:nrow(drug.data)) {
    compounds <- append(compounds, str_squish(str_trim(as.character(drug.data[i,col.compound]))))
  }
}

# Go through all csv files, gathering patient IDs and compounds present
for (data.file in filenames.csv) {
  
  # Open file
  drug.data <- read.csv(as.character(fs::path(drug.dir, data.file)), check.names = FALSE, stringsAsFactors = FALSE)
  
  # Get column that has drug compounds, if exists; else give error
  if (length(which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))) > 0) {
    col.compound <- which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))[1]
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.compound, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  # Get column that has either EC50 or IC50. Note, grabbing first instance since some have two sets, such
  #     as blood vs marrow samples. Also assuming that the EC50 value comes before the standard deviation 
  #     column that sometimes has EC50 in the name.
  if (length(which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))[1]
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.EC50, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  
  #### Check that the data is not NA for everything here. Give notice to user if it is with a reminder of
  ####     the need for static data. 
  # Check if the number of NAs is the same as number of drugs
  if (sum(is.na(suppressWarnings(as.numeric(drug.data[,col.EC50])))) == dim(drug.data)[1]) {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.static, filename=data.file, column=colnames(drug.data)[col.EC50], id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  } # else file has static data and should be fine (assuming it read in 'NA' instead of something else for irrelevant data)

  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  #### Verify that user id from filename matches the user id from the file.
  #### Checking for column with word 'file'. If doesn't exist, check for column with word 'sample'.
  ####     If 'file' column exists, check that the value has 'AML' in it. If yes, extract patient ID. 
  ####     If no, look for 'sample' column and extract patient ID from that.
  # Check that a column exists in the sheet with the word 'File' in it.
  if (length(which(str_detect(colnames(drug.data), coll("File", ignore_case = TRUE)))) > 0) {
    # Get index of column that has the word 'file' in it.
    col.index <- which(str_detect(colnames(drug.data), coll("File", ignore_case = TRUE)))[1]
    # Get the first string in the column to extract patient number
    filename.string <- drug.data[1, col.index]
    # Check to see if AML is in the string, which generally means the ID is present in the string.
    #     If AML is not in the string, then check for the 'sample' column.
    if (str_detect(filename.string, "AML")) {
      # Grab the part of the string that has AML ###, including any space, underscore, etc between AML and ###
      filename.string <- str_extract(filename.string, "AML.[0-9]+")
    } else if (length(which(str_detect(colnames(drug.data), coll("Sample", ignore_case = TRUE)))) > 0) {
      col.index <- which(str_detect(colnames(drug.data), coll("Sample", ignore_case = TRUE)))[1]
      # Get the first string in the column to extract patient number
      filename.string <- drug.data[1, col.index]
      if (str_detect(filename.string, "AML")) {
        filename.string <- str_extract(filename.string, "AML.[0-9]+")
      } # else assume that the first numbers that appear in the file column string are the patient ID numbers
    }
    # Get the patient number from the string
    filename.patient <- as.numeric(str_extract(filename.string, "[0-9]+"))
    # Check if the patient ID number from the real filename is the same as the patient ID number listed in the 
    #     filename column in spreadsheet; check that it's not NA, either.
    if (is.na(filename.patient) || filename.patient != patient) {
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(drug.data)[col.index], id1=patient, id2=filename.patient))
      # Set clean flag to false
      clean <<- FALSE
    } # else they match and the patient number should be fine
  } else if (length(which(str_detect(colnames(drug.data), coll("Sample", ignore_case = TRUE)))) > 0) {
    col.index <- which(str_detect(colnames(drug.data), coll("Sample", ignore_case = TRUE)))[1]
    # Get the first string in the column to extract patient number
    filename.string <- drug.data[1, col.index]
    if (str_detect(filename.string, "AML")) {
      filename.string <- str_extract(filename.string, "AML.[0-9]+")
    }
    # Get the patient number from the string
    filename.patient <- as.numeric(str_extract(filename.string, "[0-9]+"))
    # Check if the patient ID number from the real filename is the same as the patient ID number listed in the 
    #     filename column in spreadsheet
    if (is.na(filename.patient) || filename.patient != patient) {
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(drug.data)[col.index], id1=patient, id2=filename.patient))
      # Set clean flag to false
      clean <<- FALSE
    } # else they match and the patient number should be fine
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.id.verify, filename=data.file, column=NA, id1=patient, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  
  if (length(which(str_detect(colnames(drug.data), coll("AUC", ignore_case = TRUE)))) > 0) {
    col.AUC <- which(str_detect(colnames(drug.data), coll("AUC", ignore_case = TRUE)))[1]
  } else {
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.auc, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <<- FALSE
  }
  
  # Add patient to list
  patients <- append(patients, patient)
  
  # Append drug compounds to list
  for (i in 1:nrow(drug.data)) {
    compounds <- append(compounds, str_squish(str_trim(as.character(drug.data[i,col.compound]))))
  }
}

#### Clean up of patient and compound vectors
# There may be duplicate patient numbers (duplicate files, etc); get unique set
patients.all <- unique(patients)
patients.all <- sort(patients.all)
write.table(patients.all, file = as.character(fs::path(output.dir, "Patients_DrugSensitivity.txt")), row.names = FALSE, col.names = FALSE, sep = ",")

# Get unique set of compounds 
# Ignore NA compounds and case
compounds.all <- unique(compounds, ignore_case=TRUE)
# Remove all the excess numbered "compounds" generated by functions opening files.
#    Suppress warnings about NAs since that is the purpose of this: see if the compound
#    name can be converted to numeric without turning into an NA. If so, ignore it; if not,
#    keep it since it has letters (or at least something other than numbers).
compounds.all <- compounds.all[is.na(suppressWarnings(as.numeric(compounds.all)))]

# Sort
compounds.all <- sort(compounds.all)
#write.csv(compounds.all, paste0(output.dir, "Compound_all.csv"), row.names = FALSE)

# List of indices to remove
comp.ind <- vector()

#### Do initial synonym file stuff here
# Check if the file exists first; if not, make it.
if (fs::file_exists(as.character(fs::path(output.dir, "Compound_Synonyms.csv")))) {
  # Get information on file
  syn.info <- file.info(as.character(fs::path(output.dir, "Compound_Synonyms.csv")))
  # Make sure synonym file is not empty; R reads this in different ways sometimes, hence the multiple conditions
  if (syn.info$size > 0 && !is.na(syn.info$size)) {
    # Open synonym file
    synonyms <- read.csv(as.character(fs::path(output.dir, "Compound_Synonyms.csv")), header = FALSE, na.strings = c("", " ", NA), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    
    # Add all compounds in first column of synonym list
    #    This allows for names that users want to be shorter but never appear in compound list in files
    #   Only add the rows that have something in them
    compounds.all <- append(compounds.all, synonyms[,1])
    # Get the unique set in case synonym file added a double
    compounds.all <- unique(compounds.all, ignore_case=TRUE)
    # Sort again
    compounds.all <- sort(compounds.all)
    
    # Transform compounds to be lowercase and have no punctuation
    #compounds.all.t <- gsub("[[:punct:] ]+", "", compounds.all)
    compounds.all.t <- tolower(compounds.all)
    
    # Remove all synonym compounds; assume first column in synonyms file is primary name
    # Go through each row, ignoring first column in the row and any NA values
    for (i in 1:nrow(synonyms)) {
      # Go through each synonym and check the location of the name in the full compounds set
      for (synonym in synonyms[i,-1]) {
        # Check if the synonym is NA
        if (!is.na(synonym)) {
          # Transform synonym the same way as the compounds
          #synonym.t <- gsub("[[:punct:] ]+", "", synonym)
          synonym.t <- tolower(synonym)
          # Get index for where name appears in full compound set
          syn.index <- grep(synonym.t, compounds.all.t, fixed = TRUE)
          # Check that an appropriate index was actually found (i.e. name exists in list)
          if (length(syn.index) > 0) {
            # Add to list for removal
            comp.ind <- append(comp.ind, syn.index)
          }  
        }
      }
    } 
  } # else file is empty and move on
} else {
  # Create blank file
  file.create(as.character(fs::path(output.dir, "Compound_Synonyms.csv")))
}

# Remove compounds
if (length(comp.ind) > 0) {
  compounds.all <- compounds.all[-comp.ind] 
}

### Do initial compound removal here. This is for compounds that should not be included 
###    in the dataset at all.
# Check that the file exists first; if not, make it
if (fs::file_exists(as.character(fs::path(output.dir, "Compounds_Remove.csv")))) {
  remove.compounds.file <- read.csv(as.character(fs::path(output.dir,"Compounds_Remove.csv")), header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  # Assumption: file is not empty
  # Get the genes that should be removed (i.e. Include = 0)
  compounds.remove <- remove.compounds.file$Compound[remove.compounds.file$Include == 0]
  all.compounds.included <- data.frame(Compound = compounds.all, Include = 1)
  if (length(compounds.remove) > 0) {
    # Set genes to remove that were previously set
    for (drug in compounds.remove) {
      drug.ind <- which(all.compounds.included$Compound == drug)
      if (length(drug.ind) > 0) {
        all.compounds.included$Include[drug.ind] <- 0 
      }
    }
    #includes <- apply(all.genes.included, function(x){ifelse(x[1] %in% genes.remove$Gene, 0, 1)})
    #all.genes.included$Include <- includes  
  }
  #all.compounds.included$Include <- lapply(all.compounds.included$Compound, function(x){ifelse(x %in% compounds.remove$Compound, 0, 1)})
  write.csv(all.compounds.included, file = as.character(fs::path(output.dir, "Compounds_Remove.csv")), row.names = FALSE)
} else {
  # Make dataframe with all genes included
  all.compounds.included <- data.frame(Compound = compounds.all, Include = 1)
  # Write file to output directory
  write.csv(all.compounds.included, file = as.character(fs::path(output.dir, "Compounds_Remove.csv")), row.names = FALSE)
}

# Save a done file to use as a trigger iff there were no problems
if (clean) {
  # Clean directory resides within output directory
  clean.dir <- as.character(fs::path(output.dir, "_Clean"))
  # Check if clean output directory already exist; if no, make it
  if (!fs::dir_exists(as.character(clean.dir))) {
    fs::dir_create(as.character(clean.dir))
  }
  # Copy files to clean folder; invisible so doesn't output values
  if (length(filenames.excel) > 0 ) {
    invisible(file.copy(as.character(fs::path(drug.dir, filenames.excel)), clean.dir)) 
  }
  if (length(filenames.csv) > 0 ) {
    invisible(file.copy(as.character(fs::path(drug.dir, filenames.csv)), clean.dir))  
  }
} else { # At least one error must exist
  # Start of string message
  string <- "The following errors should be reconciled before continuing on with the analysis. Please run the cleaning code, again, to ensure there are no more errors."
  
  # First set of errors
  err.ind <- which(errors[,1] == error.sheet)
  if (length(err.ind) > 0) {
    sheet.errors.message <- "\nNo drug sensitivity was found in the following Excel files. Please ensure a sheet/tab with the data has one of the following phrases in the name: Fitted Parameters Static, FitParameters, combined, Fitted Parameters, blast param statsort, or 4Database."
    
    sheet.errors.list <- ""
    for (i in 1:length(err.ind)) {
      sheet.errors.list <- paste(sheet.errors.list, errors[err.ind[i],2], sep = "\n")
    }
    string <- paste(string, sheet.errors.message, sheet.errors.list, sep = "\n")
  }
  
  # Second set of errors
  err.ind <- which(errors[,1] == error.compound)
  if (length(err.ind) > 0) {
    compound.errors.message <- "\nNo column called Compound was found in the drug sensitivity data for the following files."
    
    compound.errors.list <- ""
    for (i in 1:length(err.ind)) {
      compound.errors.list <- paste(compound.errors.list, errors[err.ind[i],2], sep = "\n")
    } 
    string <- paste(string, compound.errors.message, compound.errors.list, sep = "\n")
  }
  
  # Third set of errors
  err.ind <- which(errors[,1] == error.EC50)
  if (length(err.ind) > 0) {
    EC50.errors.message <- "\nNo column called either EC50 or IC50 was found in the drug sensitivity data for the following files."
    
    EC50.errors.list <- ""
    for (i in 1:length(err.ind)) {
      EC50.errors.list <- paste(EC50.errors.list, errors[err.ind[i],2], sep = "\n")
    }
    string <- paste(string, EC50.errors.message, EC50.errors.list, sep = "\n")
  }
  
  # Fourth set of errors
  err.ind <- which(errors[,1] == error.static)
  if (length(err.ind) > 0) {
    static.errors.message <- "\nThe EC50 or IC50 data does not appear to be numeric in the following files. The column name checked is in parentheses."
    
    static.errors.list <- ""
    for (i in 1:length(err.ind)) {
      temp <- paste0(errors[err.ind[i], 2], " (", errors[err.ind[i], 3], ")")
      static.errors.list <- paste(static.errors.list, temp, sep = "\n")
    }
    string <- paste(string, static.errors.message, static.errors.list,sep = "\n")
  }
  
  # Fifth set of errors
  err.ind <- which(errors[,1] == error.id) 
  if (length(err.ind) > 0) {
    id.errors.message <- "\nThe patient number extracted from the filename was not the same as the patient number extracted from inside the file. Two columns were checked for, a column with the word File in it and a column with the word Sample in it. For the following files, the patient IDs from these locations did not match. "
    
    id.errors.list <- ""
    for (i in 1:length(err.ind)) {
      temp <- paste0(errors[err.ind[i], 2], " - ID from filename: ", errors[err.ind[i], 4], "; ID from column ", errors[err.ind[i], 3], ": ", errors[err.ind[i], 5])
      id.errors.list <- paste(id.errors.list, temp, sep = "\n")
    }
    string <- paste(string, id.errors.message, id.errors.list, sep = "\n")
  }
  
  # Sixth set of errors
  err.ind <- which(errors[,1] == error.id.verify)
  if (length(err.ind) > 0) {
    id.verify.errors.message <- "\nThere was no column with the word File or Sample in the name found in the following files. Due to this, the patient ID could not be verified."
    
    id.verify.errors.list <- ""
    for (i in 1:length(err.ind)) {
      id.verify.errors.list <- paste(id.verify.errors.list, errors[err.ind[i],2], sep = "\n")
    }
    string <- paste(string, id.verify.errors.message, id.verify.errors.list, sep = "\n")
  }
  
  # Seventh set of errors
  err.ind <- which(errors[,1] == error.auc)
  if (length(err.ind) > 0) {
    auc.errors.message <- "\nThere was no column called 'AUC' in the drug sensivity data for the following file(s)."

    auc.errors.list <- ""
    for (i in 1:length(err.ind)) {
      auc.errors.list <- paste(auc.errors.list, errors[err.ind[i],2], sep = "\n")
    }
    string <- paste(string, auc.errors.message, auc.errors.list, sep = "\n")
  }
  
  cat(string, file = as.character(fs::path(output.dir, "Drug_Sensitivity_File_Errors.txt")), append = FALSE)
}