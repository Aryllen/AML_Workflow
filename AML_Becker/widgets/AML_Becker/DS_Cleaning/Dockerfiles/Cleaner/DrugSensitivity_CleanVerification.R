##### 
# Drug Sensitivity: Verification of Clean Data
# 
# Author: Nicole Kauer
# Date: July 14, 2019
#
# Purpose: Read in the drug sensitivity files from the folder. Verify the following:
#     1) Excel sheet/tab exists with one of the following phrases in the name: 
#           - fitted parameters static
#           - fitparameters
#           - combined
#           - fitted parameters
#           - param statsort
#           - parameters static
#           - 4Database
#     2) A column exists in the file with either 'EC50' or 'IC50' in the name.
#     3) The EC50/IC50 column has numeric values.
#     4) Patient ID in the file name matches the patient ID found in a column with 
#        either the word 'file' or 'sample' in the name.
#
#     Program also allows user to update a synonym file and remove file. The synonym file
#     lists the extra names for drug compounds. The remove file lists the drugs that should
#     not be included in the dataset.
# 
# Execute from terminal:
#     'Rscript Variant_Source_to_CSV.R <drug sensitivity file directory in quotes>'
# Input(s): drug.dir - directory for drug sensivity files
# Output(s): Currently prints error messages if something is not correct
#
# In progress: Interactive portion for printing of drug compounds and updating synonym/remove files.
#####

# Line below is simply for RStudio. Turning off diagnostics because keeps giving
# annoying, useless warnings; apparently, this is a known issue.
# !diagnostics off

# Package(s) needed for reading in files
# For reading Excel files
# if (!require("readxl", quietly=TRUE)) {
#   install.packages("readxl", quiet=TRUE, repos='http://cran.us.r-project.org')
# }
library("readxl", quietly=TRUE)
# For string matching
# if(!require("stringr", quietly=TRUE)) {
#   install.packages("stringr", quiet=TRUE, repos='http://cran.us.r-project.org')
# }
library("stringr", quietly=TRUE)

# Get the directory folder for the drug sensitivity data
# Note: Not paying accounting for input errors. Assuming input is exactly correct.
drug.dir <- commandArgs(trailingOnly=TRUE)
drug.dir <- paste0(drug.dir, "/")

# The output directory is in the same input directory
output.dir <- paste0(drug.dir, "Output/")

# Check if output directory already exist; if no, make it
if (!dir.exists(output.dir)) {
  dir.create(output.dir)
}

### FOR TESTING ONLY ###
#drug.dir <- "C:\\Users\\nmeka\\Documents\\AML\\PamBeckerData\\DrugSensitivity\\"
#output.dir <- "C:\\Users\\nmeka\\Documents\\AML\\PamBeckerData\\DrugSensitivity\\Output\\"

# Check to see if this code has been run before and "Clean_DS.txt" was created. If yes, remove file.
# This file should only exist if this code has been run and no issues found.
if (file.exists(paste0(output.dir, "Clean_DS.txt"))) {
  file.remove(paste0(output.dir, "Clean_DS.txt"))
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

# Table for holding errors and information 
errors <- data.frame(error = numeric(), filename = character(), column = character(), 
                     id1 = character(), id2 = character(), stringsAsFactors = FALSE)

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
  sheets <- excel_sheets(paste0(drug.dir,data.file))
  
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
    # #### Need to have the program tell the user that they did not find a sheet in the file for
    # ####     the data. Have a reminder of what is expected/searched for.
    #print(paste0("The drug sensitivity data was not found in ", data.file, ". Please verify that the static data is located in a tab with one of the following phrases in the name: \"Fitted Parameters Static\", \"FitParameters\", \"combined\", \"Fitted Parameters\", \"blast param statsort\", or \"4Database\"."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.sheet, filename=data.file, column=NA, id1=NA, id2=NA))
    
    # Set clean flag to false
    clean <- FALSE
  }
  
  # Open file
  # Suppressing unimportant messages such as the notice that some columns were given names; in general,
  #    these were just "empty" columns that it decided was still a part of the file, most likely due to invisible characters.
  drug.data <- suppressMessages(read_excel(paste0(drug.dir, data.file), sheet=sheets[sensitivity.sheet]))
  
  # Get column that has drug compounds, if exists; else give error
  if (length(which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))) > 0) {
    col.compound <- which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))[1]
  } else {
    #print(paste0("No Compound column was found in ", data.file, ". Please reconcile."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.compound, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }
  # Get column that has either EC50 or IC50. Note, grabbing first instance since some have two sets, such
  #     as blood vs marrow samples. Also assuming that the EC50 value comes before the standard deviation 
  #     column that sometimes has EC50 in the name.
  if (length(which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))[1]
  } else {
    #print(paste0("No EC50 nor IC50 column was found in ", data.file, ". Please reconcile."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.EC50, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }
  
  #### Check that the data is not NA for everything here. Give notice to user if it is with a reminder of
  ####     the need for static data.
  # Check if the number of NAs is the same as number of drugs
  if (sum(is.na(suppressWarnings(as.numeric(drug.data[[col.EC50]])))) == dim(drug.data)[1]) {
    # Give message and request that the user fix the data
    #print(paste0("The ", colnames(drug.data)[col.EC50] , " column in the file ", data.file, " does not appear to have the necessary values. Please verify that the column has static, numeric values."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.static, filename=data.file, column=colnames(drug.data)[col.EC50], id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
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
      #print(paste0("The patient ID number in the actual name of the file is ", patient, ", which is not the same as the patient ID number listed (", filename.patient, ") in the ", colnames(drug.data)[col.index], " column for file ", data.file, ". Please reconcile the issue."))
      
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(drug.data)[col.index], id1=patient, id2=filename.patient))
      # Set clean flag to false
      clean <- FALSE
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
      #print(paste0("The patient ID number in the actual name of the file is ", patient, ", which is not the same as the patient ID number listed (", filename.patient, ") in the ", colnames(drug.data)[col.index], " column for file ", data.file, ". Please reconcile the issue."))
      
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(drug.data)[col.index], id1=patient, id2=filename.patient))
      # Set clean flag to false
      clean <- FALSE
    } # else they match and the patient number should be fine
  } else {
    # Was unable to find 'file' or 'sample' column and cannot verify patient ID number.
    #print(paste0("Unable to identify a patient number inside file ", data.file, ". Please check that there is either a File or Sample column with the patient ID."))

    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.id.verify, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }

  # Add patient to list
  patients <- append(patients, patient)
    
  # Append drug compounds to list
  for (i in 1:nrow(drug.data)) {
    compounds <- append(compounds, as.character(drug.data[i,col.compound]))
  }
}

# Go through all csv files, gathering patient IDs and compounds present
for (data.file in filenames.csv) {
  
  # Open file
  drug.data <- read.csv(paste0(drug.dir, data.file))
  
  # Get column that has drug compounds, if exists; else give error
  if (length(which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))) > 0) {
    col.compound <- which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))[1]
  } else {
    #rint(paste0("No Compound column was found in ", data.file, ". Please reconcile."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.compound, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }
  # Get column that has either EC50 or IC50. Note, grabbing first instance since some have two sets, such
  #     as blood vs marrow samples. Also assuming that the EC50 value comes before the standard deviation 
  #     column that sometimes has EC50 in the name.
  if (length(which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))[1]
  } else {
    #print(paste0("No EC50 nor IC50 column was found in ", data.file, ". Please reconcile."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.EC50, filename=data.file, column=NA, id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }
  
  #### Check that the data is not NA for everything here. Give notice to user if it is with a reminder of
  ####     the need for static data. 
  # Check if the number of NAs is the same as number of drugs
  if (sum(is.na(suppressWarnings(as.numeric(drug.data[,col.EC50])))) == dim(drug.data)[1]) {
    # Give message and request that the user fix the data before restarting program
    #print(paste0("The ", colnames(drug.data)[col.EC50] , " column in the file ", data.file, " does not appear to have the necessary values. Please verify that the column has static, numeric values."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.static, filename=data.file, column=colnames(drug.data)[col.EC50], id1=NA, id2=NA))
    # Set clean flag to false
    clean <- FALSE
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
      #print(paste0("The patient ID number in the actual name of the file is ", patient, ", which is not the same as the patient ID number listed (", filename.patient, ") in the ", colnames(drug.data)[col.index], " column for file ", data.file, ". Please reconcile the issue."))
      
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(drug.data)[col.index], id1=patient, id2=filename.patient))
      # Set clean flag to false
      clean <- FALSE
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
      #print(paste0("The patient ID number in the actual name of the file is ", patient, ", which is not the same as the patient ID number listed (", filename.patient, ") in the ", colnames(drug.data)[col.index], " column for file ", data.file, ". Please reconcile the issue."))
      
      # Add filename and error code to errors table
      errors <- rbind(errors, data.frame(error=error.id, filename=data.file, column=colnames(drug.data)[col.index], id1=patient, id2=filename.patient))
      # Set clean flag to false
      clean <- FALSE
    } # else they match and the patient number should be fine
  } else {
    # Was unable to find 'file' or 'sample' column and cannot verify patient ID number.
    #print(paste0("Unable to identify a patient number inside file ", data.file, ". Please check that there is either a File or Sample column with the patient ID."))
    
    # Add filename and error code to errors table
    errors <- rbind(errors, data.frame(error=error.id.verify, filename=data.file, column=NA, id1=patient, id2=NA))
    # Set clean flag to false
    clean <- FALSE
  }
  
  # Add patient to list
  patients <- append(patients, patient)
  
  # Append drug compounds to list
  for (i in 1:nrow(drug.data)) {
    compounds <- append(compounds, as.character(drug.data[i,col.compound]))
  }
}

# Read in Archive files: Run1-4 and Run6; add patients and compounds to respective vectors.
#     These two files were originally in a challenging format to read in.
#     Each file was hand-fixed to have the EC50 values in a table with the header:
#     Compound | Parameters | AML1 | AML2 | AML3 |...
#     Files are currently held in a folder called "Archived" in test drug sensitivity directory

# Run1-4
drug.data <- read.csv(paste0(drug.dir, "Archived/AMLRun1-4_Static.csv"))
# Fix the Compound column because it reads in odd on this file...
dimnames(drug.data)[[2]][1] <- "Compound" 
# Get patient column indices
patient.columns <- which(str_detect(colnames(drug.data), coll("AML", ignore_case = TRUE)))
# Get patient ids from column names
patient.ids <- dimnames(drug.data)[[2]][patient.columns]
# Go through each patient; gather id # from patient id
for (id in patient.ids) {
  # Patient id # is always after AML; currently hardcoded as characters 4-6
  patient <- as.numeric(substr(id, 4, 6))
  # Add patient to list
  patients <- append(patients, patient)
}
# Get compounds and add to list
compounds <- append(compounds, drug.data$Compound)

#Run6
drug.data <- read.csv(paste0(drug.dir, "Archived/AMLRun6_Static.csv"))
# Get patient column indices
patient.columns <- which(str_detect(colnames(drug.data), coll("AML", ignore_case = TRUE)))
# Get patient ids from column names
patient.ids <- dimnames(drug.data)[[2]][patient.columns]
# Go through each patient; gather id # from patient id
for (id in patient.ids) {
  # Patient id # is always after AML; currently hardcoded as characters 4-6
  patient <- as.numeric(substr(id, 4, 6))
  # Add patient to list
  patients <- append(patients, patient)
}
# Get compounds and add to list
compounds <- append(compounds, drug.data$Compound)

#### Clean up of patient and compound vectors
# There may be duplicate patient numbers (duplicate files, etc); get unique set
patients.all <- unique(patients)
patients.all <- sort(patients.all)
write.table(patients.all, file = paste0(output.dir, "Patients_DrugSensitivity.txt"), row.names = FALSE, col.names = FALSE, sep = ",")
#### Pop up patient list with the ability to remove patients from it.

# Get unique set of compounds 
# Ignore NA compounds and case
compounds.all <- unique(compounds, ignore_case=TRUE)
# Remove all the excess numbered "compounds" generated by functions opening files.
#    Suppress warnings about NAs since that is the purpose of this: see if the compound
#    name can be converted to numeric without turning into an NA. If so, ignore it; if not,
#    keep it since it has letters (or at least something other than numbers).
compounds.all <- compounds.all[is.na(suppressWarnings(as.numeric(compounds.all)))]

#### Do initial synonym file stuff here
# Check if the file exists first; if not, make it.
if (file.exists(paste0(output.dir, "Compound_Synonyms.csv"))) {
  # Open synonym file
  synonyms <- read.csv(paste0(output.dir, "Compound_Synonyms.csv"), header = FALSE, na.strings = c("", " ", NA), stringsAsFactors = FALSE)
  
  # Check if the file was empty
  if (length(synonyms) > 0) {
    # Add all compounds in first column of synonym list
    #    This allows for names that users want to be shorter but never appear in compound list in files
    apply(synonyms, 1, function (x) {compounds.all <- append(compounds.all, x[1])})
    # Get the unique set in case synonym file added a double
    compounds.all <- unique(compounds, ignore_case=TRUE)
    
    # Remove all synonym compounds; assume first column in synonyms file is primary name
    # Go through each row, ignoring first column in the row and any NA values
    for (i in 1:nrow(synonyms)) {
      # Go through each synonym and check the location of the name in the full compounds set
      for (synonym in synonyms[i,-1]) {
        # Check if the synonym is NA
        if (!is.na(synonym)) {
          # Get index for where name appears in full compound set
          syn.index <- grep(synonym, compounds.all, ignore.case = TRUE)
          # Check that an appropriate index was actually found (i.e. name exists in list)
          if (length(syn.index) > 0) {
            # Remove from list
            compounds.all <- compounds.all[-syn.index]
          }  
        }
      }
    }
  } # else file is empty and move on
} else {
  # Create blank file
  file.create(paste0(output.dir, "Compound_Synonyms.csv"))
}

### Do initial compound removal here. This is for compounds that should not be included 
###    in the dataset at all.
# Check that the file exists first; if not, make it
if (file.exists(paste0(output.dir, "Compounds_Remove.txt"))) {
  remove.compounds <- read.csv(paste0(output.dir,"Compounds_Remove.txt"), header = FALSE, stringsAsFactors = FALSE)
  # Check if the file is empty
  if (length(remove.compounds) > 0) {
    # Make into a list instead of a table; technically not necessary, but seems less...odd
    remove.compounds <- remove.compounds[,1]
    # Remove each compound column in list 
    for (drug in remove.compounds) {
      # Get index for drug in full dataset, if exists
      drug.index <- which(tolower(colnames(compounds.all)) == tolower(drug))
      # Check that the drug existed in dataset
      if (length(drug.index > 0)) {
        # Remove the drug from dataset via column index
        compounds.all <- compounds.all[-drug.index]
      } # else compound wasn't found for some reason and nothing should be done
    }
  } # else file is empty and move on
} else {
  # Create blank file
  file.create(paste0(output.dir, "Compounds_Remove.txt"))
}

# Sort
compounds.all <- sort(compounds.all)

# Print compound list to text; this stores the data and we could have it pop up for viewing...
write(compounds.all, file = paste0(output.dir, "Compounds.txt"))

#### Give compound list to check. 
#### Have them update the synonym file.
#### Have them update the removal file.

#### Do second stage synonym file stuff here
# Open saved and (hopefully) current and updated synonym file
# Note: if no synonyms exist, then the file will be empty at this point
synonyms <- read.csv(paste0(output.dir, "Compound_Synonyms.csv"), header = FALSE, na.strings = c("", " ", NA), stringsAsFactors = FALSE)
# Check if file was empty
if (length(synonyms) > 0) {
  # Remove all synonym compounds; assume first column in synonyms file is primary name
  # Go through each row, ignoring first column in the row and any NA values
  for (i in 1:nrow(synonyms)) {
    # Go through each synonym and check the location of the name in the full compounds set
    for (synonym in synonyms[i,-1]) {
      # Check if the synonym is NA
      if (!is.na(synonym)) {
        # Get index for where name appears in full compound set
        syn.index <- grep(synonym, compounds.all, ignore.case = TRUE)
        # Check that an appropriate index was actually found (i.e. name exists in list)
        if (length(syn.index) > 0) {
          # Remove from list
          compounds.all <- compounds.all[-syn.index]
        }  
      }
    }
  }
} # else file is empty and move on

### Do second stage compound removal here. This is for compounds that should not be included 
###    in the dataset at all.
# Note: if no compounds to remove, file may be empty
remove.compounds <- read.csv(paste0(output.dir,"Compounds_Remove.txt"), header = FALSE, stringsAsFactors = FALSE)
# Check if file was empty first
if (length(remove.compounds) > 0) {
  # Make into a list instead of a table; technically not necessary, but seems less...odd
  remove.compounds <- remove.compounds[,1]
  # Remove each compound column in list 
  for (drug in remove.compounds) {
    # Get index for drug in full dataset, if exists
    drug.index <- which(tolower(colnames(compounds.all)) == tolower(drug))
    # Check that the drug existed in dataset
    if (length(drug.index > 0)) {
      # Remove the drug from dataset via column index
      compounds.all <- compounds.all[-drug.index]
    } # else compound wasn't found for some reason and nothing should be done
  }
} # else file is empty and move on

# Print compound list to text for storage, if needed 
write(compounds.all, file = paste0(output.dir, "Compounds.txt"))

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
  err.ind <- which(errors[,1] == error.sheet)
  if (length(err.ind) > 0) {
    sheet.errors.message <- "\nNo drug sensitivity was found in the following Excel files. Please ensure a sheet/tab with the data has one of the following phrases in the name: Fitted Parameters Static, FitParameters, combined, Fitted Parameters, blast param statsort, or 4Database."
    
    sheet.errors.list <- ""
    for (i in 1:length(err.ind)) {
      sheet.errors.list <- paste(sheet.errors.list, errors[err.ind[i],2], sep = "\n")
    }
    #sheet.errors.list <- paste(errors[err.ind,2], collapse = "\n")  
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
    #compound.errors.list <- paste(errors[err.ind, 2], collapse = "\n")  
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
    #EC50.errors.list <- paste(errors[err.ind, 2], collapse = "\n")  
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
    #static.errors.list <- paste(errors[err.ind, 2], "(", errors[err.ind, 3], ")", collapse = "\n") 
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
    #id.errors.list <- paste(errors[err.ind, 2], "- ID from filename: ", errors[err.ind, 4], "; ID from column ", errors[err.ind, 3], ": ", errors[err.ind, 5])  
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
    #id.verify.errors.list <- paste(errors[err.ind, 2], sep = "\n")
    string <- paste(string, id.verify.errors.message, id.verify.errors.list, sep = "\n")
  }
  
  cat(string, file = paste0(output.dir, "Drug Sensitivity File Errors.txt"), append = FALSE)
}