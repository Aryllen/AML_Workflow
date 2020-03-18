##### 
# Drug Sensitivity and Variant Anyalysis
# 
# Author: Nicole Kauer
# Date: August 5, 2019
#
# Purpose: Read in variants and drug sensitivity data for patients. Filter based on
#     inputs from user, and do analysis.
# 
# Execute from terminal:
#     'Rscript Variant_Source_to_CSV.r -v <variant file directory in quotes>
#         + -d <drug sensitivity file directory in quotes>' 
#      Include any data filtering options.
#
# Input(s):
#     Directories
#         -v --var.dir: directory for variant files
#         -d --drug.dir: directory for drug sensitivity files
#
#     Variant Data Filters
#         -f --frequency: minimum value of frequency considered relevant for all variants (default = 2.5%)
#         -s* SIFT: logical options to include for missense SIFT types: 
#             -sd --sift.deleterious (default = TRUE)
#             -st --sift.tolerated (default = FALSE)
#         -p* Polyphen: logical options to include for missense Polyphen types: 
#             -pd --poly.prob.damaging (default = TRUE)
#             -ppd --poly.poss.damaging (default = TRUE)
#             -pb --poly.benign (default = FALSE)
#
#     Drug Sensitivity Thresholds
#         -a --active.maximum: maximum value of EC50 to be considered active; anything above is "inactive"
#         -m --active.minimum: mininum value of EC50 to be considered highly active; anything below is considered highly active
#
#     Subsets
#         -p --patient.subset: logical option to allow for editing of the patient subset
#         -c --compound.subset: logical option to allow for editing of the compound subset
#
# Output(s): 
#####

# Package(s) needed; will install if not present
# For command line arguments
# if(!require("optparse", quietly=TRUE)) {
#   install.packages("optparse", quiet=TRUE, repos='http://cran.us.r-project.org')
# }
library("optparse", quietly=TRUE)
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

# Command line options
options <- list(
  make_option(c("-v", "--var.dir"), type="character", default=NULL, help="Directory for variant files."),
  make_option(c("-d", "--drug.dir"), type="character", default=NULL, help="Directory for drug sensitivity files."),
  make_option(c("-o", "--output.dir"), type="character", default=NULL, help="Directory to output files."),
  make_option(c("-f", "--frequency"), type="double", default=2.5, help="Minimum variant frequency to be included in dataset."),
  make_option(c("--sift.deleterious"), type="logical", default=TRUE, help="Include variants with SIFT = deleterious."),
  make_option(c("-st", "--sift.tolerated"), type="logical", default=FALSE, help="Include missense variants with SIFT = tolerated."),
  make_option(c("--poly.prob.damaging"), type="logical", default=TRUE, help="Include missense variants with Polyphen = probably damaging."),
  make_option(c("--poly.poss.damaging"), type="logical", default=TRUE, help="Include missense variants with Polyphen = possibly damaging."),
  make_option(c("-pb", "--poly.benign"), type="logical", default=FALSE, help="Include missense variants with Polyphen = benign."),
  make_option(c("-a", "--active.maximum"), type="double", default=1E-6, help="Maximum EC50 value to be considered active; all values greater are considered inactive."),
  make_option(c("-m", "--active.minimum"), type="double", default=1E-12, help="Minimum EC50 value to be considered highly active; all values smaller are considered highly active."),
  make_option(c("-p", "--patient.subset"), type="logical", default=FALSE, help="Choose the patient subset, else the previously chosen subset is used."),
  make_option(c("-c", "--compound.subset"), type="logical", defaul=FALSE, help="Choose the drug compound subset, else the previously chosen subset is used.")
)
# Argument retriever
opt_parser <- OptionParser(option_list=options)
# Get arguments
opts <- parse_args(opt_parser)

# Need to have the directories input; verify that they were.
# Note: not checking that they are valid inputs; also, if I can just move the names along from the previous widget in 
# BwB, then this should be unnecessary to check.
if (is.null(opts$var.dir)) {
  print("No folder for the variant files was detected. Please input the directory in the command line option '-v <directory>'.")
} else {
  # Append / to end because the code needs it for the directory to map correctly to files
  opts$var.dir <- paste0(opts$var.dir, "/")
}
if (is.null(opts$drug.dir)) {
  print("No folder for the drug sensitivity files was detected. Please input the directory in the command line option '-d <directory>'.")
} else {
  # Append / to end because the code needs it for the directory to map correctly to files
  opts$drug.dir <- paste0(opts$drug.dir, "/")
}
if (is.null(opts$output.dir)) {
  print("No folder for the output files was detected. Please input the directory in the command line option '-v <directory>'.")
} else {
  # Append / to end because the code needs it for the directory to map correctly to files
  opts$output.dir <- paste0(opts$output.dir, "/")
}

# Output directories for variant and drug sensitivity data
output.var.dir <- paste0(opts$var.dir,"Output/")
output.drug.dir <- paste0(opts$drug.dir,"Output/")

# Output directories for variant and drug sensitivity data
#output.var.dir <- paste0(opts$var.dir, "Output/")
#output.drug.dir <- paste0(opts$drug.dir, "Output/")
# Output directory for combined files made in this program
# This directory needs to be made, if it doesn't exist. Also not sure how to put it somewhere specific.
# I've been putting it in the same directory as the folders for variants and drug sensitivity, but 
# would need to alter the string to put it there.
#output.dir <- "...\\Output\\"
#if(!dir.exists(output.dir)) {
#  dir.create(output.dir)
#}

### FOR TESTING ONLY ###
#opts$var.dir <- "C:\\Users\\Nicole\\Documents\\AML\\PamBeckerData\\Variants\\"
#opts$drug.dir <- "C:\\Users\\Nicole\\Documents\\AML\\PamBeckerData\\DrugSensitivity\\"
#opts$output.dir <- "C:\\Users\\Nicole\\Documents\\AML\\PamBeckerData\\Output\\"
#output.var.dir <- "C:\\Users\\Nicole\\Documents\\AML\\PamBeckerData\\Variants\\Output\\"
#output.drug.dir <- "C:\\Users\\Nicole\\Documents\\AML\\PamBeckerData\\DrugSensitivity\\Output\\"


##### Drug Sensitivity Data Gathering #####

# Get lists of the files in the directories
filenames.drug.excel <- list.files(opts$drug.dir, pattern = "*.xlsx")
filenames.drug.csv <- list.files(opts$drug.dir, pattern = "*.csv")

# Patient list from drug sensitivity cleaning
patients.drug <- as.numeric(readLines(paste0(output.drug.dir, "Patients_DrugSensitivity.txt")))
# Compound list
compounds <- readLines(paste0(output.drug.dir, "Compounds.txt"))

# Synonyms from file
synonyms <- read.csv(paste0(output.drug.dir, "Compound_Synonyms.csv"), header = FALSE, na.strings = c("", " ", NA), stringsAsFactors = FALSE)

##### Function(s) #####
# Synonym row index function
# Input: drug.name (string) - name of drug
# Returns: row index for matching drug name in synonym file or 0 if not found
syns.row <- function (drug.name) {
  # Initialize to default return value if drug name not found
  match.index <- 0
  # Go through each row synonym file
  for (i in 1:nrow(synonyms)) {
    # Check if drug name is in the row; save index if found
    if (drug.name %in% t(tolower(synonyms[i,]))) {
      match.index <- i
    }
  }
  return(match.index)
}
#####################

# Dataframe for the collected drug sensitivity data
data.drug <- data.frame(matrix(nrow=length(patients.drug), ncol=length(compounds)))
rownames(data.drug) <- patients.drug
colnames(data.drug) <- compounds

# Go through all Excel files and start collecting the drug sensitivity data              
for (data.file in filenames.drug.excel) {
  
  # Get names of sheets in Excel workbook
  sheets <- excel_sheets(paste0(opts$drug.dir,data.file))
  
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
    #### Need to have the program tell the user that they did not find a sheet in the file for
    ####     the data. Have a reminder of what is expected/searched for and quit the program so
    ####     they can fix and restart.
    print(paste0("The drug sensitivity data was not found in ", data.file, ". Please run the DrugSensitivity_CleanVerification widget to find any errors that should be fixed."))
  }
  
  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  patient.index <- which(rownames(data.drug) == patient)
  
  # Open file; suppress unimportant messages, again
  drug.data <- suppressMessages(read_excel(paste0(opts$drug.dir, data.file), sheet=sheets[sensitivity.sheet]))
  
  # Get column that has drug compounds
  if (length(which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))) > 0) {
    col.compound <- which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))[1]
  } else {
    col.compound <- 0
    print(paste0("No Compound column was found in ", data.file, ". Please reconcile. For now, this patient will be skipped."))
  }
  
  # Get column that has either EC50 or IC50. Note, grabbing first instance since some have two sets, such
  #     as blood vs marrow samples. Also assuming that the EC50 value comes before the standard deviation 
  #     column that sometimes has EC50 in the name.
  if (length(which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))[1]
  } else {
    col.EC50 <- 0
    print(paste0("No EC50 nor IC50 column was found in ", data.file, ". Please reconcile. For now, this patient will be skipped"))
  }
  
  # Gather the EC50 values if both they and the compound names exist in the file
  if (col.EC50 > 0 && col.compound > 0) {
    for (drug in drug.data[[col.compound]]){
      # Check to see if drug name is in compound list; if not, it's either:
      #     a) a synonym, or
      #     b) a removed compound
      #     Need to treat each category accordingly
      if (tolower(drug) %in% tolower(colnames(data.drug))) { # Drug name is in full compound dataset
        # Get index for drug compound in the full dataset
        data.drug.index <- which(tolower(colnames(data.drug)) == tolower(drug))
        # Get index for drug compound in the current patient's data
        # Note: only using value for first drug compound of the same name; some patient's have
        # duplicates, but one set will be blasts and one is something else
        drug.data.index <- which(drug.data[,col.compound] == drug)[1]
        # Get EC50 value if the column is called that; else get the IC50 value
        data.drug[patient.index, data.drug.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index, col.EC50])))
      } else { # Drug name is not in full compound dataset
        # Try to get the row index for the drug from synonym file; will be 0 if no synonyms
        syn.index <- syns.row(tolower(drug))
        if (syn.index > 0) { # Must be a synonym drug
          # Preferred drug name in full set is the first name in synonym row
          # Get index for drug compound in the full dataset
          data.drug.index <- which(tolower(colnames(data.drug)) == tolower(synonyms[syn.index, 1]))
          drug.data.index <- which(drug.data[,col.compound] == drug)[1]
          # Get EC50 value if the column is called that; else get the IC50 value
          data.drug[patient.index, data.drug.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index,col.EC50])))
        } # else must be a removed drug and should do nothing with it
      }
    }
  } # else wasn't able to find the EC50 column and the patient data gathering was skipped
}

# Go through all csv files and start collecting EC50 data
for (data.file in filenames.drug.csv) {
  
  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  patient.index <- which(rownames(data.drug) == patient)
  
  # Open file
  drug.data <- read.csv(paste0(opts$drug.dir, data.file))
  #print(data.file)
  
  # Get column that has drug compounds, if exists; else give error
  if (length(which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))) > 0) {
    col.compound <- which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))[1]
  } else {
    col.compound <- 0
    print(paste0("No Compound column was found in ", data.file, ". Please reconcile. For now, this patient will be skipped."))
  }
  
  # Get column that has either EC50 or IC50. Note, grabbing first instance since some have two sets, such
  #     as blood vs marrow samples. Also assuming that the EC50 value comes before the standard deviation 
  #     column that sometimes has EC50 in the name.
  if (length(which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))[1]
  } else {
    col.EC50 <- 0
    print(paste0("No EC50 nor IC50 column was found in ", data.file, ". Please reconcile. For now, this patient will be skipped."))
  }
  
  # Gather the EC50 values if both they and the compound names exist in the file
  if (col.EC50 > 0 && col.compound > 0) {
    for (drug in drug.data[[col.compound]]){
      # Check to see if drug name is in compound list; if not, it's either:
      #     a) a synonym, or
      #     b) a removed compound
      #     Need to treat each category accordingly
      if (tolower(drug) %in% tolower(colnames(data.drug))) { # Drug name is in full compound dataset
        # Get index for drug compound in the full dataset
        data.drug.index <- which(tolower(colnames(data.drug)) == tolower(drug))
        # Get index for drug compound in the current patient's data
        # Note: only using value for first drug compound of the same name; some patient's have
        # duplicates, but one set will be blasts and one is something else
        drug.data.index <- which(drug.data[,col.compound] == drug)[1]
        # Get EC50 value if the column is called that; else get the IC50 value
        data.drug[patient.index, data.drug.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index, col.EC50])))
      } else { # Drug name is not in full compound dataset
        # Try to get the row index for the drug from synonym file; will be 0 if no synonyms
        syn.index <- syns.row(tolower(drug))
        if (syn.index > 0) { # Must be a synonym drug
          # Preferred drug name in full set is the first name in synonym row
          # Get index for drug compound in the full dataset
          data.drug.index <- which(tolower(colnames(data.drug)) == tolower(synonyms[syn.index, 1]))
          drug.data.index <- which(drug.data[,col.compound] == drug)[1]
          # Get EC50 value if the column is called that; else get the IC50 value
          data.drug[patient.index, data.drug.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index,col.EC50])))
        } # else must be a removed drug and should do nothing with it
      }
    }
  } # else wasn't able to find the EC50 column and the patient data gathering was skipped
}

# Read in Archive files: Run1-4 and Run6; add data to table.
#     These two files were originally in a challenging format to read in.
#     Each file was hand-fixed to have the EC50 values in a table with the header:
#     Compound | Parameters | AML1 | AML2 | AML3 |...
#     Files are currently held in a folder called "Archived" in test drug sensitivity directory
# Run1-4
drug.data <- read.csv(paste0(opts$drug.dir, "Archived/AMLRun1-4_Static.csv"))
# Fix the Compound column because it reads in odd on this file...
dimnames(drug.data)[[2]][1] <- "Compound"
# Get patient column indices
patient.columns <- which(str_detect(colnames(drug.data), coll("AML", ignore_case = TRUE)))
# Get patient ids from column names
patient.ids <- dimnames(drug.data)[[2]][patient.columns]
# Vector for holding compounds and EC50 values before adding to full dataset
patient.data <- vector()
# Go through each patient; gather id # from patient id, and add EC50 values to full dataset.
for (id in patient.ids) {
  # Get patient id # 
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  # Get index of patient in Run1-4 dataset
  patient.Run.index <- which(colnames(drug.data) == id)
  # Get index of patient in full dataset
  patient.index <- which(rownames(data.drug) == patient)
  
  # Gather the EC50 values
  for (drug in drug.data$Compound){
    # Check to see if drug name is in compound list; if not, it's either:
    #     a) a synonym, or
    #     b) a removed compound
    #     Need to treat each category accordingly
    if (tolower(drug) %in% tolower(colnames(data.drug))) { # Drug name is in full compound dataset
      # Get index for drug compound in the full dataset
      data.drug.index <- which(tolower(colnames(data.drug)) == tolower(drug))
      # Get index for drug compound in the current patient's data
      # Note: only using value for first drug compound of the same name; some patient's have
      # duplicates, but one set will be blasts and one is something else
      drug.data.index <- which(drug.data$Compound == drug)[1]
      # Get EC50 value
      #     Note that the file reads in with the first patient column as a factor, instead of numeric.
      #     Fix by converting into numeric. Cannot simply use as.numeric or gets the wrong value.
      #     Have to check if the value is a factor or numeric. Alternative would be to fix column at start.
      if (is.factor(drug.data[drug.data.index, patient.Run.index])) {
        # Suppress messages about NAs
        data.drug[patient.index, data.drug.index] <- suppressMessages(suppressWarnings(as.numeric(levels(drug.data[drug.data.index, patient.Run.index])[drug.data[drug.data.index, patient.Run.index]])))  
      } else {
        data.drug[patient.index, data.drug.index] <- suppressMessages(suppressWarnings(as.numeric(drug.data[drug.data.index, patient.Run.index])))
      }
    } else { # Drug name is not in full compound dataset
      # Try to get the row index for the drug from synonym file; will be 0 if no synonyms
      syn.index <- syns.row(tolower(drug))
      if (syn.index > 0) { # Must be a synonym drug
        # Preferred drug name in full set is the first name in synonym row
        # Get index for drug compound in the full dataset
        data.drug.index <- which(tolower(colnames(data.drug)) == tolower(synonyms[syn.index, 1]))
        drug.data.index <- which(drug.data$Compound == drug)[1]
        # Get EC50 value
        #     Note that the file reads in with the first patient column as a factor, instead of numeric.
        #     Fix by converting into numeric. Cannot simply use as.numeric or gets the wrong value.
        #     Have to check if the value is a factor or numeric. Alternative would be to fix column at start.
        if (is.factor(drug.data[drug.data.index, patient.Run.index])) {
          # Suppress messages about NAs
          data.drug[patient.index, data.drug.index] <- suppressMessages(suppressWarnings(as.numeric(levels(drug.data[drug.data.index, patient.Run.index])[drug.data[drug.data.index, patient.Run.index]])))  
        } else {
          data.drug[patient.index, data.drug.index] <- suppressMessages(suppressWarnings(as.numeric(drug.data[drug.data.index, patient.Run.index])))
        }
      } # else must be a removed drug and should do nothing with it
    }
  }
}

# Run6
drug.data <- read.csv(paste0(opts$drug.dir, "Archived/AMLRun6_Static.csv"))
# Get patient column indices
patient.columns <- which(str_detect(colnames(drug.data), coll("AML", ignore_case = TRUE)))
# Get patient ids from column names
patient.ids <- dimnames(drug.data)[[2]][patient.columns]
# Vector for holding compounds and EC50 values before adding to full dataset
patient.data <- vector()
# Go through each patient; gather id # from patient id, and add EC50 values to full dataset.
for (id in patient.ids) {
  # Get patient id # 
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  # Get index of patient in Run1-4 dataset
  patient.Run.index <- which(colnames(drug.data) == id)
  # Get index of patient in full dataset
  patient.index <- which(rownames(data.drug) == patient)
  
  # Gather the EC50 values
  for (drug in drug.data$Compound){
    # Check to see if drug name is in compound list; if not, it's either:
    #     a) a synonym, or
    #     b) a removed compound
    #     Need to treat each category accordingly
    if (tolower(drug) %in% tolower(colnames(data.drug))) { # Drug name is in full compound dataset
      # Get index for drug compound in the full dataset
      data.drug.index <- which(tolower(colnames(data.drug)) == tolower(drug))
      # Get index for drug compound in the current patient's data
      # Note: only using value for first drug compound of the same name; some patient's have
      # duplicates, but one set will be blasts and one is something else
      drug.data.index <- which(drug.data$Compound == drug)[1]
      # Get EC50 value
      #     Note that the file reads in with the first patient column as a factor, instead of numeric.
      #     Fix by converting into numeric. Cannot simply use as.numeric or gets the wrong value.
      #     Have to check if the value is a factor or numeric. Alternative would be to fix column at start.
      if (is.factor(drug.data[drug.data.index, patient.Run.index])) {
        # Suppress messages about NAs
        data.drug[patient.index, data.drug.index] <- suppressMessages(suppressWarnings(as.numeric(levels(drug.data[drug.data.index, patient.Run.index])[drug.data[drug.data.index, patient.Run.index]])))  
      } else {
        data.drug[patient.index, data.drug.index] <- suppressMessages(suppressWarnings(as.numeric(drug.data[drug.data.index, patient.Run.index])))
      }
    } else { # Drug name is not in full compound dataset
      # Try to get the row index for the drug from synonym file; will be 0 if no synonyms
      syn.index <- syns.row(tolower(drug))
      if (syn.index > 0) { # Must be a synonym drug
        # Preferred drug name in full set is the first name in synonym row
        # Get index for drug compound in the full dataset
        data.drug.index <- which(tolower(colnames(data.drug)) == tolower(synonyms[syn.index, 1]))
        drug.data.index <- which(drug.data$Compound == drug)[1]
        # Get EC50 value
        #     Note that the file reads in with the first patient column as a factor, instead of numeric.
        #     Fix by converting into numeric. Cannot simply use as.numeric or gets the wrong value.
        #     Have to check if the value is a factor or numeric. Alternative would be to fix column at start.
        if (is.factor(drug.data[drug.data.index, patient.Run.index])) {
          # Suppress messages about NAs
          data.drug[patient.index, data.drug.index] <- suppressMessages(suppressWarnings(as.numeric(levels(drug.data[drug.data.index, patient.Run.index])[drug.data[drug.data.index, patient.Run.index]])))  
        } else {
          data.drug[patient.index, data.drug.index] <- suppressMessages(suppressWarnings(as.numeric(drug.data[drug.data.index, patient.Run.index])))
        }
      } # else must be a removed drug and should do nothing with it
    }
  }
}

##### Drug Sensitivity Data Thresholds & Transform #####

# Set thresholds on EC50 values 
# Any value > active.maximum is set to active.maximum/1e-3 to make large enough to be different from barely active drugs
# Any value < active.minimum is set to active.minimum
# Upper bound threshold to set inactive values to
drug.upper.threshold <- opts$active.maximum/1e-3

# Set all threshold values, if needed
for (patient in rownames(data.drug)) {
  for (drug in colnames(data.drug)) {
    # Make sure not NA
    if (!is.na(data.drug[patient, drug])) {
      # Check if value is higher than active threshold and set to upper threshold value, if so
      if (data.drug[patient, drug] > opts$active.maximum) {
        data.drug[patient, drug] <- drug.upper.threshold
      } else if (data.drug[patient, drug] < opts$active.minimum) {
        data.drug[patient, drug] <- opts$active.minimum
      } 
    } # else leave as NA
  }
}

# Remove columns (and rows, just in case) that are all NA
drugs.NA.col <- apply(data.drug, 2, function (x) sum(is.na(x)))
drugs.NA.row <- apply(data.drug, 1, function (x) sum(is.na(x)))
data.drug <- data.drug[-which(drugs.NA.row == dim(data.drug)[2]), -which(drugs.NA.col == dim(data.drug)[1])]

# -log(EC50) of all EC50 values
data.drug <- -log(data.drug)

##### Variant Data Gathering #####

# Get list of the files in the directory
filenames.var <- list.files(opts$var.dir, pattern = "*.xlsx")

# Patient list from variant cleaning
patients.var <- as.numeric(readLines(paste0(output.var.dir, "Patients_Variants.txt")))
# Vector for storing gene list
genes <- readLines(paste0(output.var.dir, "Genes.txt"))

# Dataframe for the collected variant data
data.var <- data.frame(matrix(nrow=length(patients.var), ncol=length(genes)))
rownames(data.var) <- patients.var
colnames(data.var) <- genes
# Initialize all values to 0
data.var[is.na(data.var)] <- 0

# Go through all files again and start collecting the variant data              
for (data.file in filenames.var) {

  # Get names of sheets in Excel workbook
  sheets <- excel_sheets(paste0(opts$var.dir,data.file))
  
  sheet.miss <- which(str_detect(sheets, coll("missense", ignore_case = TRUE)))
  sheet.indel <- which(str_detect(sheets, coll("indel", ignore_case = TRUE))) 

  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  patient.index <- which(rownames(data.var) == patient)
  
  # Open file to missense info
  var.data <- read_excel(paste0(opts$var.dir, data.file), sheet=sheets[sheet.miss])
  
  # Filter based on startup arguments
  # SIFT filtering
  sift.ind <- vector()
  if (opts$sift.deleterious) {
    # Get indices for deleterious SIFT
    sift.ind <- c(sift.ind, grep("delet", var.data$SIFT))
  } # Else not chosen and nothing added to sift.ind
  if (opts$sift.tolerated) {
    # Get indices for tolerated SIFT
    sift.ind <- c(sift.ind, grep("toler", var.data$SIFT))
  } # Else not chosen and nothing added to sift.ind
  # Only keep indices gathered in filters above
  # To be safe, check that sift.ind is not empty; this could happen if neither option is checked
  if (length(sift.ind) > 0) {
    var.data <- var.data[sift.ind,]
  } # else keep all data
  
  # Polyphen filtering
  poly.ind <- vector()
  if (opts$poly.prob.damaging) {
    # Get indices for probably damaging polyphen
    poly.ind <- c(poly.ind, grep("prob", var.data$Polyphen))
  } # Else not chosen and nothing added to poly.ind
  if (opts$poly.poss.damaging) {
    # Get indices for possibly damaging polyphen
    poly.ind <- c(poly.ind, grep("poss", var.data$Polyphen))
  } # Else not chosen and nothing added to poly.ind
  if (opts$poly.benign) {
    # Get indices for benign polyphen
    poly.ind <- c(poly.ind, grep("benign", var.data$Polyphen))
  }
  # Only keep indices gathered in filters above
  # To be safe, check that poly.ind is not empty; this could happen if no options are checked
  if (length(poly.ind) > 0) {
    var.data <- var.data[poly.ind,]
  } # else keep all data
  
  # Remove missense with freq < opts$Frequency (default is 2.5%)
  # Note that currently expect value in %, not decimal so / by 100
  var.data <- subset(var.data, var.data$Frequency > (opts$frequency/100) )
  
  # Count number of variations in each gene and store
  for (gene in var.data$Gene){
    # Get index for gene in the full dataset
    data.var.index <- which(colnames(data.var) == gene)
    # Add 1 to current value
    data.var[patient.index, data.var.index] <- data.var[patient.index, data.var.index] + 1
  }
  
  # Open file to indel info
  var.data <- read_excel(paste0(opts$var.dir, data.file), sheet=sheets[sheet.indel])
  # Remove indel with freq < opts$Frequency (default is 2.5%)
  # Note that currently expect value in %, not decimal so / by 100
  var.data <- subset(var.data, var.data$Frequency > (opts$frequency/100) )
  
  # Count number of variations in each gene and store
  for (gene in var.data$Gene){
    # Get index for gene in the full dataset
    data.var.index <- which(colnames(data.var) == gene)
    # Add 1 to current value
    data.var[patient.index, data.var.index] <- data.var[patient.index, data.var.index] + 1
  }
}

##### Subsetting data by which patients have both drug and variant data
patient.subset <- intersect(rownames(data.drug), rownames(data.var))
# Assuming that the subset has at least one, but put a condition just in case.
if (length(patient.subset) > 0) {
  data.var <- data.var[patient.subset,]
  data.drug <- data.drug[patient.subset,] 
} else {
  print("There appears to be no shared patients in the dataset.")
  quit()
}

##### Subsetting data by patient list #####

# Need only the subset of patients that have both drug sensitivity and variant data
# Also only need the subset of patients that the user wants
# Check to see if Patient_Subset.csv has been created, yet
#    If yes, open and use only the patient subset with 1s; only ask to edit if user has requested option with startup arg
#    If no, get patient list that has both data sets, create file, and ask user to edit subset
if (!file.exists(paste0(opts$output.dir, "Patient_Subset.csv"))) {
  # File doesn't exist, yet
  # Need to know which patients have both drug sensitivity and variant data
  patients.all.data <- intersect(rownames(data.drug), rownames(data.var))
  # Create a dataframe for patient list and whether to include
  patients.include <- data.frame(matrix(nrow=length(patients.all.data), ncol=2))
  colnames(patients.include) <- c("Patient ID", "Include")
  patients.include$`Patient ID` <- patients.all.data
  # Default is to include all patients
  patients.include$Include <- 1
  # Write to file
  write.csv(patients.include, paste0(opts$output.dir, "Patient_Subset.csv"), row.names = FALSE)
}
if (opts$patient.subset) {
  # Let user update file
  #### Add in interactive bit!
}
# Read in file
patients.include <- read.csv(paste0(opts$output.dir, "Patient_Subset.csv"))
# Get only patient IDs for which Include == 1
patient.subset <- which(patients.include$Include == 1)
# Keep only the patient subset desired in each dataset
data.drug <- data.drug[patient.subset,]
data.var <- data.var[patient.subset,]

# Remove columns that are all NA (ie 0 since I made them all 0 earlier) 
#     since this subset of patients may not have data for certain drugs
drugs.NA.col <- apply(data.drug, 2, sum)
# If there is a drug that has no data, remove it; else move on
if (length(which(drugs.NA.col == 0)) > 0) {
  data.drug <- data.drug[, -which(drugs.NA.col == 0)]  
}

##### Subset data by drug compound list #####

# Only need the subset of compounds that the user wants
# Check to see if Compound_Subset.csv has been created, yet
#    If yes, open and use only the compound subset with 1s; only ask to edit if user has requested option with startup arg
#    If no, get compound list that is not NA for the subset of patients chosen
if (!file.exists(paste0(opts$output.dir, "Compound_Subset.csv"))) {
  # File doesn't exist, yet
  # Create a dataframe for compound list and whether to include
  compounds.include <- data.frame(matrix(nrow=dim(data.drug)[2], ncol=2))
  colnames(compounds.include) <- c("Compound", "Include")
  compounds.include$`Compound` <- colnames(data.drug)
  # Default is to include all compounds
  compounds.include$Include <- 1
  # Write to file
  write.csv(compounds.include, paste0(opts$output.dir, "Compound_Subset.csv"), row.names = FALSE)
}
if (opts$compound.subset) {
  # Let user update file
  #### Add in interactive bit!
}
# Read in file
compounds.include <- read.csv(paste0(opts$output.dir, "Compound_Subset.csv"))
# Get only compounds for which Include == 1
compound.subset <- which(compounds.include$Include == 1)
# Keep only the drug subset desired in each dataset
data.drug <- data.drug[,compound.subset]

##### Need to filter out compounds that have no useful data (all NA, all inactive drug values, etc)

# Filter out drugs with all NA
# Count the number of NAs in each column
drug.na.count <- apply(data.drug, 2, function (x) {sum(is.na(x))})
# Find which columns have the same number of NAs as patients in the dataset
drug.all.na <- which(drug.na.count == dim(data.drug)[1])
# Remove drugs that are all NA, if any
if (length(drug.all.na) > 0) {
  data.drug <- data.drug[,-drug.all.na]  
}

# Filter out genes with no variants
# Count the number of 0s in each column
var.0.count <- apply(data.var, 2, function(x){length(which(x == 0))})
# Find which columns have the same number of 0s as patients in the dataset
var.all.0 <- which(var.0.count == dim(data.var)[1])
# Remove genes that have no variants, if any
if (length(var.all.0) > 0) {
  data.var <- data.var[,-var.all.0]  
}

# Filter out drugs that have no active drug values
# Count the number of values greater than the active drug minimum
drug.inactive.count <- apply(data.drug, 2, function(x){length(which(x > opts$active.maximum))})
# Find which columns have the same number of inactive values as patients in dataset
drug.all.inactive <- which(drug.inactive.count == dim(data.drug)[1])
# Remove drugs that have no active sensitivities, if any
if (length(drug.all.inactive) > 0) {
  data.drug <- data.drug[,-drug.all.inactive]  
}

##### Write to files #####

# Output filtered, subset, etc tables to files
write.csv(data.drug, file=paste0(opts$output.dir, "PatientsVsDrugs.csv"))
write.csv(data.var, file=paste0(opts$output.dir, "PatientsVsVariants.csv"))
