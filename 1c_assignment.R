#### Try It Yourself ####
# Practice Exercises 

# ----------------------------------------------------------------------------------------------------------------

# 1. Check Cholesterol level (using if) 
# Write an If statement to check cholesterol level is greater than 240, 
# if true, it will prints “High Cholesterol”

cholesterol <- 230
if(cholesterol > 240){
  print("High Cholesterol")
}

# ----------------------------------------------------------------------------------------------------------------

# 2. Blood Pressure Status (using if...else)
# Write an if…else statement to check if blood pressure is normal.
# If it’s less than 120, print: “Blood Pressure is normal”
# If false then print: “Blood Pressure is high”

Systolic_bp <- 130
if(Systolic_bp < 120){
  print("Blood Pressure is normal")
} else{
  print("Blood Pressure is high")
}

# ----------------------------------------------------------------------------------------------------------------

# 3. Automating Data Type Conversion with for loop

# Use patient_info.csv data and metadata.csv
# Perform the following steps separately on each dataset (patient_info.csv data and metadata.csv)
# Create a copy of the dataset to work on.
patient_info <- read.csv(file.choose())
patient_info_copy <- patient_info
metadata <- read.csv(file.choose())
metadata_copy <- metadata

# Identify all columns that should be converted to factor type.
# height and gender columns in metadata set
# gender, diagnosis and smoker columns in patient_info data set 

# Store their names in a variable (factor_cols).

factor_cols_patient_info <- c("gender", "diagnosis", "smoker")

factor_cols_metadata <- c("height", "gender")


# Use a for loop to convert all the columns in factor_cols to factor type.
# Pass factor_cols to the loop as a vector.
for (col in factor_cols_patient_info){
  patient_info_copy[[col]] <- as.factor(patient_info_copy[[col]])
}

str(patient_info_copy)
str(patient_info)

for (col in factor_cols_metadata){
  metadata_copy[[col]] <- as.factor(metadata_copy[[col]])
}
str(metadata_copy)
str(metadata)

# ----------------------------------------------------------------------------------------------------------------

# 4. Converting Factors to Numeric Codes

# Choose one or more factor columns (e.g., smoking_status).
# Convert "Yes" to 1 and "No" to 0 using a for loop.

# Hint:
# binary_cols <- c("smoking_status")   # store column names in a vector
# use ifelse() condition inside the loop to replace Yes with 1 and No with 0
# for (col in binary_cols) {
#   data[[col]] <- # insert your ifelse() code here
# }


binary_cols_patient <- c("smoker") 
for (col in binary_cols_patient) {
  binary_cols_patient[[col]] <- ifelse(binary_cols_patient[[col]] == "Yes", 1,
                                       ifelse(binary_cols_patient[[col]] == "No", 0, NA))
  str(binary_cols_patient)
  
  