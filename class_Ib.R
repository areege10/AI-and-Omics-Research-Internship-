getwd()
setwd("C:/Users/ASUS/Desktop")
getwd()
dir.create("AI_Omics_Internship_2025")
setwd("C:/Users/ASUS/Desktop/AI_Omics_Internship_2025")
getwd()
setwd("C:/Users/ASUS/Desktop/AI_Omics_Internship_2025/project folder")
getwd()
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("Tasks")
dir.create("plots")
patient_info <- read.csv(file.choose())
str(patient_info)
View(patient_info)

##### The following variables should ideally be converted from chr to factor: 
####gender,diagnosis and smoker
patient_info$gender <- as.factor(patient_info$gender)
str(patient_info)
patient_info$diagnosis <- as.factor(patient_info$diagnosis)
str(patient_info)
patient_info$smoker <- as.factor(patient_info$smoker)
str(patient_info)
patient_info$smoker <- factor(patient_info$smoker,
                              levels = c("Yes", "No"),
                              labels = c(1, 0))
str(patient_info)
View(patient_info)

write.csv(patient_info, file = "clean_data/patient_info_clean.csv"



