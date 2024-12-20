library(tidyverse)
library(googledrive)
library(readxl)
library(janitor)
library(magrittr)
library(stringr)
library(openxlsx)
library(dplyr)
library(OpenSpecy) #If need dev. devtools::install_github("wincowgerDEV/OpenSpecy-package")
library(googlesheets4)
library(data.table)
gs4_auth()

# Ask User the name of the folder they want to access
project_name <- readline(prompt = "Enter the Name of the Project: ")

# Read in all functions
lapply(c("R/analyze_features.R", "R/naorempty.R", "R/get_lib.R", "R/determine_thresholds.R", "R/particle_image.R", "R/library_prep.R"), source)

# Load library for analysis
if(is(tryCatch(check_lib(c("derivative")),
               error=function(e) e, warning=function(w) w), 
      "warning")){
  get_lib("derivative")
}

lib <- OpenSpecy::load_lib("derivative") |> 
  filter_spec(lib$metadata$spectrum_type == "ftir")

# Find folder names present in project folder - heavily dependent if it contains _ hyphen
folder_find <- shared_drive_find("Project") |>
  drive_ls("Customer Projects") |>
  drive_ls() |>
  filter(name %in% project_name) |>
  drive_ls() #|>
  # filter(name %in% project_name) |>
  # drive_ls()


# Determine whether to conduct individual files (.SPA) or multi file (.hdr/.dat/.jpg) or both analysis
if (any(grepl("Particle_500um", folder_find$name))) {
  source("R/individual_files.R")
} else if (any(grepl("Export_Files", folder_find$name))) {
  source("R/multi_files.R")
} else if (any(grepl("Particle_500um", folder_find$name)) | any(grepl("Export_Files", folder_find$name))) {
  source("R/individual_files.R")
  source("R/multi_files.R")
}

# Cleanup Data ----
aux_sheets <- list.files(path = "data", pattern = ".xlsx", full.names = T)
aux_csv <- list.files(path = "data", pattern = ".csv", full.names = T)

# fiber_data <- readxl::read_xlsx(aux_sheets[grepl("LabGuruParticleCount", aux_sheets, ignore.case = T)], sheet = paste0(project_name, "_FIBER")) |> 
#   clean_names()

multiplier <- readxl::read_xlsx(aux_sheets[grepl("Multiplier", aux_sheets, ignore.case = T)], sheet = "Fragments_Fibers") |>
  clean_names()

filter_data <- read_csv(aux_csv[grepl("particle_details_all", aux_csv, ignore.case = T)]) #Not currently the correct one. 
multiplier2 <- readxl::read_xlsx(aux_sheets[grepl("Multiplier", aux_sheets, ignore.case = T)], sheet = "Filter") |>
  clean_names() |> 
  mutate(sample_id = str_replace_all(sample_id, "O", "0"))

# Set MIPPR sample ID to actual project ID
# Data
labguru_df <- readxl::read_xlsx(aux_sheets[grepl("LabGuruSampleWorksheet", aux_sheets, ignore.case = T)], sheet = 1) |> 
  clean_names()

#Ideally not necessary
# sampleid <- readxl::read_xlsx(aux_sheets[grepl("LabGuruSampleWorksheet", aux_sheets, ignore.case = T)], sheet = 2) |>
#   clean_names() |> 
#   filter(!is.na(mippr_sample_id)) |> 
#   select(sfei_sample_id, mippr_sample_id)

# Read in data
full_particle <- read.csv(list.files("data", "full_particle_results.csv", full.names = T))

particle_count <- read_sheet("https://docs.google.com/spreadsheets/d/1o3uoS5JW6JDxa4UGNht3v5fVNtlY9z8s4T1z8c1kmwQ/edit?usp=sharing",
                             sheet = project_name) %>%
  mutate(across(where(is.list), ~ purrr::map_chr(., toString))) # Converts lists to strings
  

##Hardcoding fixes optional before analysis
source("R/special_cleanup.R")
#Merge datasets
source("R/data_merge.R")
#Extract control data and report metadata
source("R/sample_analysis_plan.R")

rmarkdown::render(input = "R/MicroplasticsReport.Rmd", output_file = 
                    paste0(project_name,"_Report"),
                  output_format = "word_document")

#Upload all R code at end for reproducibility.
folder_find <- shared_drive_find("Project") |>
  drive_ls("Customer Projects") |>
  drive_ls() |>
  filter(name %in% project_name)

new_folder <- drive_mkdir("RCode",
                          path = as_id(folder_find$id),
                          overwrite = TRUE
)

data_upload <- c(list.files(path = "R/", 
                          pattern = ".R", 
                          full.names = T), 
                 "data/full_particle_data.csv")

for(file in data_upload){
  drive_upload(media = file.path(file),
               path = as_id(as.character(new_folder$id))
  )
}

