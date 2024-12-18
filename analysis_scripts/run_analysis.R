library(tidyverse)
library(googledrive)
# install OpenSpecy
#devtools::install_github("wincowgerDEV/OpenSpecy-package")

library(OpenSpecy)

# Read in all functions
lapply(list.files("R", full.names = TRUE), source)

# Load library for analysis
lib <- OpenSpecy::load_lib("derivative") |> 
  filter_spec(lib$metadata$spectrum_type == "ftir")

# Ask User the name of the folder they want to access
project_name <- readline(prompt = "Enter the Name of the Project: ")

# Find folder names present in project folder - heavily dependent if it contains _ hyphen
folder_find <- shared_drive_find("Project") |>
  drive_ls("Customer Projects") |>
  drive_ls() |>
  filter(name %in% map(strsplit(project_name, split = "_"),1)) |>
  drive_ls() |>
  filter(name %in% project_name) |>
  drive_ls()

# Single project with the customer name
# folder_find <- shared_drive_find("Project") |>
#   drive_ls("Customer Projects") |>
#   drive_ls() |>
#   filter(name %in% project_name) |>
#   drive_ls("JHENG_270") |>
#   drive_ls("JHENG270") |>
#   drive_ls()

# Determine whether to conduct individual files (.SPA) or multi file (.hdr/.dat/.jpg) or both analysis
if (any(grepl("Particle_500um", folder_find$name))) {
  source("analysis_scripts/individual_files.R")
} else if (any(grepl("Export_Files", folder_find$name))) {
  source("analysis_scripts/multi_files.R")
} else if (any(grepl("Particle_500um", folder_find$name)) | any(grepl("Export_Files", folder_find$name))) {
  source("analysis_scripts/individual_files.R")
  source("analysis_scripts/multi_files.R")
}


# Generate report
# source file
source("data_cleaning/data_merge.R")
source("analysis_scripts/sample_analysis_plan.R")

rmarkdown::render("MicroplasticsReport.Rmd", output_format = "all")
