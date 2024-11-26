library(tidyverse)
library(googledrive)

# Read in all functions
lapply(list.files("R", full.names = TRUE), source)

# Ask User the name of the folder they want to access
project_name <- readline(prompt = "Enter the Name of the Project: ")

# Find folder names present in project folder 
folder_find <- shared_drive_find("Project") |> 
  drive_ls("Customer Projects") |> 
  drive_ls() |> 
  filter(name %in% map(strsplit(project_name, split = "_"),1)) |> 
  drive_ls() |> 
  filter(name %in% project_name) |> 
  drive_ls()

# Determine whether to conduct individual files (.SPA) or multi file (.hdr/.dat/.jpg) or both analysis
if (grepl(folder_find$name, "Particle_500um") ){
  source("analysis_scripts/individual_files.R")
} else if (grepl(folder_find$name, "Export_Files")){
  source("analysis_scripts/multi_files.R")
} else if (grepl(folder_find$name, "Particle_500um") | grepl(folder_find$name, "Export_Files")){
  source("analysis_scripts/individual_files.R")
  source("analysis_scripts/multi_files.R")
}


