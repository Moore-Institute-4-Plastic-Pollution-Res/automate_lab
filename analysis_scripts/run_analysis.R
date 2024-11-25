library(tidyverse)
library(googledrive)

# Ask User the name of the folder they want to access
project_name <- readline(prompt = "Enter the Name of the Google Drive Folder for this Project: ")

# Determine whether to conduct individual files (.SPA) or multi file (.hdr/.dat/.jpg) or both analysis
data_search <- data_search <- shared_drive_find("Project") |> 
  drive_ls("Customer Projects") |> 
  drive_ls() |> 
  filter(name %in% project_name) |> 
  drive_ls() 

if (data_search)

