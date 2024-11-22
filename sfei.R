library(OpenSpecy)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(data.table)
library(mmand)
library(magick)
library(cluster)
library(googledrive)
library(purrr)
library(readr)


# spa files - individual files to run for SFEI
# Ask User the name of the folder they want to access
project_name <- readline(prompt = "Enter the Name of the Google Drive Folder for this Project: ")

# Load library for analysis
lib <- load_lib("derivative")

lib <- filter_spec(lib, lib$metadata$spectrum_type == "ftir")

# Find files to download
data_search <- shared_drive_find("Project") |> 
  drive_ls("Customer Projects") |> 
  drive_ls() |> 
  filter(name %in% project_name) |> 
  drive_ls() |> 
  filter(name == "Particle_500um_Data") |> 
  drive_ls() |> 
  filter(name == "iS50 Spectra") |> 
  drive_ls() |> 
  filter(name == "Spectra") |> 
  drive_ls()
  
# Found the correct folder path - now create a function that conducts the analysis for each folder that starts with SFEI

#Set a folder to save data locally
folder_placement <- file.path(project_name, "Raw", fsep = "/")

# Create folders ----
if(!dir.exists(file.path(project_name))){
  dir.create(file.path(project_name))
  dir.create(file.path(project_name, "Raw", fsep = "/"))
  dir.create(file.path(project_name, "Results", fsep = "/"))
}

# Download files
    for (file in 1:nrow(data_search)) {
      # files in each specific folder
      file_named <- data_search$name[file]
      
      dest_path <- file.path(folder_placement,file_named, fsep = "/")
      tryCatch({
        drive_download(as_id(data_search[file, ]$id), path = dest_path,overwrite = TRUE)
      }, error = function(e) {
        message("Error downloading file: ", file_name, "\n", e)
      }) #Add a check at the end that 
    }

    # for each folder conduct the analysis"
    individual_spec <- list.files(file.path(project_name, "Raw", fsep = "/"), full.names = T, recursive = T)
    
    indiv_files <- read_any(individual_spec) |>
      c_spec("common") |>
      # 800 - 3800
      process_spec(
        restrict_range = TRUE,
        restrict_range_args = list(min = 800, max = 3800),
        conform_spec_args =
          list(
            range = lib$wavenumber,
            res = NULL,
            type = "interp"
          ),
        flatten_range = T
      )
    
    top_matches <- match_spec(x = indiv_files, 
                              library = lib, 
                              top_n = 1,
                              add_library_metadata = "sample_name",
                              add_object_metadata = "col_id") %>%
      rename(file_name = file_name.y) |> 
      mutate(sample = sub(".*_(S\\d{2})_.*", "\\1", file_name)) |> 
      select(sample,file_name,match_val,material_class,spectrum_identity,library_id)
    
    
    write.csv(
      top_matches,
      file = file.path(
        project_name,
        "Results",
        "full_particle_results.csv"
      ),
      row.names = F
    )

just_plastics <- top_matches |> 
  dplyr::filter(material_class != "mineral",
         !grepl("organic", material_class))

write.csv(
  just_plastics,
  file = file.path(
    project_name,
    "Results",
    "plastic_df.csv"
  ),
  row.names = F
)

# Upload results
# Create new folder

folder_location <- shared_drive_find("Project") |> 
  drive_ls("Customer Projects") |> 
  drive_ls() |> 
  filter(name %in% project_name) |> 
  drive_ls() |> 
  filter(name == "Particle_500um_Data") |> 
  drive_ls() |> 
  filter(name == "iS50 Spectra") |> 
  drive_ls() |> 
  filter(name == "Spectra")

new_folder <- drive_mkdir("Single_Particle_Results",
                          path = as_id(folder_location$id),
                          overwrite = TRUE
)

drive_upload(media = file.path(
  project_name,
  "Results",
  "plastic_df.csv"
), path = as_id(as.character(new_folder$id)))


drive_upload(media = file.path(
  project_name,
  "Results",
  "full_particle_results.csv"
), path = as_id(as.character(new_folder$id)))
