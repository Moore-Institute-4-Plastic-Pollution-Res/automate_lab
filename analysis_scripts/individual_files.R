library(OpenSpecy)
library(ggrepel)
library(data.table)
library(mmand)
library(magick)
library(cluster)
library(googledrive)

# spa files - individual files to run

# Load library for analysis
lib <- load_lib("derivative") |> 
  filter_spec(lib$metadata$spectrum_type == "ftir")

# Find files to download
data_search <- shared_drive_find("Project") |> 
  drive_ls("Customer Projects") |> 
  drive_ls() |> 
  filter(name %in% map(strsplit(project_name, split = "_"),1)) |> 
  drive_ls() |> 
  filter(name %in% project_name) |> 
  drive_ls() |> 
  filter(name == "Particle_500um_Data") |> 
  drive_ls() |> 
  filter(name == "iS50 Spectra") |> 
  drive_ls() |> 
  filter(name == "Spectra") |> 
  drive_ls() |> 
  filter(grepl("\\.(spa)|(csv)|(jdx)|(dx)|([0-9])$", name, ignore.case = T))
# Found the correct folder path - now create a function that conducts the analysis for each folder that starts with SFEI

#Set a folder to save data locally
local_store <- file.path("data",project_name)
local_store_raw <- file.path("data", project_name, "Raw", fsep = "/")
local_store_results <- file.path("data", project_name, "Results", fsep = "/")

# Create folders ----

if(!dir.exists(local_store)){
  dir.create(local_store, recursive = TRUE)
}
if(!dir.exists(local_store_raw)){
  dir.create(local_store_raw, recursive = TRUE)
}
if(!dir.exists(local_store_results)){
  dir.create(local_store_results, recursive = TRUE)
}

#Remove files already downloaded
data_search_to_download <- data_search |>
  filter(!name %in% list.files(local_store_raw))

# Download files
for (file in 1:nrow(data_search_to_download)) {
  # files in each specific folder
  file_named <- data_search_to_download$name[file]
  
  dest_path <- file.path(local_store_raw,file_named, fsep = "/")
  tryCatch({
    drive_download(as_id(data_search_to_download[file, ]$id), path = dest_path,overwrite = TRUE)
  }, error = function(e) {
    message("Error downloading file: ", file_named, "\n", e)
  }) #Add a check at the end that 
}

# for each folder conduct the analysis
individual_spec <- list.files(local_store_raw, full.names = T, recursive = T)

# Delete files that aren't in drive
# Loop through all files in the folder
for (file in individual_spec) {
  # Extract the file name from the full path
  file_name <- basename(file)
  
  # Check if the file is not in the list of names to keep
  if (!(file_name %in% data_search$name)) {
    # Delete the file
    file.remove(file)
    # Print a message
    print(paste0("Unmatched files deleted.\n", file))
  }
}

individual_spec <- list.files(local_store_raw, full.names = T, recursive = T)

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
  file = file.path(local_store_results,
    "full_particle_results.csv"
  ),
  row.names = F
)

just_plastics <- top_matches |> 
  dplyr::filter(material_class != "mineral",
                !grepl("organic", material_class))

write.csv(
  just_plastics,
  file = file.path(local_store_results,
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

drive_upload(media = file.path(local_store_results,
  "plastic_df.csv"
), path = as_id(as.character(new_folder$id)))


drive_upload(media = file.path(local_store_results,
  "full_particle_results.csv"
), path = as_id(as.character(new_folder$id)))