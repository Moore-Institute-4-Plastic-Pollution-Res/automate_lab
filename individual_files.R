# Ask User the name of the folder they want to access
project_name <- readline(prompt = "Enter the Name of the Google Drive Folder for this Project: ")

# Find files to download
data_search <- shared_drive_find("Project") |> 
  drive_ls("Customer Projects") |> 
  drive_ls() |> 
  filter(name %in% project_name) |> 
  drive_ls("Single_Particle_Spectra") |> 
  drive_ls()

# Download files
# # Download data
for (file in 1:nrow(data_search)){
  
  # File name
  file_name <- data_search$name[file]
  # Set download path
  dest_path <- file.path(paste0("Compost/Single_Particle/Raw"), file_name)
  
  print(dest_path)
  # Download the file
  tryCatch({
    drive_download(data_search[file, ], path = dest_path, overwrite = TRUE)
  }, error = function(e) {
    message("Error downloading file: ", file_name, "\n", e)
  })
}

wd <- "Compost/Single_Particle/Raw/"
individual_spec <- list.files(wd, full.names = T)

lib <- read_any("Library/derivative.rds")

lib <- filter_spec(lib, lib$metadata$spectrum_type == "ftir")

indiv_files <- read_any(individual_spec) |> 
    c_spec("common") |> 
    process_spec(conform_spec_args = list(range = lib$wavenumber, 
                                          res = NULL, 
                                          type = "interp"),
                 flatten_range = T
)

top_matches <- match_spec(x = indiv_files, 
           library = lib, 
           top_n = 1, 
           add_library_metadata = "sample_name",
           add_object_metadata = "col_id") %>%
    select(file_name.y, match_val, material_class, spectrum_identity, library_id)

write.csv(top_matches, 
          file.path("Compost/Single_Particle/Spectral_Results", "individual_particle_matches.csv"), row.names = F)

# Upload results
# Create new folder
new_folder <- drive_mkdir("Single_Particle_Results",
                          path = as_id("1hLtJD4gzJabIpNlG9m0xGH1X5KzPq_wO"),
                          overwrite = TRUE
)



drive_upload(media = "Compost/Single_Particle/Spectral_Results/individual_particle_matches.csv",
             path = as_id(as.character(new_folder$id)))


