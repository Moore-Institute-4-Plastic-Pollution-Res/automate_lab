# spa files - individual files to run for SFEI
# Ask User the name of the folder they want to access
project_name <- readline(prompt = "Enter the Name of the Google Drive Folder for this Project: ")

# Load library for analysis
lib <- read_any("Library/derivative.rds")

lib <- filter_spec(lib, lib$metadata$spectrum_type == "ftir")



# Find files to download
data_search <- shared_drive_find("Project") |> 
  drive_ls("Customer Projects") |> 
  drive_ls() |> 
  filter(name %in% project_name) |> 
  drive_ls() |> 
  filter(name == "SFEI_01") |> 
  drive_ls() |> 
  filter(name == "Particle_500um_Data") |> 
  drive_ls() |> 
  filter(name == "iS50 Spectra") |> 
  drive_ls() |> 
  filter(name == "SFEI_01_Spectra_Particles") |> 
  drive_ls() |> 
  filter(grepl("SFEI", name))
  
# Found the correct folder path - now create a function that conducts the analysis for each folder that starts with SFEI

# Download files
# # Download data
for (folder in 1:nrow(data_search)) {
  # folders of interest 
  folder_name <- data_search[folder,]
  
  folders_sfei <- folder_name |> 
    drive_ls() |> 
    filter(grepl("spa", name)) 
  
  folder_placement <- file.path(project_name, "Raw", folder_name$name, fsep = "/")
  
  # Create folder paths
  # Create folders ----
  if(!dir.exists(file.path(project_name))){
    dir.create(file.path(project_name))
  }
  if(!dir.exists(file.path(folder_placement))){
    dir.create(file.path(folder_placement), recursive = T)
  }
  if(!dir.exists(file.path(project_name,"Results"))){
    dir.create(file.path(project_name,"Results"))
  }
  
    for (file in 1:nrow(folders_sfei)) {
      # files in each specific folder
      file_name <- folders_sfei$name[file]
      
      dest_path <- file.path(folder_placement,file_name, fsep = "/")
      tryCatch({
        drive_download(as_id(folders_sfei[file, ]$id), path = dest_path,overwrite = TRUE)
      }, error = function(e) {
        message("Error downloading file: ", file_name, "\n", e)
      })

    }
  
    # for each folder conduct the analysis
    individual_spec <- list.files(folder_placement, full.names = T)
    
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
    
    write.csv(
      top_matches,
      file = file.path(
        "SFEI",
        "Results",
        paste0(folder_name$name,"_individual_particle_matches.csv")
      ),
      row.names = F
    )
}


# Results join all the tables together and filter for plastic 
df <- list.files(path = "SFEI/Results/", full.names = T) |> 
  lapply(read_csv) |> 
  bind_rows()

df <- df |> 
  #rename(file_name = file_name.y) |> 
  mutate(sample = sub(".*_(S\\d{2})_.*", "\\1", file_name)) |> 
  select(sample,file_name,match_val,material_class,spectrum_identity,library_id)

write_csv(df, "full_particle_results.csv")

just_plastics <- df |> 
  filter(material_class != "mineral",
         !grepl("organic", material_class))

write.csv(just_plastics, "plastic_df.csv")

 # Upload results
# Create new folder
new_folder <- drive_mkdir("Single_Particle_Results",
                          path = as_id("1_Vsn_bHxkQvT5FEXPo6Aq3ug9yheZdIe"),
                          overwrite = TRUE
)



  drive_upload(media = "plastic_df.csv",
             path = as_id(as.character(new_folder$id)))


