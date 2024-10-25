# Create necessary folders
create_folders <- function(project_name) {
  if (!dir.exists(project_name)) {
    dir.create(project_name)
  }
  for (folder in folders_of_interest) {
    if (!dir.exists(paste0(project_name, "/Raw_Data"))) {
      dir.create(paste0(project_name, "/Raw_Data"))
    }
  }
}


# Download data
data_download <- function(drive_name, project_name) {
for (folder in folders_of_interest) {
  
  data_search <- shared_drive_find("Project") |> 
    drive_ls("Customer Projects") |> 
    drive_ls(project_name) |> 
    drive_ls(folder) |> 
    drive_ls()
  
  print(data_search)
  print(nrow(data_search))
  
  
  
  # # Download data
  for (file in 1:nrow(data_search)){
    
    # File name
    file_name <- data_search$name[file]
    # Set download path
    dest_path <- file.path(paste0(project_name,"/Raw_Data"), file_name)
    
    print(dest_path)
    # Download the file
    tryCatch({
      drive_download(data_search[file, ], path = dest_path, overwrite = TRUE)
    }, error = function(e) {
      message("Error downloading file: ", file_name, "\n", e)
    })
  }
}

}

