library(tidyverse)

# download particle count locally
particle_count_data <- folder_find |> 
  filter(grepl("ParticleCount", name)) |> 
  # asummes particle count on the firt sheet of the spreadsheet
  drive_download(path = file.path(local_store, "particle_count.csv"), overwrite = TRUE, type = "csv")

# fiber
fiber_count_data <- folder_find |> 
  filter(grepl("Fiber_data", name)) |> 
  # asummes particle count on the firt sheet of the spreadsheet
  drive_download(path = file.path(local_store, "fiber_data"), overwrite = TRUE, type = "xlsx")

# download Multiplier
multiplier_data <- folder_find |> 
  filter(grepl("Multiplier", name)) |> 
  # asummes particle count on the firt sheet of the spreadsheet
  drive_download(path = file.path(local_store, "multiplier"), overwrite = TRUE, type = "xlsx")

# lab guru sheet
lab_guru <- folder_find |> 
  filter(grepl("LabGuru", name)) |> 
  drive_download(path = file.path(local_store, "lab_guru"), overwrite = TRUE, 
                 type = "xlsx")
