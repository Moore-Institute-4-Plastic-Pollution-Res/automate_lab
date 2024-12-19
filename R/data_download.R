library(tidyverse)

# download particle count locally
particle_count_data <- folder_find |> 
  filter(grepl("ParticleCount", name)) |> 
  # assumes particle count on the first sheet of the spreadsheet
  drive_download(path = file.path(local_store, "particle_count.csv"), overwrite = TRUE, type = "csv")

# fiber
fiber_count_data <- folder_find |> 
  filter(grepl("Fiber_data", name)) |> 
  # assumes particle count on the first sheet of the spreadsheet
  drive_download(path = file.path(local_store, "fiber_data"), overwrite = TRUE, type = "xlsx")

# download Multiplier
multiplier_data <- folder_find |> 
  filter(grepl("Multiplier", name)) |> 
  # asummes particle count on the firt sheet of the spreadsheet
  drive_download(path = file.path(local_store, "Multiplier"), overwrite = TRUE, type = "xlsx")

# lab guru sheet
lab_guru <- folder_find |> 
  filter(grepl("LabGuru", name)) |> 
  drive_download(path = file.path(local_store, "LabGuruSampleWorksheet"), overwrite = TRUE, 
                 type = "xlsx")

# LFB QA/QC
lfb_counts <- folder_find |> 
  filter(grepl("LFB_Counts", name)) |> 
  drive_download(path = file.path(local_store, "LFB_Counts"), overwrite = TRUE, 
                 type = "csv")

