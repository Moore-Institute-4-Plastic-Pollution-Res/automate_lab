library(tidyverse)
library(magrittr)

# Read in data
particle_count <- read.csv("SFEI/fragment_data/SFEI_01_ParticleCount - Fragment Data.csv") 
full_particle <- read.csv("SFEI/fragment_data/full_particle_results.csv")

particle_count <- particle_count |> 
  mutate(
    # Replace characters in SampleID first
    SampleID = str_replace_all(SampleID, "-", "_"),
    SampleID = str_replace_all(SampleID, "W", "0"),
    
    # Extract the last two characters from WellID
    well = str_sub(WellID, start = -2),
    # Replace 'W' in the well portion
    well = str_replace_all(well, "W", "0"),
    
    # Now create ParticleID with updated SampleID and well
    ParticleID = paste0(SampleID,"_0", well)
  ) |> 
  filter(str_length(ParticleID) > 5) |> 
  group_by(ParticleID) |> 
  distinct(ParticleLength, ParticleWidth, .keep_all = TRUE)

# weird particle count row SFEI_01_S15_040
 particle_count <- particle_count |> 
  mutate(ParticleID = 
           case_when(ParticleID == "SFEI_01_S15_00 " ~ "SFEI_01_S15_040",
                     TRUE ~ ParticleID
                     )
           )


# particle count
full_particle <- full_particle |> 
  mutate(
    ParticleID = str_extract(file_name, "^((?:[^_]+_){3}[^_]+)"),
    shape = str_match(file_name, "^(?:[^_]+_){4}([^_]+)")[,2],
    color = str_match(file_name, "^(?:[^_]+_){5}([^.]+)")[,2],
    # If a tilde is present, extract the number after it, otherwise assign "10"
    dupes = case_when(
      str_detect(file_name, "~") ~ str_extract(file_name, "(?<=~)\\d+"),
      TRUE ~ "10"
    ),
    dupes = as.numeric(dupes),
    # Create a base name by removing the ~number.spa pattern
    base_name = str_remove(file_name, "~\\d+\\.spa$")
  ) |>
  group_by(base_name) |> 
  slice_max(dupes, with_ties = FALSE) |> 
  ungroup() |> 
  select(-base_name) |> 
# check if ParticleID has .spa?
  mutate(ParticleID = 
           case_when(grepl("\\.spa", ParticleID) ~ str_replace_all(ParticleID, "\\.spa",""),
                     TRUE ~ ParticleID
                       )
           )


full_particle1 <- full_particle |>
  mutate(
    color =
      case_when(grepl("BlueBead", color) ~ "blue", TRUE ~ color),
    shape = case_when(grepl("BlueBead", shape) ~ "Bead", TRUE ~ shape),
    shape = case_when(
      grepl("EMPTY", shape) |
        grepl("Empty", shape) |
        grepl("Background", shape) ~ NA_character_,
      TRUE ~ shape
    ),
    shape =
      case_when(grepl("Pet", shape) |
                  grepl("PET", shape) ~ "PET", TRUE ~ shape),
    color =
      case_when(grepl("\\+BLU", shape) ~ "BLU", TRUE ~ color),
    shape = 
      case_when(grepl("\\+BLU", shape) ~ "FRA", TRUE ~ shape),
    # replace NA for blue in blue bead
    color = 
      case_when(
        grepl("Bead", shape) ~ "BLU", TRUE ~ color
      ),
    color = 
      case_when(
        grepl("Blu", shape) ~ "BLU", TRUE ~ color
      ),
    shape = 
      case_when(
        grepl("Blu", shape) ~ "Bead", TRUE ~ shape
      )
  ) 

# merge tables and drop the ones not present in both
merged_data <- left_join(particle_count, full_particle, by = "ParticleID")

merged_data <- merged_data |> 
  filter(!is.na(file_name)) |> 
  group_by(ParticleID)

# 'df' is your dataframe and 'col' is the column of interest
dupes <- merged_data[duplicated(merged_data$ParticleID), ]

dupes_name <- dupes$ParticleID

# remove these from dataframe
df_1 <- merged_data |> 
  filter(!(ParticleID %in% dupes_name))

# duplicate rows
dupes_df <- merged_data |> 
  filter(ParticleID %in% dupes_name) |> 
  group_by(ParticleID) |> 
  filter(grepl("~",file_name) | grepl("REDO", file_name))


fragment_data <- rbind(df_1, dupes_df)


missing <- anti_join(particle_count, fragment_data)
