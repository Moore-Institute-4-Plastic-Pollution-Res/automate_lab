library(tidyverse)
library(magrittr)
library(janitor)
library(readxl)

# Read in data
particle_count <- read.csv("SFEI/fragment_data/SFEI_01_ParticleCount - Fragment Data.csv",
                           na.strings = c("", "NA")) 
full_particle <- read.csv("SFEI/fragment_data/full_particle_results.csv")

# Particle Count ----
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
  filter(str_length(ParticleID) > 6) |> 
  group_by(ParticleID) |> 
  distinct(ParticleLength, ParticleWidth, .keep_all = TRUE)

# weird particle count row SFEI_01_S15_040
 particle_count <- particle_count |> 
  mutate(ParticleID = 
           case_when(ParticleID == "SFEI_01_S15_00 " ~ "SFEI_01_S15_040",
                     TRUE ~ ParticleID
                     )
           )

# Full Particle analyzed ----
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

# Merge tables ----
merged_data <- left_join(particle_count, full_particle, by = "ParticleID")

# Filter rows with NA value
merged_data <- merged_data |> 
  filter(!is.na(file_name)) 

# 'df' is your dataframe and 'col' is the column of interest
dupes <- merged_data[duplicated(merged_data$ParticleID), ]

dupes_name <- dupes$ParticleID

# Separate non duplicated and duplicated data
df_1 <- merged_data |> 
  filter(!(ParticleID %in% dupes_name))

# Duplicated data ----
# Selection of duplicated sample to keep 
dupes_df <- merged_data |> 
  filter(ParticleID %in% dupes_name) |> 
  group_by(ParticleID) |> 
  filter(grepl("~",file_name) | grepl("REDO", file_name))

# Recombine non duplicated and modified duplicated data 
fragment_data <- rbind(df_1, dupes_df)


# Identify missing samples form particle_count
missing <- anti_join(particle_count, fragment_data)

# Decided to keep the missing but add a note for missing ***
missing <- missing |> 
  mutate(Notes = "did not analyze")

fragment_data <- rbind(fragment_data, missing)

# Fill in particle shape and color ----
# color 
color <- c("Black", "Blue", "Brown", "Green", "Multicolor (2+ colors)", 
           "Pink", "Purple", "Red", "White", "Orange")
# particle shape
particle_shape <- c("Fragment", "Sphere", "Rubbery Fragment")


fragment_data1 <- fragment_data |> 
  mutate(
    # orange
    color =
      case_when(
        grepl("Ylo", color) | grepl("YLO", color) | 
          grepl("ORA", color) | grepl("GOLD", color) ~ "Orange",
        TRUE ~ color
      ),
    # multi color
    color =
      case_when(
        str_count(color, "_") > 1 | grepl("MULTI", color) ~ "Multicolor (2+ colors)", TRUE ~ color
      ),
    #white
    color =
      case_when(
        grepl("W", color) | grepl("CLR", color) | 
          grepl("SLVR", color) | grepl("GRY", color) | 
          grepl("Gry", color) ~ "White", TRUE ~ color
      ),
    #Blue
    color =
      case_when(
        grepl("BLU", color) ~ "Blue", TRUE ~ color
      ),
    # Black
    color = 
      case_when(
        grepl("BLK", color) | grepl("Blk", color) | grepl("BLA", color) ~ "Black", TRUE ~ color
      ),
    # Brown
    color = case_when(
      grepl("BRN", color) | grepl("TAN", color) ~ "Brown", TRUE ~ color
    ),
    # Green
    color = case_when(
      grepl("GRN", color) | grepl("Grn", color) ~ "Green", TRUE ~ color
    ),
    # Pink
    color = case_when(
      grepl("PNK", color) | grepl("Pnk", color) ~ "Pink", TRUE ~ color
    ),
    # Red
    color = case_when(
      grepl("RED", color) ~ "Red", TRUE ~ color
    ),
    # Purple
    
    #shape
    shape = 
      case_when(grepl("FRA", shape) | grepl("Fra", shape) |
                  grepl("FLM", shape) | grepl("FOM", shape) |
                  grepl("THN", shape) | grepl("CONG", shape) |
                  grepl("LG", shape) | grepl("Empty", shape) | 
                  grepl("FB", shape) | grepl("FML", shape) | 
                  grepl("COMG", shape) | grepl("XL", shape) ~ "Fragment", 
                TRUE ~ shape
                ),
    shape = case_when(grepl("SPH", shape) ~ "Sphere", TRUE ~ shape),
    ParticleShape = 
      case_when(is.na(ParticleShape) ~ shape, 
                TRUE ~ ParticleShape),
    ParticleColor = 
      case_when(is.na(ParticleColor) ~ color,
                TRUE ~ ParticleColor
                )
  )

unique(fragment_data1$shape)
unique(fragment_data1$color)

# Manual check colors and shape
# FLM3
fragment_data1 <- fragment_data1 |> 
  mutate( 
    ParticleColor = 
    case_when(grepl("FLM3", ParticleColor) ~ "White", TRUE ~ ParticleColor),
    ParticleColor = 
      case_when(grepl("corr", ParticleColor) ~ NA, TRUE ~ ParticleColor),
    ParticleShape = 
      case_when(grepl("SFEI_01_S15_066", ParticleID) ~ "Fragment", TRUE ~ ParticleShape),
    ParticleColor = 
      case_when(grepl("SFEI_01_S15_066", ParticleID) ~ "Orange", TRUE ~ ParticleColor)
  )

unique(fragment_data1$ParticleColor)
unique(fragment_data1$ParticleShape)

write.csv(fragment_data1, "SFEI/fragment_data/fragment_data_full_final.csv")


# ----------------------------- Multiplier ------------------------------------
multiplier <- read_xlsx("SFEI/SFEI_01_Multiplier.xlsx", sheet = "Fragments_Fibers") |> 
  clean_names()


# Multiplier by sample
# By each sample estimate of total particles found
fragment_total <- fragment_data1 |> 
  group_by(SampleID) |> 
  summarize(count = n()) |> 
  clean_names() |> 
  left_join(multiplier) |> 
  select(-proportions_of_sample) |> 
  mutate(total = count/multiplier)

# total fragment
