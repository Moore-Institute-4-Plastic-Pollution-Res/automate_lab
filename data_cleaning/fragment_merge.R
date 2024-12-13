library(tidyverse)
library(magrittr)
library(janitor)
library(readxl)

# Read in data
particle_count <- read.csv(
  "data_cleaning/fragment_data/SFEI_01_ParticleCount - Fragment Data.csv",
  na.strings = c("", "NA")
)
full_particle <- read.csv("data_cleaning/fragment_data/SFEI_full_particle_results.csv")

#--------------Merging Particle Count and Fragment Analysis Results--------------
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
    ParticleID = paste0(SampleID, "_0", well)
  ) |>
  filter(str_length(ParticleID) > 6) |>
  group_by(ParticleID) |>
  distinct(ParticleLength, ParticleWidth, .keep_all = TRUE)

# weird particle count row SFEI_01_S15_040
particle_count <- particle_count |>
  mutate(ParticleID =
           case_when(
             ParticleID == "SFEI_01_S15_00 " ~ "SFEI_01_S15_040",
             TRUE ~ ParticleID
           ))

# Full Particle analyzed ----
full_particle <- full_particle |>
  mutate(
    ParticleID = str_extract(file_name, "^((?:[^_]+_){3}[^_]+)"),
    shape = str_match(file_name, "^(?:[^_]+_){4}([^_]+)")[, 2],
    color = str_match(file_name, "^(?:[^_]+_){5}([^.]+)")[, 2],
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
           case_when(
             grepl("\\.spa", ParticleID) ~ str_replace_all(ParticleID, "\\.spa", ""),
             TRUE ~ ParticleID
           ))


full_particle1 <- full_particle |>
  mutate(
    color =
      case_when(
        grepl("BlueBead", color) ~ "blue",
        grepl("\\+BLU", shape) |
          grepl("Bead", shape) | grepl("Blu", shape) ~ "BLU",
        TRUE ~ color
      ),
    shape = case_when(
      grepl("BlueBead", shape) | grepl("Blu", shape) ~ "Bead",
      grepl("EMPTY", shape) |
        grepl("Empty", shape) |
        grepl("Background", shape) ~ NA_character_,
      grepl("PET", shape) ~ "PET",
      grepl("\\+BLU", shape) ~ "FRA",
      TRUE ~ shape
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
  filter(grepl("~", file_name) | grepl("REDO", file_name))

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
color <- c(
  "Black",
  "Blue",
  "Brown",
  "Green",
  "Multicolor (2+ colors)",
  "Pink",
  "Purple",
  "Red",
  "White",
  "Orange"
)
# particle shape
particle_shape <- c("Fragment", "Sphere", "Rubbery Fragment")


fragment_data1 <- fragment_data |>
  mutate(
    # orange
    color =
      case_when(
        grepl("Ylo", color) | grepl("YLO", color) |
          grepl("ORA", color) | grepl("GOLD", color) ~ "Orange",
        # multi color
        str_count(color, "_") > 1 |
          grepl("MULTI", color) ~ "Multicolor (2+ colors)",
        #white
        grepl("W", color) | grepl("CLR", color) |
          grepl("SLVR", color) | grepl("GRY", color) |
          grepl("Gry", color) ~ "White",
        #Blue
        grepl("BLU", color) ~ "Blue",
        # Black
        grepl("BLK", color) |
          grepl("Blk", color) | grepl("BLA", color) ~ "Black",
        # Brown
        grepl("BRN", color) | grepl("TAN", color) ~ "Brown",
        # Green
        grepl("GRN", color) | grepl("Grn", color) ~ "Green",
        # Pink
        grepl("PNK", color) | grepl("Pnk", color) ~ "Pink",
        # Red
        grepl("RED", color) ~ "Red",
        TRUE ~ color
      ),
    # Purple
    #shape
    shape =
      case_when(
        grepl("FRA", shape) | grepl("Fra", shape) |
          grepl("FLM", shape) | grepl("FOM", shape) |
          grepl("THN", shape) | grepl("CONG", shape) |
          grepl("LG", shape) | grepl("Empty", shape) |
          grepl("FB", shape) | grepl("FML", shape) |
          grepl("COMG", shape) |
          grepl("XL", shape) ~ "Fragment",
        grepl("SPH", shape) ~ "Sphere",
        TRUE ~ shape
      ),
    ParticleShape =
      case_when(is.na(ParticleShape) ~ shape, TRUE ~ ParticleShape),
    ParticleColor =
      case_when(is.na(ParticleColor) ~ color, TRUE ~ ParticleColor)
  )

unique(fragment_data1$shape)
unique(fragment_data1$color)

# Manual check colors and shape
# FLM3
fragment_data1 <- fragment_data1 |>
  mutate(
    ParticleColor =
      case_when(
        grepl("FLM3", ParticleColor) ~ "White",
        grepl("SFEI_01_S15_066", ParticleID) ~ "Orange",
        grepl("corr", ParticleColor) ~ NA,
        TRUE ~ ParticleColor
      ),
    ParticleShape =
      case_when(
        grepl("SFEI_01_S15_066", ParticleID) ~ "Fragment",
        TRUE ~ ParticleShape
      )
  ) |>
  select(-dupes)  |>
  # match_val > 0.6
  filter(match_val > 0.6)

unique(fragment_data1$ParticleColor)
unique(fragment_data1$ParticleShape)

write.csv(fragment_data1,
          "data_cleaning/fragment_data/SEI_fragment_data_full_final.csv")


# ----------------------------- Multiplier ------------------------------------
multiplier <- read_xlsx("data_cleaning/SFEI_01_Multiplier.xlsx", sheet = "Fragments_Fibers") |>
  clean_names()


# Multiplier by sample
# By each sample estimate of total particles found
# total of each type with the multiplier
fragment_total <- fragment_data1 |>
  group_by(SampleID, ParticleShape) |>
  summarize(count = n()) |>
  pivot_wider(names_from = ParticleShape, values_from = count) |>
  clean_names() |>
  left_join(multiplier) |>
  select(-proportions_of_sample) |>
  mutate(across(.cols = where(is.numeric) & !all_of("multiplier"), ~ .x /
                  multiplier)) |>
  select(-multiplier) |>
  mutate(total_count = rowSums(across(everything()), na.rm = TRUE))


# count of plastics
plastics_total <- fragment_data1 |>
  filter(!(
    material_class %in% c(
      "mineral",
      "organic matter",
      "cellulose derivatives (ether cellulose)"
    )
  ), !is.na(ParticleShape)) |>
  group_by(SampleID) |>
  summarize(count_plastic = n()) |>
  clean_names()


# break down of plastics
plastics_breakdown <- fragment_data1 |>
  filter(!(
    material_class %in% c(
      "mineral",
      "organic matter",
      "cellulose derivatives (ether cellulose)"
    )
  ), !is.na(ParticleShape)) |>
  mutate(
    match_rate =
      case_when(
        match_val > 0.60 ~ "above 60%",
        match_val < 0.60 & match_val > 0.30 ~ "above 30%",
        TRUE ~ "0"
      )
  ) |>
  filter(match_rate != "0") |>
  group_by(SampleID, match_rate) |>
  summarize(count_plastic = n()) |>
  pivot_wider(names_from = match_rate, values_from = count_plastic) |>
  clean_names()

plastics_df <- left_join(plastics_total, plastics_breakdown) |>
  left_join(multiplier) |>
  select(-proportions_of_sample) |>
  mutate(across(.cols = where(is.numeric) & !all_of("multiplier"), ~ .x /
                  multiplier)) |>
  select(-multiplier) |>
  rename(SampleID = sample_id, "Particle Count" = count_plastic)

# Particle count by polymer
plastics_breakdown <- fragment_data1 |>
  filter(!(
    material_class %in% c(
      "mineral",
      "organic matter",
      "cellulose derivatives (ether cellulose)"
    )
  ), !is.na(ParticleShape)) |>
  group_by(SampleID, material_class) |>
  summarize(count = n()) |>
  left_join(multiplier, by = c("SampleID" = "sample_id")) |>
  mutate(across(.cols = where(is.numeric) & !all_of("multiplier"), ~ .x /
                  multiplier)) |>
  select(-c(multiplier, proportions_of_sample)) |>
  group_by(material_class) |>
  summarize(count = sum(count))



write.csv(plastics_df,
          "data_cleaning/final/SFEI_fragment_data_final.csv")
write.csv(plastics_breakdown,
          "data_cleaning/final/SFEI_fragment_polymer_count.csv")
