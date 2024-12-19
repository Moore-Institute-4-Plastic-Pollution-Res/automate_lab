
#--------------Merging Particle Count and Fragment Analysis Results--------------
# Particle Count ----
particle_count2 <- particle_count |>
  mutate(
    # Replace characters in SampleID first
    SampleID = str_replace_all(SampleID, "-", "_"),
    SampleID = str_replace_all(SampleID, "W", "0"),
    # Extract the last two characters from WellID
    well = ifelse(!grepl("^F", WellID), str_sub(WellID, start = -2),  WellID),
    # Replace 'W' in the well portion
    well = str_replace_all(well, "W", "0"),
    
    # Now create ParticleID with updated SampleID and well
    ParticleID = paste0(SampleID, "_0", well)
  ) |>
  filter(str_length(ParticleID) > 6) |>
  group_by(ParticleID) |>
  distinct(ParticleLength, ParticleWidth, .keep_all = TRUE)

# # weird particle count row SFEI_01_S15_040
# particle_count <- particle_count |>
#   mutate(ParticleID =
#            case_when(
#              ParticleID == "SFEI_01_S15_00 " ~ "SFEI_01_S15_040",
#              TRUE ~ ParticleID
#            ))

# Full Particle analyzed ----
full_particle2 <- full_particle |>
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


full_particle3 <- full_particle2 |>
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
merged_data <- left_join(particle_count2, full_particle3, by = "ParticleID")

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
missing <- anti_join(particle_count2, fragment_data, by = "ParticleID")

# Decided to keep the missing but add a note for missing ***
# Good compromise here but also might want to flag these for a second look.
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
  mutate(color = ParticleColor, shape = ParticleShape) |>
  mutate(
    # orange
    color =
      case_when(
        grepl("^(YLO|ORA|GOLD)$", color, ignore.case = T)  ~ "Orange",
        # multi color
        str_count(color, "_") > 1 |
          grepl("MULTI", color) ~ "Multicolor (2+ colors)",
        #white
        grepl("^(W|WHT|CLR|SLVR|GRY)$", color, ignore.case = T)  ~ "White",
        #Blue
        grepl("^(BLU)$", color) ~ "Blue",
        # Black
        grepl("^(BLK|BLA)$", color, ignore.case = T) ~ "Black",
        # Brown
        grepl("^(BRN|TAN)$", color) ~ "Brown",
        # Green
        grepl("^GRN$", color, ignore.case = T) ~ "Green",
        # Pink
        grepl("^PNK$", color, ignore.case = T) ~ "Pink",
        # Red
        grepl("^RED$", color, ignore.case = T) ~ "Red",
        TRUE ~ color
      ),
    # Purple
    #shape
    shape =
      case_when(
        grepl("^(FRA|FLM|FOM|THN|CONG|LG|EMPTY|FB|FML|COMG|XL)$", 
              shape, ignore.case = T)  ~ "Fragment",
        grepl("^SPH$", 
              shape, ignore.case = T) ~ "Sphere",
        grepl("^(FBR|FIBER)$", 
              shape, ignore.case = T) ~ "Fiber Bundle",
        TRUE ~ shape
      )
    #ParticleShape =
    #  case_when(is.na(ParticleShape) ~ shape, TRUE ~ ParticleShape),
    #ParticleColor =
    #  case_when(is.na(ParticleColor) ~ color, TRUE ~ ParticleColor)
  ) |> 
  mutate(ParticleShape = shape, 
         ParticleColor = color)

# Manual check colors and shape
#Should add a test to the main code for this to double check things. 
unique(fragment_data1$shape)
unique(fragment_data1$color)

write.csv(fragment_data1,
          "data/fragment_data_full_final.csv", 
          row.names = F)


# ----------------------------- Multiplier ------------------------------------
# Multiplier by sample
# By each sample estimate of total particles found
# total of each type with the multiplier
# Should we do this now or once at end?
# Remove matches less than 0.6
fragment_data1 <- fragment_data1 |>
  # match_val > 0.6
  filter(match_val > 0.6)

unique(fragment_data1$ParticleColor)
unique(fragment_data1$ParticleShape)

fragment_total <- fragment_data1 |>
  group_by(SampleID, ParticleShape) |>
  summarize(count = n()) |>
  pivot_wider(names_from = ParticleShape, values_from = count) |>
  clean_names() |>
  left_join(multiplier) |>
  #select(-proportions_of_sample) |>
  mutate(across(.cols = where(is.numeric) & !all_of(c("subsample_ratio","multiplier")), ~ (.x /
                  subsample_ratio)/multiplier)) |>
  select(-multiplier) |>
  mutate(total_count = rowSums(across(everything()), na.rm = TRUE))


# count of plastics
plastics_total <- fragment_data1 |>
  filter(!(
    material_class %in% c(
      "mineral",
      "organic matter"
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
      "organic matter"
    )
  ), !is.na(ParticleShape)) |>
  group_by(SampleID) |>
  summarize(count_plastic = n()) |>
  clean_names()

fragment_particle_count <- left_join(plastics_total, plastics_breakdown) |>
  left_join(multiplier) |>
  select(-proportions_of_sample) |>
  mutate(across(.cols = where(is.numeric) &
                  !all_of(c("subsample_ratio","multiplier")), ~ (.x /
                  subsample_ratio)/multiplier)
         ) |>
  select(-c(multiplier, subsample_ratio)) |>
  mutate(count_plastic = floor(count_plastic)) |> 
  rename(SampleID = sample_id, "Particle Count" = count_plastic) 

# Particle count by polymer
fragment_breakdown <- fragment_data1 |>
  filter(!(
    material_class %in% c(
      "mineral",
      "organic matter"
    )
  ), !is.na(ParticleShape)) |>
  group_by(SampleID, material_class) |>
  summarize(count = n()) |>
  left_join(multiplier, by = c("SampleID" = "sample_id")) |>
  mutate(across(
    .cols = where(is.numeric) &
      !all_of(c("subsample_ratio", "multiplier")),
    ~ (.x/subsample_ratio) / multiplier
  )
  ) |>
  select(-c(multiplier, proportions_of_sample)) |>
  group_by(material_class) |>
  summarize(count = sum(count)) |> 
  mutate(count = floor(count))