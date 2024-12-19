# Merge plastic count data for fragment, fiber, and filter 
# Read in data

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
  ) #|>
  #filter(str_length(ParticleID) > 6) |>
  #group_by(ParticleID) |>
  #distinct(ParticleLength, ParticleWidth, .keep_all = TRUE)

# Full Particle analyzed ----
full_particle2 <- full_particle |>
  mutate(
    ParticleID = str_extract(file_name, "^((?:[^_]+_){3}[^_]+)"),
    #shape = str_match(file_name, "^(?:[^_]+_){4}([^_]+)")[, 2],
    #color = str_match(file_name, "^(?:[^_]+_){5}([^.]+)")[, 2],
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

# Merge tables ----
merged_data <- left_join(particle_count2, full_particle2, by = "ParticleID")

# Filter rows with NA value
merged_data <- merged_data |>
  filter(!is_empty_or_na(file_name))

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
fragment_data <- bind_rows(df_1, dupes_df)

# Identify missing samples form particle_count
# Decided to keep the missing but add a note for missing ***
# Good compromise here but also might want to flag these for a second look.

missing <- anti_join(particle_count2, fragment_data, by = "ParticleID") |>
  mutate(Notes = ifelse(!grepl("^F", WellID) & !grepl("MIPPR", SampleID), "did not analyze", "")) 

fragment_data <- bind_rows(fragment_data, missing)

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
    ParticleColor =
      case_when(
        grepl("^(YLO|ORA|GOLD|yellow)$", ParticleColor, ignore.case = T)  ~ "Orange",
        # multi ParticleColor
        str_count(ParticleColor, "_") > 1 |
          grepl("MULTI", ParticleColor, ignore.case = T) ~ "MultiParticleColor (2+ colors)",
        #white
        grepl("^(W|WHT|CLR|SLVR|GRY|clear)$", ParticleColor, ignore.case = T)  ~ "White",
        #Blue
        grepl("^(BLU|Dk Blue|Blue)$", ParticleColor, ignore.case = T) ~ "Blue",
        # Black
        grepl("^(BLK|BLA)$", ParticleColor, ignore.case = T) ~ "Black",
        # Brown
        grepl("^(BRN|TAN)$", ParticleColor, ignore.case = T) ~ "Brown",
        # Green
        grepl("^GRN$", ParticleColor, ignore.case = T) ~ "Green",
        # Pink
        grepl("^PNK$", ParticleColor, ignore.case = T) ~ "Pink",
        # Red
        grepl("^RED$", ParticleColor, ignore.case = T) ~ "Red",
        TRUE ~ ParticleColor
      ),
    # Purple
    #shape
    ParticleShape =
      case_when(
        grepl("^(FRA|FLM|FOM|THN|CONG|LG|EMPTY|FB|FML|COMG|XL)$", 
              ParticleShape, ignore.case = T)  ~ "Fragment",
        grepl("^SPH$", 
              ParticleShape, ignore.case = T) ~ "Sphere",
        grepl("^(FBR|FIBER)$", 
              ParticleShape, ignore.case = T) & !grepl("^F", WellID) ~ "Fiber Bundle",
        TRUE ~ ParticleShape
      )
  ) |> 
  mutate(ParticlePolymerType = tolower(ParticlePolymerType),
         ParticlePolymerType = 
           case_when(
             grepl("cellulose", ParticlePolymerType, ignore.case = T) ~ "cellulose derivative (ether cellulose)",
             grepl("organic", ParticlePolymerType, ignore.case = T) ~ "organic matter",
             grepl("(polycrylon)", ParticlePolymerType, ignore.case = T) ~ 
               "polyacrylonitrile",
             grepl("polypro", ParticlePolymerType, ignore.case = T) ~ "polypropylene",
             grepl("un", ParticlePolymerType, ignore.case = T) ~ "unknown",
             grepl("polyester|terephthalates", ParticlePolymerType, ignore.case = T) ~ "poly(esters/ethers/diglycidylethers/terephthalates)s",
             is_empty_or_na(ParticlePolymerType) ~ material_class,
             TRUE ~ ParticlePolymerType
           )
  ) |>
  mutate(match_val = ifelse(is_empty_or_na(match_val), Correlation, match_val)) |>
  mutate(match_val = ifelse(grepl("^F", WellID) & is_empty_or_na(match_val), 0.66, match_val)) |>
  ungroup() |>
  select(SampleID, WellID, ParticleShape, ParticleColor, 
         ParticlePolymerType, ParticleLength, ParticleWidth,
         SieveSizeRange, Notes, library_id, match_val) |>
  clean_names() |>
  rename(particle_id = well_id,
         max_cor_val = match_val, 
         max_cor_name = library_id, 
         max_length_um = particle_length, 
         min_length_um = particle_width, 
         material_class = particle_polymer_type)

# Manual check colors and shape
#Should add a test to the main code for this to double check things. 
print(unique(fragment_data1$particle_shape))
print(unique(fragment_data1$particle_color))

# Plastic count merge

all_data <- bind_rows(filter_data |> 
                        mutate(particle_id = as.character(particle_id)) |>
                        mutate(
                          sample_id = str_replace_all(sample_id, "_REDO",""),
                          sample_id = str_replace_all(sample_id, "O", "0"),
                          sample_name = str_replace(sample_id, "^([^_]+_[^_]+)_(S|W)(\\d+)_.*", "\\1_\\2\\3"),
                          sample_name = gsub("(_[0-9]of[0-9])|(_AQ[0-9])", "", sample_id)
                        ), 
                      fragment_data1 |> mutate(max_length_um= as.numeric(max_length_um),
                                              min_length_um= as.numeric(min_length_um)
                                              )) |> 
  mutate(count = 1) |>
  left_join(multiplier, by = "sample_id")|>
  left_join(multiplier2, by = c("sample_name" = "sample_id")) |>
  mutate(sample_name = ifelse(is_empty_or_na(sample_name), sample_id, sample_name)) |>
  mutate(count = ifelse(is_empty_or_na(area_um2), count / multiplier.x ,  count / multiplier.y))  

#Major client output
fwrite(all_data, "full_particle_data.csv")

#cleanup and analysis. 
confident_plastic <- all_data |>
  filter(
    max_cor_val >= 0.6,
    !(
      material_class %in% c(
        "other material",
        "mineral",
        "organic matter",
        "unknown"
      )
    )) 

confident_plastic_nocontrol <- confident_plastic |>
  filter(!str_detect(sample_id, "MIPPR")) 

#assuming that the large and small particle multipliers have different IDs here. 
sample_plastic <- confident_plastic_nocontrol |>
  group_by(sample_name) |> 
  summarise(count = round(sum(count), 0)) |>
  ungroup()    

#assuming that the large and small particle multipliers have different IDs here. 
material_plastic <- confident_plastic_nocontrol |>
  group_by(material_class) |> 
  summarise(count = round(sum(count), 0)) |>
  ungroup()

print(sum(material_plastic$count) == sum(sample_plastic$`Particle Count`))
# clean environment
# vars_keep <- c("plastic_count", "polymer_count", "MIPPR_breakdown", "filter_data", "project_name", "local_store_results")
# rm(list = setdiff(ls(), vars_keep))

# ** Need to find a better way to check for similarities of material class to group them

