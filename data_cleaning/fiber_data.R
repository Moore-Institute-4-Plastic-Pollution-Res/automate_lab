
# Load in fiber data and merge 
# df <- fiber_data$sheet_names |> 
#   lapply(function(x) read.xlsx(fiber_data, sheet = x)) |> 
#   do.call(what=rbind)

# Data cleaning total
fiber_df <- fiber_data |> 
  #filter(grepl("JHENG", sample_id)) |> 
  # Variation of Fiber Polymer Type 
  mutate(polymer_id = tolower(polymer_id),
         polymer_id = 
           case_when(
             grepl("cellulose", polymer_id, ignore.case = T) ~ "cellulose derivative (ether cellulose)",
             grepl("organic", polymer_id, ignore.case = T) ~ "organic matter",
             grepl("(polycrylon|polya)", polymer_id, ignore.case = T) ~ 
               "polyacrylonitrile",
             grepl("polypro", polymer_id, ignore.case = T) ~ "polypropylene",
             grepl("un", polymer_id, ignore.case = T) ~ "unknown",
             grepl("polyester|terephthalates", polymer_id, ignore.case = T) ~ "poly(esters/ethers/diglycidylethers/terephthalates)s",
             TRUE ~ polymer_id
           )
  ) |> 
  # clean extra text in ID
  mutate(sample_id = word(sample_id, 1)) |>
  mutate(particle_color = case_when(
    grepl("^(red)$", particle_color, ignore.case = T)  ~ "Red",
    #white
    grepl("^(Dk Blue|Blue)$", particle_color, ignore.case = T)  ~ "Blue",
    #Blue
    grepl("^(clear)$", particle_color, ignore.case = T) ~ "White",
    # Black
    grepl("^(yellow)$", particle_color, ignore.case = T) ~ "Orange",
    TRUE ~ particle_color
  )) 

#Check results
unique(fiber_df$polymer_id)
unique(fiber_df$particle_shape)
unique(fiber_df$particle_color)

write.csv(fiber_df,
          "data/fiber_data_full_final.csv", 
          row.names = F)

# Total count of measured object in each sample
fiber_total <- fiber_df |> 
  group_by(sample_id) |> 
  summarize(total_count = n()) 

# Total count of identified plastic in each sample 
fiber_type <- fiber_df |> 
  group_by(sample_id,FiberPolymerType) |> 
  mutate(
    type = 
    case_when(
      grepl("P", FiberPolymerType) ~ "Plastic",
      TRUE ~ FiberPolymerType
    )
  ) |> 
  group_by(sample_id,type) |> 
  summarize(count_type = n()) |> 
  pivot_wider(names_from = type, values_from = count_type) |> 
  left_join(multiplier) |> 
  mutate(across(
    .cols = where(is.numeric) & !all_of(c("subsample_ratio", "multiplier")),
    ~ (.x/subsample_ratio)/multiplier)
    ) |> 
  select(-c(proportions_of_sample, multiplier)) |> 
  mutate(
    total_count = rowSums(across(where(is.numeric)), na.rm = TRUE))

fiber_particle_count <- fiber_type |> 
  select(sample_id, Plastic) |> 
  mutate(Plastic = floor(Plastic)) |> 
  rename(SampleID = sample_id,
         "Particle Count" = Plastic
         ) 

# Polymer type count
fiber_breakdown <- fiber_df |>
  group_by(sample_id, FiberPolymerType) |>
  summarize(count = n()) |>
  left_join(multiplier) |>
  mutate(across(
    .cols = where(is.numeric) &
      !all_of(c("subsample_ratio", "multiplier")),
    ~ (.x / subsample_ratio) / multiplier
  )) |>
  select(-c(proportions_of_sample, multiplier)) |>
  group_by(FiberPolymerType) |>
  summarize(count = sum(count)) |>
  # remove non plastics
  filter(!(
    FiberPolymerType %in% c("Organic", "Unknown")
  )) |>
  # format polymer name as fragment polymer material class
  mutate(
    FiberPolymerType =
      case_when(
        grepl("Poly", FiberPolymerType) &
          !(str_detect(FiberPolymerType, "Polyvinylalcohols")) ~
          str_replace(FiberPolymerType, "^Poly([A-Za-z]+)$", "poly(\\1)"),
        TRUE ~ FiberPolymerType
         )) |>
  rename(material_class = FiberPolymerType)

# Hard coded changes
fiber_breakdown <- fiber_breakdown |> 
  mutate(material_class = 
           case_when(
             material_class == "poly(ester)" ~ 
               "poly(esters/ethers/diglycidylethers/terephthalates)s",
             material_class == "poly(acrylonitrile)" ~
               "polyacrylonitriles (nitriles)",
             TRUE ~ material_class
           ),
         count = floor(count),
         material_class = tolower(material_class)
           )
