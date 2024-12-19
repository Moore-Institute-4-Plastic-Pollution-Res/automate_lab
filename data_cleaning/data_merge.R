# Merge plastic count data for fragment, fiber, and filter 
# Read in data
source("data_cleaning/fragment_merge.R")
#source("data_cleaning/fiber_data.R")
#source("data_cleaning/filter_500.R")

# Plastic count merge

all_data <- bind_rows(filter_data |> mutate(particle_id = as.character(particle_id)), 
                      fragment_data1|> mutate(max_length_um= as.numeric(max_length_um),
                                              min_length_um= as.numeric(min_length_um)
                                              ))

fwrite(all_data, "full_particle_data.csv")

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
    )) |>
  filter(!str_detect(sample_id, "MIPPR")) |> 
  mutate(
    sample_id = str_replace_all(sample_id, "_REDO",""),
    sample_id = str_replace_all(sample_id, "O", "0"),
    sample_name = str_replace(sample_id, "^([^_]+_[^_]+)_(S|W)(\\d+)_.*", "\\1_\\2\\3")
  ) |>
  mutate(count = 1) |>
  left_join(multiplier2) |>
  left_join(multiplier, by = "sample_id")

#assuming that the large and small particle multipliers have different IDs here. 
sample_plastic <- confident_plastic |>
  group_by(sample_id) |> 
  summarize(count = ifelse(is.na(multiplier.x), sum(count) * multiplier.x,  sum(count) * multiplier.y)) |> 
  rename(SampleID = sample_name,
         "Particle Count" = count
  )

#assuming that the large and small particle multipliers have different IDs here. 
material_plastic <- confident_plastic |>
  group_by(material_class) |> 
  summarize(count = ifelse(is.na(multiplier.x), sum(count) * multiplier.x,  sum(count) * multiplier.y)) |> 
  rename(SampleID = sample_name,
         "Particle Count" = count
  )

# Merge polymer count data for fragment, fiber, and filter ----
# Read in data
polymer_count <- rbind(individual_breakdown, filter_breakdown) |> 
  group_by(material_class) |> 
  summarize(count = sum(count)) |> 
  mutate(percent = round((count/(sum(count, na.rm = TRUE)) * 100),1)
         ) |> 
  arrange(desc(percent)) |> 
  na.omit()

# clean environment
# vars_keep <- c("plastic_count", "polymer_count", "MIPPR_breakdown", "filter_data", "project_name", "local_store_results")
# rm(list = setdiff(ls(), vars_keep))

# ** Need to find a better way to check for similarities of material class to group them

