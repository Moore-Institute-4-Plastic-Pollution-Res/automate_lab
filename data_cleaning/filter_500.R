library(tidyverse)
library(readxl)
library(stringr)

filter_data <- read_csv("data_cleaning/filter_500um/filter_500_full_data.csv")
multiplier <- readxl::read_xlsx("data_cleaning/SFEI_01_Multiplier.xlsx", sheet = "Filter") |>
  clean_names()

# total count to then focus on the good spectra and plastics
total_count <- filter_data |> 
  group_by(sample_id) |> 
  summarize(total = n())

# good spectra and plastic
total_plastic <- filter_data |>
  filter(bad_spectra == "TRUE",
         max_cor_val > 0.6,
         !(
           material_class %in% c(
             "other material",
             "mineral",
             "organic matter",
             "cellulose derivatives (ether cellulose)"
           )
         )) |>
  mutate(
    sample_id = str_replace_all(sample_id, "O", "0"),
    sample_name = str_replace(sample_id, "^([^_]+_[^_]+)_(S|W)(\\d+)_.*", "\\1_\\2\\3")
    ) 

# Particle count
particle_count <- total_plastic |> 
  group_by(sample_name) |> 
  summarize(count = n()) |> 
  left_join(multiplier, by = c("sample_name" = "sample_id")) |> 
  mutate(across(.cols = where(is.numeric) & !all_of("sample_multiplier"), ~ .x /
                  sample_multiplier),
         count = floor(count)
         ) |> 
  select(-sample_multiplier) |> 
  rename(SampleID = sample_name,
         "Particle Count" = count
         )

write.csv(particle_count, "data_cleaning/final/SFEI_filter_plastic_count.csv")

# Particle count by polymer
polymer_count <- total_plastic |> 
  group_by(sample_name,material_class) |> 
  summarize(count = n()) |> 
  left_join(multiplier, by = c("sample_name" = "sample_id")) |> 
  mutate(across(.cols = where(is.numeric) & !all_of("sample_multiplier"), ~ .x /
                  sample_multiplier),
         count = floor(count)
  ) |> 
  group_by(material_class) |> 
  summarize(count = sum(count))

write.csv(polymer_count, "data_cleaning/final/SFEI_filter_polymer_count.csv")

