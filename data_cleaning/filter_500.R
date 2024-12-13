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
         !(material_class %in% c("other material", "mineral", 
                                 "organic matter", "cellulose derivatives (ether cellulose)"))
         ) |> 
  group_by(sample_id) |> 
  summarize(plastic_count = n())

filter_data_df <- left_join(total_count, total_plastic) |> 
  mutate(
    sample_name = str_replace(
      sample_id, 
      "^([^_]+_[^_]+)_(S|W)(\\d+)_.*", 
      "\\1_\\2\\3"
    )
  ) |> 
  group_by(sample_name) |>
  summarise(across(where(is.numeric), ~sum(.x,na.rm = TRUE))) |> 
  rename(sample_id = sample_name) |> 
  # modify the O as 0
  mutate(sample_id = 
           str_replace_all(sample_id, "O", "0")
           ) |> 
  left_join(multiplier) |> 
  mutate(total_multi = plastic_count/sample_multiplier)

write.csv(filter_data_df, "data_cleaning/final/SFEI_filter_final.csv")

# Particle count (plastic count)
particle_count <- filter_data |> 
  filter(bad_spectra == "TRUE",
         !(material_class %in% c("other material", "mineral", 
                                 "organic matter", "cellulose derivatives (ether cellulose)"))
  )
