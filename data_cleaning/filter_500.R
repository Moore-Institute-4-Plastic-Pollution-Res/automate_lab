library(tidyverse)
library(readxl)
library(stringr)

filter_data <- read_csv("data_cleaning/data/filter_500_full_data.csv")
multiplier <- readxl::read_xlsx("data_cleaning/data/SFEI_01_Multiplier.xlsx", sheet = "Filter") |>
  clean_names() |> 
  mutate(
    sample_id = str_replace_all(sample_id, "O", "0"))


# total count to then focus on the good spectra and plastics
total_count <- filter_data |> 
  group_by(sample_id) |> 
  summarize(total = n())

# good spectra and plastic
total_plastic <- filter_data |>
  filter(bad_spectra != "TRUE",
         max_cor_val >= 0.6,
         !(
           material_class %in% c(
             "other material",
             "mineral",
             "organic matter"
           )
         )) |>
  mutate(
    sample_id = str_replace_all(sample_id, "_REDO",""),
    sample_id = str_replace_all(sample_id, "O", "0"),
    sample_name = str_replace(sample_id, "^([^_]+_[^_]+)_(S|W)(\\d+)_.*", "\\1_\\2\\3")
    ) 

# Particle count
filter_particle_count <- total_plastic |> 
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

# Particle count by polymer (exclude MIPPR)
filter_breakdown <- total_plastic |>
  filter(!str_detect(sample_id, "MIPPR")) |> 
  group_by(sample_name,material_class) |> 
  summarize(count = n()) |> 
  left_join(multiplier, by = c("sample_name" = "sample_id")) |> 
  mutate(across(.cols = where(is.numeric) & !all_of("sample_multiplier"), ~ .x /
                  sample_multiplier),
         count = floor(count)
  ) |> 
  group_by(material_class) |> 
  summarize(count = sum(count)) |> 
  mutate(material_class =
           case_when(
             str_detect(material_class, "_") ~ str_replace_all(material_class, "_","/"),
             TRUE ~ material_class
           ))

# MPIRR breakdown of polymer type
MIPPR_breakdown <-  total_plastic |>
  filter(str_detect(sample_id, "MIPPR")) |> 
  group_by(sample_name,material_class) |> 
  summarize(count = n()) |> 
  left_join(multiplier, by = c("sample_name" = "sample_id")) |> 
  mutate(across(.cols = where(is.numeric) & !all_of("sample_multiplier"), ~ .x /
                  sample_multiplier),
         count = floor(count)
  ) |> 
  select(-sample_multiplier)
