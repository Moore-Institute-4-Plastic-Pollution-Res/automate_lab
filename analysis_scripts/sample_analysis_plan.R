library(tidyverse)
library(readxl)
library(dplyr)
library(janitor)

# Data
labguru_df <- read_xlsx("EXAMPLE_01 LabGuruSampleWorksheet.xlsx", sheet = 1) |> 
  clean_names()

labguru_df <- read_xlsx("JHENG_270_LabGuruSampleWorksheet.xlsx", sheet = 1) |> 
  clean_names()


# Extract distinct sieves used as a string---
sieves <- labguru_df |> 
  select(pre_filtered) |> 
  distinct() |>
  drop_na() |> 
  # if combined sieve values in a row
  separate_rows(pre_filtered, sep = ', ') |> 
  mutate(values = str_sub(`pre_filtered`, end = -3),
         values = str_sub(`values`, start = 2),
         values = paste0(values, "µm")
         ) |> 
  distinct(values) 

if (length(sieves$values) > 2) {
  sieves_string <- sieves |> 
    pull(values) |> 
    as.character() |> 
    paste(collapse = ", ")
} else if (length(sieves$values) == 2) {
  sieves_string <- sieves |> 
    pull(values) |> 
    as.character() |> 
    paste(collapse = " and ")
} else {
  sieves_string <- ""
}

# Digestion and exceptions ---- assume each sample gets the same concentration for now
digestion <- labguru_df |> 
  select(digestion_concentration, digestion_solution_composition) |> 
  mutate(digestion_solution_composition = tolower(digestion_solution_composition)) |> 
  distinct() |> 
  drop_na() 

# digestion concentration multiple option
if (length(digestion$digestion_concentration) > 2){
  max <- max(digestion$digestion_concentration)
  min <- min(digestion$digestion_concentration)
  solution <- unique(digestion$digestion_solution_composition)
  digestion <- paste0(
    min*100, "% to ",
    max*100, "% ", solution
  )
} 

# digestion exception 
digestion_exception <- labguru_df |> 
  filter(is.na(digestion_solution_composition),
         is.na(digestion_concentration),
         !is.na(sample_id)
         ) |> 
  pull(sample_id) |> 
  as.character()

# Density solutions ----
density <- labguru_df |> 
  select(density_separation_density, density_separation_composition) |> 
  distinct() |> 
  drop_na()

# dependent on the chemical - assuming we always use CaCl2
if (density$density_separation_composition == "CaCl2") {
  density <- density |>
    mutate(
      value = paste0(
        density_separation_density,
        " g/ml density solution of calcium chloride (",
        density_separation_composition,
        ")"
      )
    ) |> 
    pull(value)
} else {
  density <- density |>
    mutate(
      value = paste0(
        density_separation_density,
        " g/ml density solution of ",
        density_separation_composition
      )
    ) |> 
    pull(value)
}

# filter 
filter <- labguru_df |> 
  select(filter_diameter, filter_pore_size) |> 
  distinct() |> 
  drop_na() |> 
  mutate(value = paste0(filter_diameter, 
                        " stainless steel filters with ", 
                        filter_pore_size)) |> 
  pull(value)

# added volume for sieving
volume_added <-labguru_df |> 
  select(SampleID, `Volume solution added`)




  


