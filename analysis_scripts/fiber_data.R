library(tidyverse)
library(readxl)
library(openxlsx)

# Read and merge spreadsheets
fiber_data <- loadWorkbook("data_cleaning/fiber_data/SFEI_01_FIBER_Meas.xlsx")

multiplier <- readxl::read_xlsx("data_cleaning/SFEI_01_Multiplier.xlsx", sheet = "Fragments_Fibers") |>
  clean_names()



df <- fiber_data$sheet_names |> 
  lapply(function(x) read.xlsx(fiber_data, sheet = x)) |> 
  do.call(what=rbind)

# Data cleaning total
fiber_df <- df |> 
  mutate(
    sample_id = 
      str_replace(
        Source, 
        "^([^_]+_[^_]+)_(S|W)(\\d+)_.*", 
        "\\1_\\2\\3"
  )) |> 
  filter(grepl("SFEI", sample_id)) |> 
# Variation of fiberpolymer type
  mutate(FiberPolymerType =
           case_when(
             grepl("Cellulose", FiberPolymerType) ~ "Cellulose Derivative",
             grepl("Organic", FiberPolymerType) ~ "Organic",
             grepl("Polycrylon", FiberPolymerType) | grepl("Polyacry", FiberPolymerType) ~ 
               "Polyacrylonitrile",
             grepl("Polypro", FiberPolymerType) | grepl("polypro", FiberPolymerType)~ "Polypropylene",
             grepl("Un", FiberPolymerType) | grepl("un", FiberPolymerType) ~ "Unknown",
             grepl("polyester", FiberPolymerType) ~ "Polyester",
             grepl("Polya", FiberPolymerType) ~ "Polyacrylonitrile",
             TRUE ~ FiberPolymerType
           )
           )


fiber_total <- fiber_df |> 
  group_by(sample_id) |> 
  summarize(total_count = n()) 

# fiber_plastics 
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
    .cols = where(is.numeric) & !all_of("multiplier"),
    ~ .x/multiplier)
    ) |> 
  select(-c(proportions_of_sample, multiplier)) |> 
  mutate(
    total_count = rowSums(select(where(is.numeric)), na.rm = TRUE)
  )
  




