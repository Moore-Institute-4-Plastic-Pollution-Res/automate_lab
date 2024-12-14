library(tidyverse)
library(readxl)
library(openxlsx)
library(janitor)

# Read and merge spreadsheets
fiber_data <- loadWorkbook("data_cleaning/fiber_data/SFEI_01_FIBER_Meas.xlsx")

multiplier <- readxl::read_xlsx("data_cleaning/SFEI_01_Multiplier.xlsx", sheet = "Fragments_Fibers") |>
  clean_names()

# Load in fiber data and merge 
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
# Variation of Fiber Polymer Type 
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

fiber_type <- fiber_type |> 
  select(sample_id, Plastic) |> 
  mutate(Plastic = floor(Plastic)) |> 
  rename(SampleID = sample_id,
         "Particle Count" = Plastic
         ) 

write.csv(fiber_type, "data_cleaning/final/SFEI_fiber_plastic_count.csv")    


# Polymer type count
fiber_breakdown <- fiber_df |> 
  group_by(sample_id, FiberPolymerType) |> 
  summarize(count = n()) |> 
  left_join(multiplier) |> 
  mutate(across(
    .cols = where(is.numeric) & !all_of(c("subsample_ratio","multiplier")),
    ~ (.x/subsample_ratio)/multiplier)
  ) |> 
  select(-c(proportions_of_sample, multiplier)) |> 
  group_by(FiberPolymerType) |> 
  summarize(count = sum(count)) |> 
  # remove non plastics
  filter(!(FiberPolymerType %in% c("Cellulose Derivative", "Organic", "Unknown"))) |> 
  # format polymer name as fragment polymer material class
  mutate(FiberPolymerType = 
           case_when(
             grepl("Poly", FiberPolymerType) & !(str_detect(FiberPolymerType,"Polyvinylalcohols")) ~ 
               str_replace(FiberPolymerType, "^Poly([A-Za-z]+)$", "poly(\\1)"),
             TRUE ~ tolower(FiberPolymerType)
           )
           ) |> 
  mutate(count = floor(count))

write.csv(fiber_breakdown, "data_cleaning/final/SFEI_fiber_polymer_count.csv")  
