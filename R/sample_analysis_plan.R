
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

if (length(digestion$digestion_concentration) > 2) {
  max <- max(digestion$digestion_concentration)
  min <- min(digestion$digestion_concentration)
  solution <- unique(digestion$digestion_solution_composition)
  digestion <- paste0("digested in ", min * 100, "% to ", max * 100, "% ", solution)
} else if(is.na(digestion)){
    digestion <- "not digested "
  } else{
  digestion <- paste0("digested in ",
    digestion$digestion_concentration * 100,
    "% ",
    digestion$digestion_solution_composition
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
if (nrow(density)!= 0) {
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
} else {
    density <- "not density separated"
}

# filter 
filter <- labguru_df |> 
  select(filter_diameter, filter_pore_size) |> 
  distinct() |> 
  drop_na() |> 
  mutate(
    filter_pore_size = str_replace_all(filter_pore_size, "u", "µ"),
    value = paste0(filter_diameter, 
                        " stainless steel filters with ", 
                        filter_pore_size)) |> 
  pull(value)

# added volume for sieving
if(is.na(unique(labguru_df$volume_solution_added))){
  volume_added = "alcojet"
} else{
  volume_added <- labguru_df |> 
    select(sample_id, volume_solution_added) |> 
    filter(!str_detect(sample_id, "MIPPR")) |> 
    distinct(volume_solution_added) |> 
    separate(
      volume_solution_added,
      into = c("a", "b","c","d", "e", "f", "g"),
      sep = " "
    ) |> 
    mutate(
      value = paste(a, b, "of", c, "and", e, f, "of", g)
    ) |> 
    pull(value)  
}


## ------------------------- QA/QC --------------------------

# Procedural Blank
# MIPPR Procedual blanks
# plastic_summary_all for this
mippr_pb <- confident_plastic |> 
  filter(str_detect(sample_name, "PB")) |> 
  group_by(material_class) |> 
  summarize(count = sum(count)) |> 
  mutate(percent = (round(count/sum(count) * 100, 1))) |> 
  arrange(desc(percent))

# LFB table prep - assuming one LFB present
mippr_lfb_name <- confident_plastic |> 
  filter(str_detect(sample_id, "LFB")) |> 
  distinct(sample_id) |> 
  pull(sample_id)
  
mippr_lfb <- confident_plastic |>
  filter(str_detect(sample_id, "LFB")) |> 
  group_by(material_class) |>
  summarize(count = sum(count)) |> 
  mutate(type_of_plastic = 
           case_when(
             material_class == "poly(ethylene)" | 
               material_class == "poly(propylene)" ~ "Polyolefins",
             TRUE ~ material_class
           )
           ) |> 
  group_by(type_of_plastic) |> 
  summarize(post_count = sum(count)) |> 
  mutate(type_of_plastic = str_to_title(type_of_plastic)) 

lfb <- read_csv("data/LFB_Counts.csv") |> 
  rename(`Pre-Count` = Count,
         `Type of Plastic` = Type,
         `Size Class` = Size
         ) |> 
  clean_names()

lfb_df <- left_join(lfb, mippr_lfb) |>
  mutate(
    across(everything(), ~replace(., is.na(.), 0)),
    percent = round((post_count / pre_count) * 100, 1))

# spike recovery
pa66 <- lfb_df |> 
  filter(type_of_plastic == "Polyamides") |> 
  pull(percent)

cellu <- lfb_df |> 
  filter(type_of_plastic == "Cellulose Derivatives (Ether Cellulose)") |> 
  pull(percent)

poly <- lfb_df |> 
  filter(type_of_plastic == "Polyolefins") |> 
  pull(percent)

# Is the recovery acceptable or not
recovery  <-
  lfb_df |>
  mutate(check =
           case_when(percent >= 50 & percent <= 150 ~ "TRUE",
                     TRUE ~ "FALSE"
                       ))
if (any(str_detect(recovery$check, "FALSE"))) {
  acceptable <- "not within the acceptable range"
} else{
  acceptable <- "within the acceptable range"
}

# Which of the recoveries is not acceptable
recover_na <- 
  recovery |> 
  filter(check == "FALSE") |> 
  mutate(text = paste0(
    type_of_plastic, " in the ", size_class
  )) |>  pull()

if (length(recover_na) > 1){
  recover <- paste0(recover_na, collapse = " and ")
} else{
  recover <- recover_na
}

# Notes
notes <- all_data |>
  filter(!is.na(notes) & notes != "") |>
  select(sample_id, notes, particle_id) |>
  distinct() |>
  mutate(text = paste(sample_id, particle_id, notes)) |>
  pull(text) |>
  paste(collapse = " and ")
