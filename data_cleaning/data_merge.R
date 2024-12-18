# Merge plastic count data for fragment, fiber, and filter 
# Read in data
source("data_cleaning/fragment_merge.R")
source("data_cleaning/fiber_data.R")
source("data_cleaning/filter_500.R")

# Plastic count merge
plastic_df <- rbind(fiber_particle_count, fragment_particle_count, filter_particle_count) 
plastic_count <- plastic_df |> 
  filter(!str_detect(SampleID, "MIPPR")) |> 
  group_by(SampleID) |> 
  summarize(`Particle Count` = sum(`Particle Count`))

# Merge polymer count data for fragment, fiber, and filter ----
# Read in data
polymer_count <- rbind(fiber_breakdown, fragment_breakdown, filter_breakdown) |> 
  group_by(material_class) |> 
  summarize(count = sum(count)) |> 
  mutate(percent = round((count/(sum(count)) * 100),1)
         ) |> 
  arrange(desc(percent))

# set sample ID with particle count
plastic_count <- left_join(sampleid, plastic_count, by = c("mippr_sample_id" = "SampleID")) |> 
  select(-mippr_sample_id) |> 
  drop_na()

# clean environment
vars_keep <- c("plastic_count", "polymer_count", "MIPPR_breakdown", "filter_data")
rm(list = setdiff(ls(), vars_keep))

# ** Need to find a better way to check for similarities of material class to group them

