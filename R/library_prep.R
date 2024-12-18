# Set library ----
lib <- OpenSpecy::load_lib("derivative")
lib <- filter_spec(lib, lib$metadata$spectrum_type == "ftir") #This improves accuracy 10%

# Material metadata rename and group similar material class
lib$metadata <- lib$metadata %>%
  mutate(material_class = ifelse(
    material_class %in% c("polyacrylamides", "polyamides (polylactams)"),
    "poly(acrylamide/amid)s",
    material_class
  )) %>%
  mutate(
    material_class = ifelse(
      material_class %in% c(
        "polyterephthalates",
        "polyesters",
        "polyethers (polyglycols, oxyalkylenes, glycidyl ethers & cyclic ethers)",
        "polydiglycidyl ethers (polyepoxides, polyhydroxyethers, phenoxy)"
      ),
      "poly(esters/ethers/diglycidylethers/terephthalates)s",
      material_class
    )
  )
