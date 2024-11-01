library(OpenSpecy)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(data.table)
library(mmand)
library(magick)
library(cluster)
library(googledrive)
library(purrr)

# Load all R scripts in the folder
walk(list.files(path = "R", pattern = "\\.R$", full.names = TRUE), source)

#Run once
# check if it has been installed or not 
get_lib(c("mediod", "derivative"))

#------------------- Google Drive File Download ------------------------------
# Authorize Google Drive Connection 
drive_auth(cache = ".secrets/")

# Ask User the name of the folder they want to access
project_name <- readline(prompt = "Enter the Name of the Google Drive Folder for this Project: ")

# Access ID of all files in each folder --------------------------------------
folders_of_interest <- c("Mosaic_Images", "Export_Files")

# Create folders to store data any project  name


# ----------------------------- Functions ---------------------------------
particle_image <- function(proc_map_metadata, map_metadata, material_class, pixel_length = 25, origin, title){
    ggplot(data = proc_map_metadata, aes(x = centroid_x*pixel_length + origin[1], y = centroid_y*pixel_length + origin[2], label = feature_id)) +
        geom_raster(data = map_metadata, aes(x = x*pixel_length + origin[1], y = y*pixel_length + origin[2], fill = .data[[material_class]])) +
        geom_text_repel(size = 1.5, alpha = 0.5) +
        coord_fixed() +
        scale_fill_viridis_d() +
        theme_bw(base_size = 5) +
        labs(x = "X (um)", y = "Y (um)", title = title) +
        #guides(fill = guide_legend(override.aes = list(size = 5)))  +
        theme(legend.position = "none", plot.title = element_text(size=5))
}

determine_thresholds <- function(x, features, test_values, true_count){
    values <- function(y){
        def_features(x, features = features > y)[["metadata"]][["feature_id"]] |>
            unique() |>
            length() - 1
    }  
    
   plots <- ggplot()+
        geom_line(aes(x = test_values, y = values)) +
        geom_hline(aes(yintercept = true_count)) +
        #scale_x_log10() +
        scale_y_log10() +
        theme_bw()
   
    error <- values - true_count
    
    ideal_threshold <- test_values[which.min(abs(error))]
    
    ideal_threshold_error <- error[which.min(abs(error))]
    
    list(values, plots, ideal_threshold, ideal_threshold_error)
}

########### -----------------------------------
walk(list.files(path = "R", pattern = "\\.R$", full.names = TRUE), source)



#Spectral Library ----
#get_lib("derivative") #Run your first time. 
#"C:/Users/winco/OneDrive/Documents/Positive_Controls" controls to fine tune
#If finding poor id accuracy, can try removing other plastic from filter to see if the material is a known novel type.
#path = "C:\\Users\\winco\\OneDrive\\Documents\\OpenSpecy_offline", 
lib <- OpenSpecy::load_lib(type = "derivative") 
lib <- filter_spec(lib, 
                   #!lib$metadata$material_class %in% c("other plastic", 
                    #                                   "other material"#, 
                    #                                   #"mineral", 
                                                       #"organic matter"
                     #                                  ) & 
                       lib$metadata$spectrum_type == "ftir") #This improves accuracy 10%

# Material metadata rename and group similar material class
lib$metadata <- lib$metadata %>%
    mutate(material_class = ifelse(material_class %in% c("polyacrylamides", "polyamides (polylactams)"), "poly(acrylamide/amid)s", material_class)) %>%
    mutate(material_class = ifelse(material_class %in% c("polyterephthalates", "polyesters", "polyethers (polyglycols, oxyalkylenes, glycidyl ethers & cyclic ethers)", "polydiglycidyl ethers (polyepoxides, polyhydroxyethers, phenoxy)"), "poly(esters/ethers/diglycidylethers/terephthalates)s", material_class)) 

#Set this ----
#You'll set this wherever you put your zip folders. C:/Users/winco/OneDrive/Documents/Positive_Controls
wd <- file.path(project_name, "Spectral_Results")

#Unzip files in batch ----
# zip_files <- list.files(path = wd,
#            pattern = "\\.zip$",
#            full.names = T)
# 
# for(zip_file in zip_files){
#     print(zip_file)
#     unzip(zipfile = zip_file,exdir=wd)  # unzip your file
# }

# check images in folder
# list.files(path = wd,
#            pattern = ".JPG$", 
#            full.names = TRUE)

#Remove already processed data
#files <- files[!file.exists(paste0(dirname(files), "/particles_", gsub("\\.dat$", ".rds", basename(files))))]

#files <- files[11]
#files <- files[grepl("GL", files)]
#origins <- list(x = c(-3476, 4324), y = c(-7289, -7217))#list(x = c(-1866, 0), y = c(-9987, 0))
#transform <- c(1:4, 14, 16,17)
#files <- list.files(path = wd, pattern = ".dat", full.names = T)
#for(file in files){
#    print(file)
#    read_envi(file) |>
#        adj_intens(type = "transmittance", make_rel = F) |>
#        write_spec(gsub(".dat", ".rds", file))
#}

#This is the main function. 

#debugging
#debug(analyze_features)
#options(error = browser)
#files = files[12]

start = Sys.time()

analyze_features(drive_name = "Projects", 
                     project_name = project_name,
                     folders_of_interest = folders_of_interest,
                     lib = lib,
                     #img = img,#"C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\Mosaic Image files\\Recovery_red_beads_75-90um_5um-screen.JPG",
                     #bottom_left = list(c(171, 472)),
                     #top_right = list(c(632, 10)),
                     #origins = NULL, 
                     spectral_smooth = F, 
                     #sigma1 = c(0.001,2),
                     #sigma2 = c(0.001,2),
                     #spatial_smooth = ,
                     close = F,
                     adj_map_baseline = F, #Probably needs to be updated somehow because currently very small change and decreasing id accuracy. 
                     sn_threshold = 0.01, #0.01 the lowest signal considered "particles"
                     cor_threshold = 0.66, #the lowest level you would consider something known
                     area_threshold = 1, #to remove artefacts from background. 
                     label_unknown = F, #whether to scrub information when below correlation threshold.
                     remove_nonplastic = F, #whether to remove all nonplastic ided particles
                     remove_unknown = F,
                     pixel_length = 25, #25 conversion from pixels to length, set up for microns
                     metric = "sig_times_noise", #technique for the signal
                     abs = F, #F This determines whether to take the absolute value of the metric.
                     particle_id_strategy = "collapse", #This describes how the matching happens. 
                     vote_count = 10,
                     collapse_function = median,
                     k = 1,
                     k_weighting = "mean", #or multiple
                     wd = paste0(project_name,"/Spectral_Results"), #will put results here. 
                     types = c(
                        "particle_heatmap_thresholded",
                                "particle_heatmap", 
                               "particle_details",
                               "particle_summary",
                                "particle_image",
                                "spectra_processed"
                               ),
                     by = c("sample", 
                            "all"),
                     width = 1000, 
                     height = 1000, 
                     units = "px")


end.time <- Sys.time()
time.taken <- end.time - start.time



    # All study figures ----
    study_results <- fread(paste0(wd, "/particle_details_all.csv")) %>%
        mutate(sample_id = gsub("([0-9]{1}of[0-9]{1}_.*)|(_lft)|(_rt)", "", sample_id))

    #Crowding factor
    study_results %>%
        group_by(sample_id) %>%
        summarise(crowd_factor = sum(area_um2[sqrt(area_um2) > 500])/sum(area_um2))
    
    #Percent comp
    study_results_percent <- study_results %>%
        group_by(sample_id, material_class) %>%
        summarise(count = n()) %>%
        ungroup() %>%
        group_by(sample_id) %>%
        mutate(percent = count/sum(count) * 100) 
    
    all <- study_results_percent %>%
        group_by(material_class) %>%
        summarise(count = sum(count)) %>%
        ungroup() %>%
        mutate(percent = count/sum(count) * 100) %>%
        mutate(sample_id = "All") %>% 
        as.data.frame()
    
    study_results_percent_all <- bind_rows(study_results_percent, all)
    
    ggplot(study_results_percent_all, aes(x = percent, y = material_class, fill = material_class)) +
        geom_bar(stat = "identity") +
        facet_wrap(sample_id ~.) +
        theme_bw(base_size = 8) +
        geom_text(aes(label = count), size = 2, hjust = 0) +
        scale_fill_viridis_d() +
        theme(legend.position = "none")
    
    
    study_results %>%
        rbind(study_results %>% mutate(sample_id = "All")) %>%
        ggplot(aes(x = sqrt(area_um2), color = material_class)) +
        stat_ecdf() +
        facet_wrap(sample_id ~ .) +
        scale_x_log10() +
        scale_color_viridis_d() +
        theme_bw(base_size = 8) +
        labs(x = "Nominal Particle Size", y = "Proportion Smaller")
    
    study_results %>%
        rbind(study_results %>% mutate(sample_id = "All")) %>%
        ggplot(aes(x = aspect_ratio, y = material_class)) +
        geom_boxplot() +
        facet_wrap(sample_id ~ .) +
        scale_color_viridis_d() +
        theme_bw(base_size = 8) +
        labs(x = "Aspect Ratio", y = "Material Class")
    
    #correlation between known and unknown
    correlation_test <- study_results_percent %>%
        pivot_wider(names_from = material_class, values_from = count, id_cols = sample_id) %>%
        as.data.table()
    
    correlation_test2 <- correlation_test[,-"sample_id"][,lapply(.SD, mean_replace)]
    cor_matrix <- cor(correlation_test2)
    #This shouldn't have a positive correlation to all plastics, would suggest that 
    #The plastics are just due to having a large number of particles. 
    
    #Check for most common id'ed spectra
    sort(table(study_results$max_cor_name)/sum(nrow(study_results)))
    
    #Manually assess any greater than 1 %
    check_these <- sort(table(study_results$max_cor_name)/sum(nrow(study_results)))[sort(table(study_results$max_cor_name)/sum(nrow(study_results))) > 0.01]

    lib_check <- filter_spec(lib, names(check_these))
    
    for(n in 1:ncol(lib_check$spectra)){
        plot(filter_spec(lib_check, n))
    }
#Compare with qa samples ----
threshold = 0.66

valid_results <- fread(paste0(wd, "/valid_results.csv")) %>%
    mutate(sample_id = gsub("FLB", "LFB", sample_id)) %>%
    filter(!grepl("PrcdrlBlnk", sample_id))
    

actual_results <- fread(paste0(wd, "/particle_details_all.csv")) %>%
    mutate(sample_id = gsub("(_lft$)|(_rt$)", "", sample_id)) %>%
    #filter(max_cor_val >= threshold) %>% #Acceptable accuracy doing it like this too. 
    filter(!grepl("(PrcdrlBlnk)|(dup)", sample_id)) %>%
    group_by(sample_id, material_class) %>%
    summarise(count = n(), 
              area = sum(area_um2)) 

raw_results_to_share <- fread(paste0(wd, "/particle_details_all.csv")) %>%
    filter(!grepl("PrcdrlBlnk", sample_id)) %>%
    filter(!material_class %in% c("mineral", "organic matter", "unknown")) %>%
    rename(ParticleID = particle_id,
           SpectraMatchValue = max_cor_val,
           Polymer = material_class,
           Length = max_length_um, 
           Width = min_length_um) %>%
    mutate(MethodologyID = ifelse(Polymer == "cellulose derivatives (ether cellulose)", "sop.pdf_2", "sop.pdf_3"), 
           PhotoID = paste0(sample_id, ".JPG"), 
           SpectraID = paste0(sample_id, ".rds"), 
           FinalAnalysisDate = "2024-04-20",
           SampleID = gsub("(_lft$)|(_rt$)", "", sample_id)) %>%
    dplyr::select(ParticleID, MethodologyID, SampleID, PhotoID, SpectraMatchValue, FinalAnalysisDate, Polymer, Length, Width)

fwrite(raw_results_to_share, paste0(wd, "\\", "raw_particle_results_to_share.csv"))

joined <- left_join(valid_results, 
                    actual_results, 
                    by = c("material_class", 
                           "sample_id"))

joined2 <- joined %>%
    mutate(across(everything(), ~ifelse(is.na(.), 0,.)))%>%
    mutate(recovery_automated = (count+recovered_in_pre)/total_spiked) %>%
    mutate(avg_area_est = area/(count+recovered_in_pre), 
           ratio = in_ftir_view/count)

mean(joined2$recovery_automated, na.rm = T)

sd(joined2$recovery_automated, na.rm = T)/mean(joined2$recovery_automated, na.rm = T) * 100

#mean(joined2$recovery_manual, na.rm = T)

#sd(joined2$recovery_manual, na.rm = T)/mean(joined2$recovery_manual, na.rm = T) *100

#ratio auto count compared to in view
#mean(joined2$automated_accuracy, na.rm = T)
#sd(joined2$automated_accuracy, na.rm = T)/mean(joined2$automated_accuracy, na.rm = T) *100

#mean(joined2$recovery_just_in_view, na.rm = T)
#sd(joined2$recovery_just_in_view, na.rm = T)/mean(joined2$recovery_just_in_view, na.rm = T) *100

#average size to see if errors are from particle clumping
joined2 %>%
    group_by(material_class) %>%
    summarise(mean = sqrt(mean(avg_area_est, na.rm = T)))

joined2 %>%
    group_by(material_class) %>%
    summarise(mean_spiked = mean(total_spiked, na.rm = T))


#Per size, this is really what needs to be accurate for automation to be accepted. 
#Benchmark to beat, 65% Recovery for  CA with 26% RSD, 61% recovery for PE with 26% RSD. 
recovery <- joined2 %>%
    mutate(corrected_recovery = ifelse(material_class == "polyolefins (polyalkenes)", (count *  1.77 +recovered_in_pre)/total_spiked, recovery_automated), 
           visual_recovery = in_ftir_view/total_spiked) %>%
    group_by(material_class) %>%
    summarise(
        mean_spiked = mean(total_spiked), 
        mean = mean(recovery_automated, na.rm = T), 
              ratio_mean = mean(ratio),
              rsd = sd(recovery_automated, na.rm = T)/mean(recovery_automated, na.rm = T) * 100, 
              corrected_mean = mean(corrected_recovery), 
              corrected_rsd = sd(corrected_recovery, na.rm = T)/mean(corrected_recovery, na.rm = T) * 100, 
              visual_recovery_accuracy = mean(visual_recovery),
              visual_recovery_rsd = sd(visual_recovery)/mean(visual_recovery))

#Per size
#joined2 %>%
#    group_by(material_class) %>%
#    summarise(mean = mean(recovery_manual, na.rm = T), 
#              rsd = sd(recovery_manual, na.rm = T)/mean(recovery_manual, na.rm = T) * 100)

#MDA
blank_results_raw <- fread(paste0(wd, "/particle_details_all.csv")) %>%
    filter(grepl("LRB", sample_id)) %>%
    #filter(max_cor_val >= threshold) %>% #Acceptable accuracy doing it like this too. 
    mutate(sample_id = gsub("(lft)|(rt)|(a$)|(b$)", "", sample_id)) %>%
    mutate(area_bins = cut(sqrt(area_um2), c(0, 212, 500, 5000))) %>%
    group_by(sample_id, material_class, area_bins) %>%
    summarise(count = n()) %>%
    ungroup()

files[grepl("LRB", files)]

MDA <- function(x){
    (mean(x, na.rm = T) + 3 + 3.29*sd(x, na.rm = T)*sqrt(1+1/4)) |> as.integer()
}

BDL <- function(x){
    (x + 3 + 4.65*sqrt(x))|> as.integer()
}

#Manual MDA Calculation
count_0_212 <- c(0,0,0,11)
MDA(count_0_212)

count_500_5000 <- c(0,0,0,1)
MDA(count_500_5000)

count_0_212 <- c(0,0,0,11)
mean(count_0_212, na.rm = T) + 3 + 3.29*sd(count_0_212, na.rm = T)*sqrt(1+1/4)

unique(blank_results_raw$sample_id) |> length()

BDL <- blank_results_raw %>%
    filter(!material_class %in% c("mineral", "organic matter", "other materials", "unknown")) %>%
    group_by(area_bins, sample_id) %>%
    summarise(count = sum(count)) %>%
    mutate(BDL_FTIR = (count + 3 + 4.65*sqrt(count))|> as.integer())


#Current benchmark to beat, 20 particles
join_empty <- expand.grid(sample_id = unique(blank_results_raw$sample_id), 
            area_bins = unique(blank_results_raw$area_bins))

MDA <- blank_results_raw %>%
    filter(!material_class %in% c("mineral", "organic matter", "other material", "unknown")) %>%
    group_by(sample_id, area_bins) %>%
    summarise(count = sum(count)) %>%
    ungroup() %>%
    right_join(join_empty) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    group_by(area_bins) %>%
    summarise(MDA_FTIR = mean(count, na.rm = T) + 3 + 3.29*sd(count, na.rm = T)*sqrt(1+1/4))

fwrite(recovery, paste0(wd, "/recovery_results.csv"))
fwrite(joined2, paste0(wd, "/analysis_results.csv"))
fwrite(MDA, paste0(wd, "/mda_results.csv"))
fwrite(blank_results_raw, paste0(wd, "/blank_results_raw.csv"))

#Test Accuracy of Routine ----
#Benchmarks to beat
#Best is 0.04 with mean*sd and range restriction to avoid co2 and dead region
#Not avoiding the dead region for matching though, that seems to be useful. 
#Using full library, adding overlap reduction and redundancy reduction,
#and being a little more laxed on naming to not require highest granularity we can get additional gains. 
#Count accuracy mean 124, cv 0.37.
#ID accuracy benchmark to beat 97 mean, cv 0.03. 
#Area accuracy mean 88, cv 0.47.
#particles <- lapply(list.files("C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\Mosaic Image files", pattern = ".csv", full.names = T), 
#                    function(x){
#                        y = data.table::fread(x)
#                        y$name = gsub(".csv", "", basename(x))
#                        y
#                    }
#                    ) |>
#    data.table::rbindlist() 

#particles_summary <- particles %>%
#    group_by(name) %>%
#    summarise(median_area = median(Area), 
#           median_feret = median(Feret),
#           count = n())

#true_values <- fread(paste0(wd, "/True_Values.csv")) %>%
#    left_join(particles_summary, by = c("Map" = "name"))

#fwrite(true_values, "C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\True_Values2.csv")

true_values <- fread("C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\true_values2.csv")

for(row in 1:nrow(true_values)){
    print(true_values[[row,"Map"]])
    test_accuracy <- fread(paste0(wd, "/particle_details_", true_values[[row,"Map"]], ".csv"))
    if(!is.na(true_values[[row,"count"]]) & true_values[[row,"count"]] != ""){
        print("Count % of Total Spiked")
        count <- nrow(test_accuracy)/true_values[[row,"count"]] * 100
        true_values[row,"Count Accuracy"] <- count
        print(count) 
    }
    
    if(!is.na(true_values[[row,"SpecID"]]) & true_values[[row,"SpecID"]] != ""){
        print("Accurate SpecID % of Total Spiked")
        categories <- test_accuracy$material_class
        #remove uncertain particles 
        #categories <- categories[test_accuracy$max_cor_val > 0.7]
        id <- (1 - sum(!grepl(true_values[[row,"SpecID"]], categories))/length(categories))*100
        true_values[row,"SpecID Accuracy"] <- id
        print(id)  
    }
    
    if(!is.na(true_values[[row,"PNID"]]) & true_values[[row,"PNID"]] != ""){
        print("Accurate PNID % of Total Spiked")
        categories <- test_accuracy$material_class
        #remove uncertain particles 
        #categories <- categories[test_accuracy$max_cor_val > 0.7]
        id <- (1 - sum(!grepl(true_values[[row,"PNID"]], categories))/length(categories))*100
        true_values[row,"PNID Accuracy"] <- id
        print(id)  
    }
    
    #if(!is.na(true_values[[row,"ID"]]) & true_values[[row,"ID"]] != ""){
    #    print("Nonplastic Contamination % of Total Spiked")
    #    categories <- test_accuracy$material_class
        #remove uncertain particles 
        #categories <- categories[test_accuracy$max_cor_val > 0.7]
    #    np_id <- (sum(!grepl("poly", categories))/length(categories))*100
    #    true_values[row,"NP ID Percent"] <- np_id
    #    print(np_id)  
    #}
    
    if(!is.na(true_values[[row,"median_area"]]) & true_values[[row,"median_area"]] != ""){
        print("Accurate mean size % of mean")
        area <- median(test_accuracy$area_um2)/true_values[row,"median_area"] * 100
        true_values[row,"Area Accuracy"] <- area
        print(area)
    }
    
    if(!is.na(true_values[[row,"median_feret"]]) & true_values[[row,"median_feret"]] != ""){
        print("Accurate mean size % of mean")
        area <- median(test_accuracy$max_length_um)/true_values[row,"median_feret"] * 100
        true_values[row,"Size Accuracy"] <- area
        print(area)
    }
}

true_values$`Count Accuracy` |> mean(na.rm = T)
true_values$`Count Accuracy` |> sd(na.rm = T)/true_values$`Count Accuracy` |> mean(na.rm = T)

true_values$`SpecID Accuracy` |> mean(na.rm = T)
true_values$`SpecID Accuracy` |> sd(na.rm = T)/true_values$`SpecID Accuracy` |> mean(na.rm = T)

true_values$`PNID Accuracy` |> mean(na.rm = T)
true_values$`PNID Accuracy` |> sd(na.rm = T)/true_values$`PNID Accuracy` |> mean(na.rm = T)

true_values$`Area Accuracy` |> mean(na.rm = T)
true_values$`Area Accuracy` |> sd(na.rm = T)/true_values$`Area Accuracy` |> mean(na.rm = T)

true_values$`Size Accuracy` |> mean(na.rm = T)
true_values$`Size Accuracy` |> sd(na.rm = T)/true_values$`Size Accuracy` |> mean(na.rm = T)

##HumanHair_15Nov223_control, 214 particles. 
##Loofah Sponge #2_15Nov2023_control, 67 particles.
##Loofah Sponge_15Nov223_control, 96 particles
##PMMA_15Nov223_control, 109
##RedPETFibers_15Nov223_control, 44

#Check results single map ----
map <- read_envi("C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls/RedBrick_21Nov2023_control.dat", spectral_smooth = T)
map <- read_any("C:\\Users\\winco\\OneDrive\\Documents\\CSULB MP Analysis\\100-300_A_LT.rds")
map$metadata$sig <- sig_noise(map |> restrict_range(min = c(750, 2700),
                                                    max = c(1800, 4000),
                                                    make_rel = F), 
                              metric = "sig_times_noise", abs = F)

heatmap_spec(map, sn = map$metadata$sig)

plot(sample_spec(map, 1))


#RedPETFibers_15Nov223_control
#BlueSiliconeGasket_21Nov2023_control
#Recovery_red_beads_75-90um_5um-screen
baseline <- median(as.matrix(map$spectra))

map_proc <- map |> restrict_range(min = c(750, 2700),
                   max = c(1800, 4000),
                   make_rel = F) 

map_proc$spectra <- map_proc$spectra - baseline

median(as.matrix(map_proc$spectra))
#map_proc <- process_spec(map, 
#                         restrict_range = T, 
#                         smooth_intens = F,
#                         #smooth_intens_args = list(polynomial = 3, window = 11, derivative = 1, abs = TRUE),
#                         restrict_range_args = list(min = 750, max = 4000), 
#                         flatten_range = T, #Improved ids by 7 %
#                         make_rel = F)

distances <- vapply(map_proc$spectra, function(x){
    sum(dist(x)/dist(map_proc$wavenumber)) 
}, numeric(1))

maxs <- vapply(map_proc$spectra, function(x){
    max(x) 
}, numeric(1))

mins <- vapply(map_proc$spectra, function(x){
    min(x) 
}, numeric(1))

rss <- vapply(map_proc$spectra, function(x){
  sqrt(sum(x^2))  
}, numeric(1))

rms <- vapply(map_proc$spectra, function(x){
    sqrt(mean(x^2))  
}, numeric(1))

ss <- vapply(map_proc$spectra, function(x){
    sum(x^2)  
}, numeric(1))

mean_t <- vapply(map_proc$spectra, function(x){
    mean(x)  
}, numeric(1))

sd_t <- vapply(map_proc$spectra, function(x){
    sd(x) 
}, numeric(1))

ggplot() +
    geom_bin2d(bins = 500, aes(x = mean_t, y = sd_t)) +
    scale_fill_continuous(type = "viridis") +
    theme_bw() +
    coord_fixed()

test <- sd_t * mean_t > 0.01

ggplot() +
    geom_point(aes(x = mean_t, y = sd_t, color = test)) +
    scale_color_discrete(type = "viridis") +
    theme_bw() +
    coord_fixed()

test_val =  156315
rss[test_val]
sd_t[test_val]
sum_t[test_val]

plot(filter_spec(map_proc, test_val), make_rel = F)

matches <- match_spec(process_spec(filter_spec(map, 170835), 
                                   conform_spec = T,
                                   conform_spec_args = list(range = lib$wavenumber, res = NULL), 
                                   restrict_range = T, 
                                   restrict_range_args = list(min = 750, max = 4000), 
                                   flatten_range = T #Improved ids by 7 %
        ),
lib, add_library_metadata = "sample_name")

plot(filter_spec(lib, "e04d65c1d75a75f51da809bc036e9e3b"))
 
plotly_spec(filter_spec(lib, "e04d65c1d75a75f51da809bc036e9e3b"), filter_spec(map, 170835) |> process_spec())

#Merge Maps ----



#Read in the file
    

#Image overlay ----
#library(raster)
library(imager)
library(magick)

t1 <- "C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\Mosaic Image files\\Silky Terrier Hair_21Nov2023_control_b.JPG"
t2 <- "C:\\Users\\winco\\OneDrive\\Documents\\QC3\\mosaic images\\MIPPR_IDC_LFB_22Feb24_GL_33_lft.JPG"
t3 <- "C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\Mosaic Image files\\Recovery_red_beads_75-90um_5um-screen.JPG"

mosaic <- jpeg::readJPEG(t1)

# Identify red pixels based on the specified condition
red_pixels <- mosaic[,,1]*255 > 50 & mosaic[,,1] > 2 * mosaic[,,2] & mosaic[,,1] > 2 * mosaic[,,3]

# Create a black and white image for visualization
bw_image <- array(0, dim = c(dim(red_pixels), 3))
bw_image[,,1] <- bw_image[,,2] <- bw_image[,,3] <- 1 # Set all pixels to white
bw_image[red_pixels] <- 0 # Set red pixels to black

# Plot the black and white image
plot(as.raster(bw_image))

# Get the x and y coordinates of the red pixels
red_coords <- which(red_pixels, arr.ind = TRUE)
x_coords <- red_coords[,2]
y_coords <- red_coords[,1]

# Calculate histograms of the x and y coordinates
x_hist <- sort(table(x_coords), decreasing = T)
y_hist <- sort(table(y_coords), decreasing = T)

# Find the two x and y values with the highest counts
x_boundaries <- as.integer(names(x_hist[c(TRUE, abs(diff(as.numeric(names(x_hist)))) > 10)][1:2]))
y_boundaries <- as.integer(names(y_hist[c(TRUE, abs(diff(as.numeric(names(y_hist)))) > 10)][1:2]))

# Sort the boundaries
x_min <- min(x_boundaries)
x_max <- max(x_boundaries)
y_min <- max(y_boundaries)
y_max <- min(y_boundaries)

bl = c(x_min, y_min)
tr = c(x_max, y_max)

plot(as.raster(mosaic))
rect(xleft = x_min, ybottom = y_min, xright = x_max, ytop = y_max, border = "black", lwd = 1)

map <- read_envi("C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\Recovery_red_beads_75-90um_5um-screen.dat")

snr <- sig_noise(map |> restrict_range(min = c(750, 2700),
                                       max = c(1800, 4000),
                                       make_rel = F) |> flatten_range(make_rel = F), 
                 metric = "sig_times_noise")

map_dim <- c(length(unique(map$metadata$x)), 
             length(unique(map$metadata$y)))
top_right = c(x_max, y_max)
bottom_left = c(x_min, y_min)
xscale = (top_right[1]-bottom_left[1])/map_dim[1]
yscale = (bottom_left[2]-top_right[2])/map_dim[2]
x_vals = as.integer(map$metadata$x*xscale+bottom_left[1])
y_vals = as.integer(bottom_left[2] - map$metadata$y*yscale)
particles  = snr >= 0.04

# Convert the image to a raster
image_raster <- as.raster(mosaic2)

# Create a matrix of coordinates for indexing
coords <- cbind(y_vals[particles], x_vals[particles])

image_raster[coords] <- "#00FF00"


png('file.png', height=nrow(image_raster), width=ncol(image_raster)) 
plot(image_raster)
dev.off()

# Convert the raster to an array for saving
raster_array <- col2rgb(image_raster) / 255
raster_array <- array(raster_array, dim = c(dim(image_raster)[1], dim(image_raster)[2], 3))

raster_array <- col2rgb(image_raster)
raster_array <- array(raster_array, dim = c(3, nrow(image_raster), ncol(image_raster)))

# Transpose the array to match the format expected by writePNG
raster_array <- aperm(raster_array, c(2, 3, 1)) / 255

# Save the updated image as a PNG file
save_path <- "updated_image.png"
png::writePNG(raster_array, save_path)

jpeg::writeJPEG(raster_array, 'test.jpg')

plot(image_raster)

# Print the boundaries
cat("x_min:", x_min, "x_max:", x_max, "y_min:", y_min, "y_max:", y_max, "\n")

# Visualize the boundaries on the original image
plot(as.raster(mosaic2))
rect(xleft = x_min, ybottom = y_min, xright = x_max, ytop = y_max, border = "black", lwd = 1)

image_raster2 <- as.raster(mosaic2)

red_colors <- c("#CC343C", 
                "#A13938")

image_b <- image_read(t1)

image_b_imager <- as.cimg(t1)

image_b_hsv <- RGBtoHSV(image_b_imager)


image_b_hsv <- image_convert(image_b, colorspace = "HSV")



red_channel <- channel(image_b_hsv, "red")

red_threshold <- red_channel > 0.5


image_b_hsv <- image_convert(image_b, colorspace = "HSV")

mosaic <- image_read(t3)
image_raster <- as.raster(mosaic)
mosaic2 <- jpeg::readJPEG(t3)
image_raster2 <- as.raster(mosaic2)

hex <- 
image_matrix <- matrix(data = as.character(image_raster), nrow = nrow(image_raster), byrow = TRUE)
check <- c(image_matrix)

mosaic <- image_draw(mosaic)
height <- image_info(mosaic)$height
#offset 42 y, 170 x t1, image dim = 766x515, map dim = 466x466
#offset 39 y, 288 x t2, image dim = 771x538, map dim = 301x640
#t3, bottom left 179,480, top right 624,21. map = 196x202
bottom_left = c(179, 480)
top_right = c(624, 21)
map_dim = c(196,202)
xscale = (top_right[1]-bottom_left[1])/map_dim[1]
yscale = (bottom_left[2]-top_right[2])/map_dim[2]
particle_centroid = c(875, 4675)/25
x = particle_centroid[1]*xscale+bottom_left[1]
y = bottom_left[2] - particle_centroid[2]*yscale

points(x = x, y = y, col = "black", pch = 5, cex = 0.5)  # x and y are the coordinates for the dot
dev.off()
image_raster <- as.raster(mosaic)
pixel_value = image_raster[height-235, 420]

dt <- data.frame(offset = c(170, 288), img = c(766, 771), map = c(466, 301))

model <- lm(offset ~ img + map, data = dt)
summary(model)

detection <- image_read("C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\particle_heatmap_thresholdedPMMA_15Nov223_control.jpg")

mosaic2 <- image_scale(mosaic, geometry = format("%dx%d!", image_info(mosaic)$width, image_info(mosaic)$height))

image_composite(mosaic, 
                image_scale(
                    image_crop(detection, 
                               paste0(image_info(detection)$width - 100, 
                                      "x", 
                                      image_info(detection)$height - 100, 
                                      "+50+50")
                               ),
                    paste0(image_info(mosaic)$width-280, "x", image_info(mosaic)$width-280)), 
                operator = "Overlay", 
                offset = geometry_point(x = 157, y = 0))

image_crop(detection, "900X900+50+50")

dim(m)
dim(detection)
resized_detection <- resize(detection, size_x = 100, size_y = 100)
test <- implot(mosaic, plot(detection * 0.5))

plot(test)

# Crop the image
cropped_image <- crop(mosaic, x1 = 100, x2 = 300, y1 = 50, y2 = 200)

# Resize image2 to match the dimensions of image1
image2_resized <- imresize(mosaic, dim(detection)[1:2])

# Blend the images
# alpha controls the transparency of the second image
overlay_image <- detection * 0.5 + mosaic * 0.5

plot(overlay_image)

mosaic <- brick("C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\Mosaic Image files\\PMMA_15Nov223_control.JPG", values = T)
detection <- brick("C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\particle_heatmap_thresholdedPMMA_15Nov223_control.jpg")

plotRGB(mosaic)
plotRGB(mosaic_resampled)
plotRGB(detection, add = TRUE, alpha = 0.5)

mosaic_resampled <- resample(mosaic, detection, method="bilinear")

image_raster2 <- projectRaster(from=mosaic, to=detection)

plot(image_raster2)

detection[,,4] = 0.5  # set alpha to semi-transparent
mosaic[,,4] = 1  # set alpha to 1 (opaque)

png('test.png', width = 2, height = 2, units = 'in', res = 150)
par(mai=c(0,0,0,0))
plot.new()
rasterImage(detection, 0, 0, 1, 1)
rasterImage(mosaic,   0, 0, 1, 1)
dev.off()


#Library reformatting ----
cors <- cor_spec(lib, lib)
diag(cors) <- NA
change = T
while(change){
    print(nrow(cors))
    find_top <- apply(cors, 2, which.max)
    true_identities <- lib$metadata$material_class[match(colnames(cors), lib$metadata$sample_name)]
    print(sort(table(true_identities), decreasing = T)/length(true_identities))
    top_identities <- lib$metadata$material_class[match(rownames(cors)[find_top], lib$metadata$sample_name)]
    accurate = top_identities == true_identities
    cors2 <- cors[accurate, accurate]
    change = nrow(cors) != nrow(cors2)
    cors = cors2     
}

#This reduction seems to be working similarly to raw library but below is decreasing accuracy by about 5%
lib_filtered <- filter_spec(lib, colnames(cors))

wavenumbers <- lib_filtered$wavenumber

#Check this, make sure its what you want to test. 
wavenumber_filter <- (wavenumbers >= 800 & wavenumbers <= 3200) & !(wavenumbers >= 2200 & wavenumbers <= 2420)#(wavenumbers >= 1000 & wavenumbers <= 2000)#
# Removing CO2 region entirely.

new_wavenumbers = wavenumbers[wavenumber_filter]

nobaseline <- transpose(lib_filtered$spectra[wavenumber_filter,][, lapply(.SD, function(x){make_rel(x, na.rm = T)})][,wavenumbers := new_wavenumbers], keep.names = "id", make.names = "wavenumbers") %>%
    na.omit()  
#normalizes everything 0-1

ids <- lib_filtered$metadata %>%
    filter(sample_name %in% nobaseline$id) %>%
    dplyr::select(sample_name, spectrum_type, material_class) %>%
    mutate(material_class = gsub("α", "a", tolower(material_class))) %>%
    mutate(material_class = ifelse(material_class %in% c("polyacrylamides", "polyamides (polylactams)"), "poly(acrylamide/amid)s", material_class)) %>%
    mutate(material_class = ifelse(material_class %in% c("polyterephthalates", "polyesters", "polyethers (polyglycols, oxyalkylenes, glycidyl ethers & cyclic ethers)", "polydiglycidyl ethers (polyepoxides, polyhydroxyethers, phenoxy)"), "poly(esters/ethers/diglycidylethers/terephthalates)s", material_class)) %>%
    #filter(!material_class %in% c("other material", "other plastic")) %>%
    mutate(material_class = paste0(tolower(spectrum_type), "_", tolower(material_class))) %>%
    mutate(sample_level = as.numeric(as.factor(material_class))) %>%
    select(-spectrum_type)

k = 100
count <- table(ids$material_class)
#count2 <- ifelse(count > k, k, ceiling(count/2))

to_reduce <- count[count > k]
#How many of the classes are getting reduced?
length(to_reduce)/length(unique(ids$material_class))

ids_to_keep <- lapply(1:length(to_reduce), function(x){
    print(names(to_reduce)[x])
    
    id_to_reduce <- ids %>%
        dplyr::filter(material_class == names(to_reduce)[x])
    
    spectra_to_reduce2 <- transpose(nobaseline %>%
                                        filter(id %in% id_to_reduce$sample_name),
                                    make.names = "id", 
                                    keep.names = "wavenumber") %>%
        select(-wavenumber) 
    
    spectra_not_to_reduce <- transpose(nobaseline %>%
                                           filter(!id %in% id_to_reduce$sample_name),
                                       make.names = "id", 
                                       keep.names = "wavenumber") %>%
        select(-wavenumber)
    
    internal_cors <- cor(spectra_to_reduce2)
    
    distance_cors <- as.dist(1-internal_cors)
    
    pam.res <- pam(distance_cors, k = k, diss = TRUE, pamonce = 6)
    
    names(spectra_to_reduce2)[pam.res$id.med]
})


not_to_reduce <- ids %>%
    filter(!material_class %in% names(to_reduce)) %>%
    pull(sample_name)

ids_to_keep_with_small <- c(unlist(ids_to_keep), not_to_reduce)

compressed_library <- filter_spec(lib_filtered, ids_to_keep_with_small)

#Use Blanks For Correction ----
blanks <- list.files(path = "C:\\Users\\winco\\OneDrive\\Documents\\MB2\\Export Files",pattern =  "PrcdrlBlnk.*\\.rds$", full.names = T)
other <- list.files(path = "C:\\Users\\winco\\OneDrive\\Documents\\MB2\\Export Files",pattern =  "\\.rds$", full.names = T)
other <- other[!other %in% blanks]

blank_spec <- read_any(blanks) |> c_spec(res = NULL)
other_spec <- read_any(other) |> c_spec(res = NULL)
lib <- load_lib("derivative") %>% filter_spec(., 
                                                 !.$metadata$material_class %in% c("other plastic", 
                                                                                     "other material"#, 
                                                                                     #"mineral", 
                                                                                     #"organic matter"
                                                 ) & 
                                                     .$metadata$spectrum_type == "ftir")

lib_no_blanks <- filter_spec(lib, !lib$metadata$sample_name %in% blank_spec$metadata$max_cor_name)

names(other_spec$metadata) <- paste0("sample_", names(other_spec$metadata))

is_OpenSpecy(other_spec)

matches <- match_spec(other_spec, blank_spec,top_n = 1, add_library_metadata = "col_id", add_object_metadata = "sample_col_id")
matches2 <- match_spec(other_spec, lib_no_blanks,top_n = 1, add_library_metadata = "sample_name", add_object_metadata = "sample_col_id")

matches_analyzed <- matches %>%
    mutate(more_similar = match_val > sample_max_cor_val,
           similarity_ratio = match_val/sample_max_cor_val)

sum(matches_analyzed$more_similar)/nrow(matches_analyzed)

hist(matches_analyzed$similarity_ratio)

#test results
summary_matches_precent_all <- matches_analyzed %>%
    mutate(more_similar = as.character(more_similar)) %>%
    group_by(more_similar, sample_material_class)%>%
    summarise(count = n()) %>%
    ungroup() %>%
    group_by(more_similar) %>%
    dplyr::mutate(percent = count/sum(count) * 100) 

all <- summary_matches_precent_all %>%
    group_by(sample_material_class) %>%
    summarise(count = sum(count)) %>%
    ungroup() %>%
    mutate(percent = count/sum(count) * 100) %>%
    mutate(more_similar = "All") %>% 
    as.data.frame()

blank_summary <- blank_spec$metadata %>%
    group_by(material_class) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    dplyr::mutate(percent = count/sum(count) * 100) %>%
    mutate(more_similar = "Blanks") %>%
    rename(sample_material_class = material_class)

noblanks_summary <- matches2 %>%
    mutate(material_class = ifelse(match_val > 0.66, material_class, "unknown")) %>%
    dplyr::filter(!material_class %in% c("organic matter",  "mineral")) %>%
    group_by(material_class) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    dplyr::mutate(percent = count/sum(count) * 100) %>%
    mutate(more_similar = "NoBlanks") %>%
    rename(sample_material_class = material_class)

lib_summary <- lib$metadata %>%
    dplyr::filter(!material_class %in% c("organic matter",  "mineral")) %>%
    group_by(material_class) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    dplyr::mutate(percent = count/sum(count) * 100) %>%
    mutate(more_similar = "Lib") %>%
    rename(sample_material_class = material_class)

study_results_percent_all <- bind_rows(summary_matches_precent_all, 
                                       all,
                                       blank_summary,
                                       noblanks_summary,
                                       lib_summary)

ggplot(study_results_percent_all, aes(x = percent, y = sample_material_class, fill = sample_material_class)) +
    geom_bar(stat = "identity") +
    facet_wrap(more_similar ~.) +
    theme_bw(base_size = 8) +
    geom_text(aes(label = count), size = 2, hjust = 0) +
    scale_fill_viridis_d() +
    theme(legend.position = "none")

(table(lib_no_blanks$metadata$material_class)/length(lib_no_blanks$metadata$material_class)) |> 
    sort()


#Test individual overlap
material = "polyolefins (polyalkenes)"

other_plastic <- filter_spec(other_spec, other_spec$metadata$sample_material_class == material)

blanks_plastic <- filter_spec(blank_spec, blank_spec$metadata$material_class == material)

lib_blanks <- filter_spec(lib, blanks_plastic$metadata$max_cor_name)
lib_other <- filter_spec(lib, other_plastic$metadata$sample_max_cor_name)

table(blanks_plastic$metadata$max_cor_name)
plot(blanks_plastic)
plot(lib_blanks)
plotly_spec(sample_spec(blanks_plastic, 5), sample_spec(lib_blanks, 5))

lib_blanks$metadata$organization |> table()
lib_other$metadata$organization |> table()




