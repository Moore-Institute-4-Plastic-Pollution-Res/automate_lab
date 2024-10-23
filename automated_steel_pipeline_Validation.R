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

#Run once
# check if it has been installed or not 
get_lib(c("mediod", "derivative"))

#------------------- Google Drive File Download ------------------------------
# Authorize Google Drive Connection 
drive_auth()

# Access ID of all files in each folder
files <- c("spectra", "sample data", "images")

for (file in files) {}
spectra <- shared_drive_find("test_data") |> 
    drive_ls("spectra") |> 
    drive_ls() 

sample_data <- shared_drive_find("test_data") |>
    drive_ls("sample data") |> 
    drive_ls()

images <- shared_drive_find("test_data") |> 
    drive_ls("images") |> 
    drive_ls()






# Create folders to store files downloaded
if (!dir.exists("spectra_data_download")) {
    dir.create("spectra_data_download")
} else if (!dir.exists("sample_data_download")) {
    dir.create("sample_data_download")
} else if (!dir.exists("images_download")){
    dir.create("images_download")
}

# Download data
for (i in 1:nrow(spectra)){
    
    # File name
    file_name <- spectra$name[i]
    # Set download path
    dest_path <- file.path("spectra_data_download", file_name)
    # Download
    drive_download(spectra[i,], path = dest_path, overwrite = TRUE)
}


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


analyze_features <- function(files, 
                             img = NULL, 
                             bottom_left = NULL, 
                             top_right = NULL, 
                             lib,
                             adj_map_baseline = F,
                             material_class = "material_class",
                             origins = list(x = rep(0, 
                                                    length(files)), 
                                            y = rep(0, length(files))),
                             spectral_smooth = FALSE, 
                             sigma1 = c(1,1,1), 
                             spatial_smooth = FALSE, 
                             sigma2 = c(3,3),
                             close = F,
                             close_kernel = c(4,4),
                             sn_threshold = 0.04,
                             cor_threshold = 0.7, 
                             area_threshold = 1,
                             label_unknown = F,
                             remove_nonplastic = F,
                             remove_unknown = F,
                             pixel_length = 25, 
                             metric = "sig_times_noise",
                             abs = F, 
                             k = 1, 
                             k_weighting = "mean",
                             wd, 
                             particle_id_strategy = "collapse",
                             vote_count = Inf, 
                             collapse_function = median,
                             types = c("particle_image",
                                       "particle_heatmap",
                                       "particle_heatmap_thresholded",
                                       "median_spec_plot",
                                       "median_spec",
                                       "particle_details",
                                       "particle_summary",
                                       "spectra_raw", 
                                       "spectra_processed",
                                       "time"), 
                             by = c("material", "sample", "all"),
                             width = 850, height = 850, units = "px"){
    #Variable input tests. 
    if(!all(particle_id_strategy %in% c("collapse", "all cell id", "particle cell vote"))) stop("Incorrect particle_id_strategy")
    if(!all(is.character(files))) stop("Incorrect files")
    if(!all(is.character(img)|is.null(img))) stop("Incorrect img")
    if(!all(is.list(bottom_left)|is.null(bottom_left))) stop("Incorrect bottom_left")
    if(!all(is.list(top_right)|is.null(top_right))) stop("Incorrect top_right")
    if(!all(is_OpenSpecy(lib)|is.list(lib))) stop("Incorrect lib")
    if(!all(is.logical(adj_map_baseline))) stop("Incorrect adj_map_baseline")
    if(!all(is.character(material_class))) stop("Incorrect material_class")
    if(!all(is.list(origins))) stop("Incorrect origins")
    if(!all(is.logical(spectral_smooth))) stop("Incorrect spectral_smooth")
    if(!all(is.numeric(sigma1))) stop("Incorrect sigma1")
    if(!all(is.logical(spatial_smooth))) stop("Incorrect spatial_smooth")
    if(!all(is.numeric(sigma2))) stop("Incorrect sigma2")
    if(!all(is.logical(close))) stop("Incorrect close")
    if(!all(is.numeric(close_kernel))) stop("Incorrect close_kernel")
    if(!all(is.numeric(sn_threshold))) stop("Incorrect sn_threshold")
    if(!all(is.numeric(cor_threshold))) stop("Incorrect cor_threshold")
    if(!all(is.numeric(area_threshold))) stop("Incorrect area_threshold")
    if(!all(is.logical(label_unknown))) stop("Incorrect label_unknown")
    if(!all(is.logical(remove_nonplastic))) stop("Incorrect remove_nonplastic")
    if(!all(is.logical(remove_unknown))) stop("Incorrect remove_nonplastic")
    if(!all(is.numeric(pixel_length))) stop("Incorrect pixel_length")
    if(!all(is.character(metric))) stop("Incorrect metric")
    if(!all(is.logical(abs))) stop("Incorrect abs")
    if(!all(is.numeric(k))) stop("Incorrect k")
    if(!all(is.character(k_weighting))) stop("Incorrect k_weighting")
    if(!all(is.character(wd))) stop("Incorrect wd")
    if(!all(is.character(particle_id_strategy))) stop("Incorrect particle_id_strategy")
    if(!all(is.character(types))) stop("Incorrect types")
    if(!all(is.character(by))) stop("Incorrect by")
    if(!all(is.numeric(width))) stop("Incorrect width")
    if(!all(is.numeric(height))) stop("Incorrect height")
    if(!all(is.character(units))) stop("Incorrect units")

    #Extracts the file names
    file_names <- gsub("(.*/)|(\\..{1,3}$)", "", files)
    
    #Corrects names to not throw error with some functions
    lib$metadata[[material_class]] <- gsub("/", "_", lib$metadata[[material_class]])
    
    #Loops through all files
    for(file in 1:length(files)){
        #Prints file names
        print(files[file])
        
        #Begins timer
        time_start = Sys.time()
        
        #Read in the file
        if(grepl(".dat", files[file])){
            map <- read_envi(files[file], spectral_smooth = spectral_smooth, sigma = sigma1) |>
                       adj_intens(type = "transmittance", make_rel = F) #If converting from transmittance, otherwise hash out
        }
        
        else{
            map <- read_any(files[file]) 
        }
        
        #Correct baseline if specified
        if(adj_map_baseline){

            baseline <- median(as.matrix(map$spectra), na.rm = T)
            print(baseline)
            map$spectra <- map$spectra - baseline
            
        }

        #extracts origins if available. 
        if(!is.null(origins)){
            originx <- origins[[1]][file]
            originy <- origins[[2]][file]            
        }
        
        #Auto attempts to extract origins if description available. 
        else if("description" %in% names(map$metadata)){
            origin = unique(map$metadata$description)
            originx = gsub(",.*", "", gsub(".*X=", "",  origin)) |> as.numeric()
            originy = gsub(".*Y=", "",  origin) |> as.numeric()
        }
        
        #Set origins to zero if neither of the above work
        else{
            originx <- list(x = rep(0, length(files)), y = rep(0, length(files)))[[1]][file]
            originy <- list(x = rep(0, length(files)), y = rep(0, length(files)))[[2]][file] 
        }

        #Estimate the sig noise value on the raw data
        map$metadata$snr <- sig_noise(map |>
                                          restrict_range(min = c(800, 2420),
                                                         max = c(2200, 3200),
                                                         make_rel = F), #Making comparable to lda and mediod
                                    metric = metric,
                                    spatial_smooth = spatial_smooth,
                                    sigma = sigma2,
                                    abs = abs)
        
        #Create a heatmap of  the sig noise values
        if("particle_heatmap" %in% types){
            plot <- ggplot() +
                geom_raster(data = map$metadata, 
                            aes(x = x * pixel_length + originx, 
                                y = y * pixel_length + originy, 
                                fill = snr)) +
                coord_fixed() +
                theme_bw() +
                labs(x = "X (um)", y = "Y (um)") +
                scale_fill_viridis_c()
            
            ggsave(filename = paste0("particle_heatmap_", 
                                     file_names[file], 
                                     ".png"), 
                   plot = plot, 
                   path = wd, 
                   width = width, 
                   height = height, 
                   units = units)
            
            if(all(types == "particle_heatmap")) next
        }
        
        #Create a thresholded heatmap using the sig noise values and threshold
        if("particle_heatmap_thresholded" %in% types){
            plot <- ggplot() +
                geom_raster(data = map$metadata, 
                            aes(x = x * pixel_length + originx, 
                                y = y * pixel_length + originy, 
                                fill = snr > sn_threshold)) +
                scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +  # Assign black to TRUE and white to FALSE
                coord_fixed() +
                theme_bw() +
                theme(legend.position = "none") +
                labs(x = "X (um)", y = "Y (um)")
            
            ggsave(filename = paste0("particle_heatmap_thresholded", 
                                     file_names[file], 
                                     ".jpg"), 
                   plot = plot, 
                   path = wd, 
                   width = width, 
                   height = height, 
                   units = units)
            
            #Skip to the next map if only requesting thresholded outputs.
            if(all(types %in% c("particle_heatmap", "particle_heatmap_thresholded"))) next
        }
        
        #Check if all or no features are thresholded and skip to the next if so. 
        if(sum(map$metadata$snr > sn_threshold) == 0 | all(map$metadata$snr > sn_threshold)) next
    
        #Auto detect the in the images using the red boxes. 
        if(is.null(bottom_left[[file]]) & !is.null(img[file])){
            mosaic <- jpeg::readJPEG(img[file])
            
            # Identify red pixels based on the specified condition
            red_pixels <- mosaic[,,1]*255 > 50 & mosaic[,,1] > 2 * mosaic[,,2] & mosaic[,,1] > 2 * mosaic[,,3]
            
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
            
        }
         else{
             bl = bottom_left[[file]]
             tr = top_right[[file]]
         }
        
        
        #Identify every pixel and use that to identify the features. 
        if(particle_id_strategy == "all cell id"){

                proc_map  <- process_spec(filter_spec(map, map$metadata$snr > sn_threshold), 
                                 conform_spec = T,
                                 conform_spec_args = list(range = if(!is_OpenSpecy(lib)){lib$all_variables} else {lib$wavenumber}, res = NULL), 
                                 restrict_range = T, 
                                 restrict_range_args = list(min = c(800, 2420), max = c(2200, 3200))#, #Making comparable to lda and mediod), 
                                 #flatten_range = T #Improved ids by 7 %
                ) 
            
        }
        
        #Just use the particle thresholds for identifying where the particles are. 
        else{
            map$metadata$threshold <- map$metadata$snr > sn_threshold
            
            map <- def_features(map, 
                                features = map$metadata$threshold,
                                shape_kernel = sigma2,
                                shape_type = "box",
                                close = close, 
                                close_kernel = close_kernel,
                                close_type = "box",
                                img[file], 
                                bl,
                                tr)
        }
        
       
        #Collapse the particles to their median. Remove -88 and areas with only 1 pixel. Process the spectra. 
        if(particle_id_strategy == "collapse"){
            if(sum(map$metadata$feature_id != "-88" & map$metadata$area > area_threshold) == 0) next
            proc_map  <- collapse_spec(map, fun = collapse_function) %>%
                filter_spec(., .$metadata$feature_id != "-88" & .$metadata$area > area_threshold) %>%
                process_spec(., 
                             conform_spec = T,
                             conform_spec_args = list(range = if(!is_OpenSpecy(lib)){lib$all_variables} else {lib$wavenumber}, res = NULL), 
                             restrict_range = T, 
                             restrict_range_args = list(min = c(800, 2420), max = c(2200, 3200))#)#, #Making comparable to lda and mediod), 
                             #flatten_range = T #Improved ids by 7 %
                )    
        }
        
        #Process data for particle cell voting, don't compress. 
        if(particle_id_strategy == "particle cell vote"){
            if(sum(map$metadata$feature_id != "-88" & map$metadata$area > area_threshold) == 0) next
            
            proc_map  <- map %>%
                filter_spec(., .$metadata$feature_id != "-88" & .$metadata$area > area_threshold) %>%
                process_spec(., 
                             conform_spec = T,
                             conform_spec_args = list(range = if(!is_OpenSpecy(lib)){lib$all_variables} else {lib$wavenumber}, res = NULL), 
                             restrict_range = T, 
                             restrict_range_args = list(min = c(800, 2420),
                                                        max = c(2200, 3200))#, #Making comparable to lda and mediod),
                             #flatten_range = T #Improved ids by 7 %
                ) 
        }
    
        #Run AI based classification
        if(!is_OpenSpecy(lib)){
            cors <- match_spec(proc_map, lib)
            proc_map$metadata$max_cor_val <- cors$value
            proc_map$metadata$material_class <- cors$name
        }
    
    #Run multi-k classification
    else if(k > 1){
        if(k_weighting == "raw"){
            cors <- match_spec(proc_map, 
                               lib, 
                               top_n = k,
                               add_library_metadata = "sample_name")  %>%
                group_by(object_id, material_class) %>%
                summarise(prob = n()/k) %>%
                group_by(object_id) %>%
                top_n(prob, n = 1) %>%
                ungroup() %>%
                rename(max_cor_val = prob)    
        }
        
        if(k_weighting == "sum"){
            cors <- match_spec(proc_map, 
                               lib, 
                               top_n = k,
                               add_library_metadata = "sample_name")  %>%
                group_by(object_id, material_class) %>%
                summarise(prob = (n()/k+mean(match_val))/2) %>%
                group_by(object_id) %>%
                top_n(prob, n = 1) %>%
                ungroup() %>%
                rename(max_cor_val = prob)    
        }
        
        if(k_weighting == "multiple") {
            cors <- match_spec(proc_map, 
                               lib, 
                               top_n = k,
                               add_library_metadata = "sample_name")  %>%
                group_by(object_id, material_class) %>%
                summarise(prob = n()/k*mean(match_val)) %>%
                group_by(object_id) %>%
                top_n(prob, n = 1) %>%
                ungroup() %>%
                rename(max_cor_val = prob)
        }
        
            
        proc_map$metadata <- left_join(proc_map$metadata,
                                       cors, 
                                       by = c("feature_id" = "object_id"))
        
        
    }
        
    #Run standard top 1 k classification
    else{
        cors <- cor_spec(proc_map, lib)
        
        max_cors <- max_cor_named(cors)
        
        proc_map$metadata$max_cor_val <- max_cors
        
        proc_map$metadata$max_cor_name <- names(max_cors)
        
        proc_map$metadata <- left_join(proc_map$metadata, 
                                       lib$metadata %>% 
                                           dplyr::select(sample_name, material_class), 
                                       by = c("max_cor_name" = "sample_name"))
        
        
    }
        #If using all cell id set the thresholds and identify the features.  
        if(particle_id_strategy == "all cell id"){

            map$metadata <- left_join(map$metadata, 
                                      dplyr::select(proc_map$metadata, 
                                                    x, 
                                                    y, 
                                                    max_cor_val, 
                                                    material_class),
                                      by = c("x", "y")) %>%
                mutate(threshold = snr > sn_threshold & max_cor_val > cor_threshold) %>%
                mutate(material_class = ifelse(threshold, material_class, "background"))
            
            map <- def_features(x = map, 
                                 features = map$metadata$material_class,
                                 close = T, 
                                 close_kernel = c(4,4), 
                                 close_shape = "box") 
            
            map$metadata <- map$metadata %>%
                mutate(material_class = ifelse(!grepl("background", feature_id) & area > area_threshold, gsub("_.*", "", feature_id), NA))
            
            proc_map <- filter_spec(map, !is.na(map$metadata$feature_id)) %>%
                collapse_spec(fun = collapse_function)
            
        }
        
        #If using particle cell voting conduct the voting process
        if(particle_id_strategy == "particle cell vote"){
            ts <- transpose(proc_map$spectra)
            ts$id <- proc_map$metadata$feature_id
            proc_map$spectra <- ts[, lapply(.SD, median, na.rm = T), by = "id"] |>
                transpose(make.names = "id")
            
            # Step 1: Subset proc_map$metadata and then select unique rows based on certain columns
            subset_metadata <- proc_map$metadata[, .(feature_id, area, perimeter, feret_max, feret_min, centroid_y, centroid_x, first_x, first_y, rand_x, rand_y)]
            unique_metadata <- unique(subset_metadata)
            
            # Step 2: Order by Proportion and select the top entry per feature_id
            proportions <- proportions[order(-max_cor_val), .SD[1:min(vote_count, .N)], by = .(feature_id)]
            
            proportions <- proc_map$metadata[, .(sum_cor = sum(max_cor_val)), by = c("feature_id", material_class)]
            proportions[, max_cor_val := sum_cor / sum(sum_cor), by = .(feature_id)]

            # Step 3: Order by Proportion and select the top entry per feature_id
            top_proportions <- proportions[order(-max_cor_val), .SD[1], by = .(feature_id)]
            
            # Step 4: Join with unique_metadata
            proc_map$metadata <- unique_metadata[top_proportions, on = .(feature_id)]
            
        }
        
        #Relabel spectra as unknown if below the threshold
        if(label_unknown){
            proc_map$metadata[[material_class]] <- ifelse(proc_map$metadata$max_cor_val < cor_threshold, "unknown",proc_map$metadata[[material_class]])
        }
        
        #Remove any nonplastic particles from further analysis
        if(remove_nonplastic){
            if(sum(!proc_map$metadata[[material_class]] %in% c("mineral", "organic matter", "other material")) == 0) next
            proc_map <- filter_spec(proc_map, !proc_map$metadata[[material_class]] %in% c("mineral", "organic matter", "other material"))
        }
        
        #Remove any nonplastic particles from further analysis
        if(remove_unknown){
            if(sum(!proc_map$metadata[[material_class]] %in% c("unknown")) == 0) next
            proc_map <- filter_spec(proc_map, !proc_map$metadata[[material_class]] %in% c("unknown"))
        }
        
        if(particle_id_strategy %in% c("particle cell vote", "collapse")){
            map$metadata <- left_join(map$metadata, dplyr::select(proc_map$metadata, feature_id, material_class), by = "feature_id")
        }
        
        #Plot overlay on mosaic
        if(!is.null(img[file])){
            map_dim <- c(length(unique(map$metadata$x)), 
                         length(unique(map$metadata$y)))
            xscale = (tr[1]-bl[1])/map_dim[1]
            yscale = (bl[2]-tr[2])/map_dim[2]
            x_vals = as.integer(map$metadata$x*xscale+bl[1])
            y_vals = as.integer(bl[2] - map$metadata$y*yscale)
            particles  = map$metadata$threshold
            
            # Convert the image to a raster
            mosaic <- jpeg::readJPEG(img[file])
            image_raster <- as.raster(mosaic)
            
            # Create a matrix of coordinates for indexing
            coords <- cbind(y_vals[particles], x_vals[particles])
            
            image_raster[coords] <- "#00FF00"
            
            png(paste0(wd, "/overlay_", file_names[file], ".png"), 
                height=nrow(image_raster), width=ncol(image_raster)) 
            
            plot(image_raster)
            dev.off()
        }
        
        #Plot sample level figures for particle types
        if("particle_image" %in% types & "sample" %in% by){
            plot <- particle_image(proc_map$metadata, map$metadata, material_class, pixel_length, origin = c(originx, originy),  title = "All Particles")
            ggsave(filename = paste0("particle_image_", file_names[file], ".png"), plot = plot, path = wd, width = width, height = height, units = units)
        }
        
        #Plot material-sample level figures for particle types
        if("particle_image" %in% types & "material" %in% by){
            for(material in unique(proc_map$metadata[[material_class]])){
                map_subset <- map$metadata %>%
                    mutate(!!material_class := ifelse(.data[[material_class]] == material, .data[[material_class]], NA))
                proc_map_subset <- proc_map$metadata %>%
                    dplyr::filter(.data[[material_class]] == material)
                plot <- particle_image(proc_map_subset, map_subset, material_class, pixel_length, origin = c(originx, originy), title = paste0(material, " Particles"))
                ggsave(filename = paste0(material, "_particle_image_", file_names[file], ".png"), plot = plot, path = wd, width = width, height = height, units = units)
            }
        }
        
        #Plot some example spectra for each material type in each sample
        if("median_spec_plot" %in% types & "sample" %in% by){
            metadata_to_use <- proc_map$metadata %>%
                group_by(material_class) %>%
                sample_n(10, replace = T) %>%
                distinct() %>%
                ungroup()
            
            spectra_sample <- proc_map$spectra %>%
                mutate(wavenumber = proc_map$wavenumber) %>%
                pivot_longer(cols = -wavenumber, names_to = "feature_id", values_to = "intensity") %>%
                inner_join(metadata_to_use %>% dplyr::select(feature_id, material_class), by = "feature_id") 
            
            spectra_plot <- ggplot(spectra_sample) +
                geom_line(aes(x = wavenumber, y = intensity, group = feature_id)) +
                facet_grid(rows = vars(.data[[material_class]])) +
                theme_bw(base_size = 5) +
                labs(title = paste0("Example spectra in ", file_names[file]), x = "Wavenumbers (1/cm)", y = "Absorbance (Min-Max Norm)") +
                scale_x_reverse() +
                theme(plot.title = element_text(size=5))
            
            ggsave(filename = paste0("median_spec_plot_", file_names[file], ".png"), plot = spectra_plot, path = wd, width = width, height = height, units = units)
        }
        
        #Create a unique figure for each material-sample spectra examples. 
        if("median_spec_plot" %in% types & "material" %in% by){
            for(material in unique(proc_map$metadata[[material_class]])){
                metadata_to_use <- proc_map$metadata %>%
                    filter(.data[[material_class]] == material) %>%
                    sample_n(10, replace = T) %>%
                    distinct() %>%
                    ungroup()
                
                spectra_sample <- proc_map$spectra %>%
                    mutate(wavenumber = proc_map$wavenumber) %>%
                    pivot_longer(cols = -wavenumber, names_to = "feature_id", values_to = "intensity") %>%
                    inner_join(metadata_to_use %>% dplyr::select(feature_id, material_class), by = "feature_id") 
                
                spectra_plot <- ggplot(spectra_sample) +
                    geom_line(aes(x = wavenumber, y = intensity, group = feature_id)) +
                    facet_grid(rows = vars(feature_id)) +
                    theme_bw(base_size = 5) +
                    labs(title = paste0("Example Spectra ", material, " in ", file_names[file]), x = "Wavenumbers (1/cm)", y = "Absorbance (Min-Max Norm)") +
                    scale_x_reverse() +
                    theme(plot.title = element_text(size=5))
                
                ggsave(filename = paste0(material, "_median_spec_plot_", file_names[file], ".png"), plot = spectra_plot, path = wd, width = width, height = height, units = units)
            }
        }
        
        #Create a dataset for the output. 
        if("particle_details" %in% types){
            particle_info <- proc_map$metadata %>%
                rename(particle_id = feature_id) %>%
                mutate(area = area * pixel_length^2, 
                       perimeter = perimeter * pixel_length, 
                       feret_min = feret_min * pixel_length, 
                       feret_max = feret_max * pixel_length,
                       centroid_x = centroid_x * pixel_length + originx, 
                       centroid_y = centroid_y * pixel_length + originy,
                       first_x = first_x * pixel_length + originx, #dont work with vote option yet. 
                       first_y = first_y * pixel_length + originy, 
                       #rand_x = rand_x * pixel_length + originx, 
                       #rand_y = rand_y * pixel_length + originy, 
                       bad_spectra = max_cor_val < cor_threshold,
                       acc_analy_conf = ifelse(max_cor_val > 0.6, "confident", ifelse(max_cor_val < 0.3, "undetermined", "possible")), 
                       sample_id = file_names[file]) %>%
                mutate(aspect_ratio = feret_max/feret_min, 
                       circularity = (perimeter^2)/(4*pi*area))
            
            selected_columns <- c("particle_id", "sample_id", "max_cor_val", 
                                  "bad_spectra", "material_class", "area", 
                                  "perimeter", "feret_max", "feret_min", 
                                  "aspect_ratio", "circularity", "centroid_x", 
                                  "centroid_y", "first_x", "first_y", "acc_analy_conf"#, 
                                  #"rand_x", "rand_y"#dont work with vote option yet. 
                                  )
            
            if ("max_cor_name" %in% names(particle_info)) {
                selected_columns <- c(selected_columns, "max_cor_name")
            }
            
            if ("mean_cor" %in% names(particle_info)) {
                selected_columns <- c(selected_columns, "mean_cor")
            }
            
            if ("mean_snr" %in% names(particle_info)) {
                selected_columns <- c(selected_columns, "mean_snr")
            }
            
            if (all(c("r", "g", "b") %in% names(particle_info))) {
                selected_columns <- c(selected_columns, c("r", "g", "b"))
            }
            
            particle_info <- particle_info %>%
                dplyr::select(all_of(selected_columns)) %>%
                rename(max_length_um = feret_max, 
                       min_length_um = feret_min, 
                       perimeter_um = perimeter,
                       area_um2 = area,
                       centroid_y = centroid_y, 
                       centroid_x = centroid_x)
            
            fwrite(particle_info, paste0(wd, "/particle_details_", file_names[file], ".csv"))
            
        }
        
        #Create summary tables for each sample
        if("particle_summary" %in% types){
            particle_summary <- proc_map$metadata %>%
                rename(particle_id = feature_id) %>%
                group_by(material_class) %>%
                summarise(count = n()) %>%
                mutate(sample_id = file_names[file])
            
            fwrite(particle_summary, paste0(wd, "/particle_summary_", file_names[file], ".csv"))
        }
        
        #Save the raw spectra
        if("spectra_raw" %in% types){
            saveRDS(map, file = paste0(wd, "/particles_raw_", file_names[file], ".rds"))
        }
        
        #Save the processed spectra
        if("spectra_processed" %in% types){
            saveRDS(proc_map, file = paste0(wd, "/particles_", file_names[file], ".rds"))
        }
        
        #Save the median spectra as a csv
        if("median_spec" %in% types){
            spectra_sample <- proc_map$spectra %>%
                mutate(wavenumber = proc_map$wavenumber) %>%
                pivot_longer(cols = -wavenumber, names_to = "feature_id", values_to = "intensity") %>%
                inner_join(metadata_to_use %>% dplyr::select(feature_id, material_class), by = "feature_id") %>%
                mutate(sample_id = file_names[file])  
            
            fwrite(spectra_sample, paste0(wd, "/median_spec_", file_names[file], ".csv"))
        }
        #Print the final run time. 
        time_diff = Sys.time() - time_start
        if("time" %in% types){
            saveRDS(time_diff, file = paste0(wd, "/time_", file_names[file], ".rds"))
        }
        print(time_diff)
    }
    
    #After all samples are finished, concatenate the summary, details, adn spectra csv tables. 
    if("particle_summary" %in% types & "all" %in% by){
        if(file.exists(paste0(wd, "/particle_summary_all.csv"))) file.remove(paste0(wd, "/particle_summary_all.csv"))
        
        fwrite(rbindlist(lapply(list.files(path = wd, pattern = "^particle_summary.*\\.csv$", full.names = T), fread), fill = T), 
               paste0(wd, "/particle_summary_all.csv"))
    }
    if("particle_details" %in% types & "all" %in% by){
        if(file.exists(paste0(wd, "/particle_details_all.csv"))) file.remove(paste0(wd, "/particle_details_all.csv"))
            
        fwrite(rbindlist(lapply(list.files(path = wd, pattern = "^particle_details.*\\.csv$", full.names = T), fread), fill = T), 
               paste0(wd, "/particle_details_all.csv"))
    }
    
    if("median_spec" %in% types & "all" %in% by){
        fwrite(rbindlist(lapply(list.files(path = wd, pattern = "^median_spec.*\\.csv$", full.names = T), fread), fill = T), 
               paste0(wd, "/median_spec_all.csv"))
    }
}

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

lib$metadata <- lib$metadata %>%
    mutate(material_class = ifelse(material_class %in% c("polyacrylamides", "polyamides (polylactams)"), "poly(acrylamide/amid)s", material_class)) %>%
    mutate(material_class = ifelse(material_class %in% c("polyterephthalates", "polyesters", "polyethers (polyglycols, oxyalkylenes, glycidyl ethers & cyclic ethers)", "polydiglycidyl ethers (polyepoxides, polyhydroxyethers, phenoxy)"), "poly(esters/ethers/diglycidylethers/terephthalates)s", material_class)) 

#Set this ----
#You'll set this wherever you put your zip folders. C:/Users/winco/OneDrive/Documents/Positive_Controls
wd <- "C:\\Users\\winco\\OneDrive\\Documents\\HPU\\Win" 

#Unzip files in batch ----
# zip_files <- list.files(path = wd,
#            pattern = "\\.zip$",
#            full.names = T)
# 
# for(zip_file in zip_files){
#     print(zip_file)
#     unzip(zipfile = zip_file,exdir=wd)  # unzip your file
# }


files <- list.files(path = wd, 
                    pattern = "(\\.dat$)|(\\.img$)",  
                    full.names = T)

files_hdr <- list.files(path = wd, 
                        pattern = ".hdr$",  
                        full.names = T)

img <- gsub(".dat", ".JPG", files)

gsub(".dat", ".hdr", files[!file.exists(gsub(".dat", ".hdr", files))])
gsub(".hdr", ".dat", files_hdr[!file.exists(gsub(".hdr", ".dat", files_hdr))])
img[!file.exists(img)]

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
    analyze_features(files = files, 
                     lib = lib,
                     img = img,#"C:\\Users\\winco\\OneDrive\\Documents\\Positive_Controls\\Mosaic Image files\\Recovery_red_beads_75-90um_5um-screen.JPG",
                     #bottom_left = list(c(171, 472)),
                     #top_right = list(c(632, 10)),
                     #origins = NULL, 
                     spectral_smooth = F, 
                     #sigma1 = c(0.001,2),
                     #sigma2 = c(0.001,2),
                     #spatial_smooth = T,
                     close = F,
                     adj_map_baseline = F, #Probably needs to be updated somehow because currently very small change and decreasing id accuracy. 
                     sn_threshold = 0.6, #0.04 the lowest signal considered "particles"
                     cor_threshold = 0.66, #the lowest level you would consider something known
                     area_threshold = 1, #to remove artefacts from background. 
                     label_unknown = F, #whether to scrub information when below correlation threshold.
                     remove_nonplastic = F, #whether to remove all nonplastic ided particles
                     remove_unknown = F,
                     pixel_length = 25, #25 conversion from pixels to length, set up for microns
                     metric = "sig_times_noise", #technique for the signal
                     abs = T, #F This determines whether to take the absolute value of the metric.
                     particle_id_strategy = "collapse", #This describes how the matching happens. 
                     vote_count = 10,
                     collapse_function = median,
                     k = 1,
                     k_weighting = "mean", #or multiple
                     wd = wd, #will put results here. 
                     types = c("particle_heatmap_thresholded",
                               "particle_heatmap", 
                               "particle_details",
                               "particle_summary",
                               "particle_image",
                               "spectra_processed"),
                     by = c("sample", 
                            "all"),
                     width = 1000, 
                     height = 1000, 
                     units = "px")

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




