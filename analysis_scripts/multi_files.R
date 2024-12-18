library(OpenSpecy)
library(ggrepel)
library(data.table)
library(mmand)
library(magick)
library(cluster)
library(googledrive)


#Set a folder to save data locally ---
local_store <- file.path("data",project_name)
local_store_results <- file.path("data", project_name, "Results", fsep = "/")
# Set base directory of files downloaded and analyzed
local_store_raw <- file.path("data", project_name, "Raw", fsep = "/")


# Create folders ----
if(!dir.exists(local_store)){
  dir.create(local_store, recursive = TRUE)
}
if(!dir.exists(local_store_results)){
  dir.create(local_store_results, recursive = TRUE)
}

# Download data 
analyze_features <- function(project_name,
                             folder = "Export_Files",
                             bottom_left = NULL,
                             top_right = NULL,
                             lib,
                             adj_map_baseline = F,
                             material_class = "material_class",
                             spectral_smooth = F,
                             sigma1 = c(1, 1, 1),
                             spatial_smooth = F,
                             sigma2 = c(3, 3),
                             close = F,
                             close_kernel = c(4, 4),
                             sn_threshold = 0.01,
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
                             wd = local_store_results,
                             particle_id_strategy = "collapse",
                             vote_count = Inf,
                             collapse_function = median,
                             types = c(
                               "particle_image",
                               "particle_heatmap",
                               "particle_heatmap_thresholded",
                               "median_spec_plot",
                               "median_spec",
                               "particle_details",
                               "particle_summary",
                               "spectra_raw",
                               "spectra_processed",
                               "time"
                             ),
                             by = c("material", "sample", "all"),
                             width = 850,
                             height = 850,
                             units = "px"
) {
  
  # Find spectra files
  # data_search <- shared_drive_find("Project") |> 
  #   drive_ls("Customer Projects") |> 
  #   drive_ls() |> 
  #   filter(name %in% project_name) |> 
  #   drive_ls(folder) |> 
  #   drive_ls()
  
  data_search <- shared_drive_find("Project") |> 
    drive_ls("Customer Projects") |> 
    drive_ls() |> 
    filter(name %in% project_name) |> 
    drive_ls("JHENG_270") |> 
    drive_ls("JHENG270") |> 
    drive_ls(folder) |> 
    drive_ls()
  
  # Delimit file names to group
  data_temp <- data_search |> 
    mutate(temp_name = name) |> 
    separate(col = temp_name, into = c("name_new", "file"), sep = "[.]") |> 
    distinct(name, .keep_all = TRUE)
  
  missing_files <- character()

  # Identify missing files
  for (unique_name in unique(data_temp$name_new)){
    
    matching_files <- data_temp |> 
      filter(name_new %in% unique_name)
    # Here is where we check if all the necessary files are available 
    # without downloading based on matching files
    
    # Check if they are present
    hdr_present <- any(grepl(".hdr$", matching_files$name))
    dat_present <- any(grepl(".dat$", matching_files$name))
    img_present <- any(grepl("\\.jpg$|\\.JPG$", matching_files$name))
    
    # Skip and log if any required file is missing
    if (!hdr_present || !dat_present || !img_present) {
      missing_files <- c(missing_files, unique_name)
      cat("Warning:", unique_name,
          if(!hdr_present) "missing .hdr file.\n",
          if(!dat_present) "missing .dat file.\n",
          if(!img_present) "missing image file (.JPG).\n")
    }
  }
  
  # Stop if the missing files list is not empty
  if (length(missing_files) > 0){
    stop("Missing files")
  }

    # Group files by name, download, and do analysis
  for (unique_name in unique(data_temp$name_new)){
    
    matching_files <- data_temp |> 
      filter(name_new %in% unique_name)
  
    # Create temp file 
    if (dir.exists(local_store_raw)){
      unlink(local_store_raw, recursive = TRUE)
    }
    dir.create(local_store_raw, recursive = TRUE)
    
    # Download every file in the grouped matching file pairs ---
    for (file in 1:nrow(matching_files)) {
      tryCatch({
        drive_download(matching_files[file, ], path = file.path(local_store_raw, matching_files$name[file]), overwrite = TRUE)
      }, error = function(e) {
        message("Error downloading file: ", matching_files$name[file], "\n", e)
        stop()
      })
    }
    
    # List spectra files from Raw Data folder (local_store_raw)
    file <- list.files(path = local_store_raw,
                       pattern = "(\\.dat$)|(\\.img$)",
                       full.names = T)
    
    img <- list.files(path = local_store_raw,
                      pattern = "(\\.JPG$)|(\\.jpg$)",
                      full.names = T)
    
    
    #Extracts the file names
    file_names <- gsub("(.*/)|(\\..{1,3}$)", "", file)
    
    
    origins = list(x = rep(0, length(file)), y = rep(0, length(file)))
    
    #Begins timer
    time_start = Sys.time()
    
    
    #Variable input tests ----
    if(!all(particle_id_strategy %in% c("collapse", "all cell id", "particle cell vote"))) stop("Incorrect particle_id_strategy")
    # if(!all(is.character(data_all))) stop("Incorrect files")
    # if(!all(is.character(img)|is.null(img))) stop("Incorrect img")
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
    
    #Read in the file ----  
    if (grepl(".dat", file)) {
      map <- read_envi(file, spectral_smooth = spectral_smooth, sigma = sigma1)
    }
    else{
      map <- read_any(file)
    }

    #Correct baseline if specified
    if (adj_map_baseline) {
      baseline <- median(as.matrix(map$spectra), na.rm = T)
      map$spectra <- map$spectra - baseline
    }
    
    #extracts origins if available.
    if (!is.null(origins)) {
      originx <- origins[[1]]
      originy <- origins[[2]]
    } else if ("description" %in% names(map$metadata)) {
      #Auto attempts to extract origins if description available
      origin = unique(map$metadata$description)
      originx = gsub(",.*", "", gsub(".*X=", "", origin)) |> as.numeric()
      originy = gsub(".*Y=", "", origin) |> as.numeric()
    } else{
      originx <- rep(0, length(file))
      originy <- rep(0, length(file))
    }
    
    #Estimate the sig noise value on the raw data
    map$metadata$snr <- sig_noise(
      map |>
        restrict_range(
          min = c(750, 2420),
          max = c(2200, 4000),
          make_rel = F
        ),
      #Making comparable to lda and mediod
      metric = metric,
      spatial_smooth = spatial_smooth,
      sigma = sigma2,
      abs = abs
    )
    
    #Create a heatmap ----
    if ("particle_heatmap" %in% types) {
      plot <- ggplot() +
        geom_raster(
          data = map$metadata,
          aes(
            x = x * pixel_length + originx,
            y = y * pixel_length + originy,
            fill = snr
          )
        ) +
        coord_fixed() +
        theme_bw() +
        labs(x = "X (um)", y = "Y (um)") +
        scale_fill_viridis_c()
      
      # save heatmap
      ggsave(
        filename = paste0("particle_heatmap_", file_names, ".png"),
        path = wd,
        plot = plot,
        width = width,
        height = height,
        units = units
      )
      
      print(plot)
      
      if (all(types == "particle_heatmap"))
        next
    }
    
    
    #Create a thresholded heatmap ----
    if ("particle_heatmap_thresholded" %in% types) {
      plot <- ggplot() +
        geom_raster(
          data = map$metadata,
          aes(
            x = x * pixel_length + originx,
            y = y * pixel_length + originy,
            fill = snr > sn_threshold
          )
        ) +
        scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +  # Assign black to TRUE and white to FALSE
        coord_fixed() +
        theme_bw() +
        theme(legend.position = "none") +
        labs(x = "X (um)", y = "Y (um)")
      
      ggsave(
        filename = paste0("particle_heatmap_thresholded", file_names, ".jpg"),
        path = wd,
        plot = plot,
        width = width,
        height = height,
        units = units
      )
      print(plot)
      
      #Skip to the next map if only requesting thresholded outputs.
      if (all(types %in% c("particle_heatmap", "particle_heatmap_thresholded")))
        next
    }
    
    #Check if all or no features are thresholded and skip to the next if so.
    if (sum(map$metadata$snr > sn_threshold) == 0 |
        all(map$metadata$snr > sn_threshold))
      next
    
    #Auto detect the in the images using the red boxes.
    if (is.null(bottom_left[[file]]) & !is.null(img)) {
      
      # Catch images that cannot be read
      mosaic <- tryCatch(
        {
          jpeg::readJPEG(img)
        },
        error = function(e){
          message("Error": e$message)
          NULL
        },
        warning = function(w){
          message("Warning:", w$message)
          NULL
        }
      )
      
      # Identify red pixels based on the specified condition
      red_pixels <- mosaic[, , 1] * 255 > 50 &
        mosaic[, , 1] > 2 * mosaic[, , 2] & mosaic[, , 1] > 2 * mosaic[, , 3]
      
      # Get the x and y coordinates of the red pixels
      red_coords <- which(red_pixels, arr.ind = TRUE)
      x_coords <- red_coords[, 2]
      y_coords <- red_coords[, 1]
      
      # Calculate histograms of the x and y coordinates
      x_hist <- sort(table(x_coords), decreasing = T)
      y_hist <- sort(table(y_coords), decreasing = T)
      
      # Find the two x and y values with the highest counts
      x_boundaries <- as.integer(names(x_hist[c(TRUE, abs(diff(as.numeric(
        names(x_hist)
      ))) > 10)][1:2]))
      y_boundaries <- as.integer(names(y_hist[c(TRUE, abs(diff(as.numeric(
        names(y_hist)
      ))) > 10)][1:2]))
      
      # Sort the boundaries
      x_min <- min(x_boundaries)
      x_max <- max(x_boundaries)
      y_min <- max(y_boundaries)
      y_max <- min(y_boundaries)
      
      bl = c(x_min, y_min)
      tr = c(x_max, y_max)
      
    } else{
      bl = bottom_left[[file]]
      tr = top_right[[file]]
    }
    
    #Identify every pixel and use that to identify the features.
    if (particle_id_strategy == "all cell id") {
      proc_map  <- process_spec(
        filter_spec(map, map$metadata$snr > sn_threshold),
        conform_spec = T,
        conform_spec_args = list(range = if (!is_OpenSpecy(lib)) {
          lib$all_variables
        } else {
          lib$wavenumber
        }, res = NULL),
        restrict_range = T,
        restrict_range_args = list(min = c(750, 2420), max = c(2200, 4000)),#, #Making comparable to lda and mediod),
        #flatten_range = T #Improved ids by 7 %
        
        print(proc_map)
      )
      
    } else{
      #Just use the particle thresholds for identifying where the particles are.
      
      map$metadata$threshold <- map$metadata$snr > sn_threshold
      
      map <- def_features(
        map,
        features = map$metadata$threshold,
        shape_kernel = sigma2,
        shape_type = "box",
        close = close,
        close_kernel = close_kernel,
        close_type = "box",
        img,
        bl,
        tr
      )
    }
    
    # Check if 'feature_id' exists in 'map$metadata'
    if (!"feature_id" %in% names(map$metadata)) {
      stop("feature_id column is missing in map$metadata")
    }
    
    #Collapse the particles to their median. Remove -88 and areas with only 1 pixel. Process the spectra.
    if (particle_id_strategy == "collapse") {
      # need to def features and then collapse spec?
      # map <- def_features(map)
      
      if (sum(map$metadata$feature_id != "-88" &
              map$metadata$area > area_threshold) == 0)
        next
      proc_map  <- collapse_spec(map, fun = collapse_function) %>%
        filter_spec(.,
                    .$metadata$feature_id != "-88" &
                      .$metadata$area > area_threshold) %>%
        process_spec(
          .,
          conform_spec = T,
          conform_spec_args = list(range = if (!is_OpenSpecy(lib)) {
            lib$all_variables
          } else {
            lib$wavenumber
          }, res = NULL),
          restrict_range = T,
          restrict_range_args = list(min = c(800, 2420), max = c(2200, 3200))#)#, #Making comparable to lda and mediod),
          
          #flatten_range = T #Improved ids by 7 %
        )
    }
    
    
    #Process data for particle cell voting, don't compress.
    if (particle_id_strategy == "particle cell vote") {
      if (sum(map$metadata$feature_id != "-88" &
              map$metadata$area > area_threshold) == 0)
        next
      
      proc_map  <- map %>%
        filter_spec(.,
                    .$metadata$feature_id != "-88" &
                      .$metadata$area > area_threshold) %>%
        process_spec(
          .,
          conform_spec = T,
          conform_spec_args = list(range = if (!is_OpenSpecy(lib)) {
            lib$all_variables
          } else {
            lib$wavenumber
          }, res = NULL),
          restrict_range = T,
          restrict_range_args = list(min = c(800, 2420), max = c(2200, 3200))#, #Making comparable to lda and mediod),
          #flatten_range = T #Improved ids by 7 %
        )
    }
    
    #Run AI based classification
    if (!is_OpenSpecy(lib)) {
      cors <- match_spec(proc_map, lib)
      proc_map$metadata$max_cor_val <- cors$value
      proc_map$metadata$material_class <- cors$name
    } else if (k > 1) {
      #Run multi-k classification
      
      if (k_weighting == "raw") {
        cors <- match_spec(proc_map,
                           lib,
                           top_n = k,
                           add_library_metadata = "sample_name")  %>%
          group_by(object_id, material_class) %>%
          summarise(prob = n() / k) %>%
          group_by(object_id) %>%
          top_n(prob, n = 1) %>%
          ungroup() %>%
          rename(max_cor_val = prob)
      }
      
      if (k_weighting == "sum") {
        cors <- match_spec(proc_map,
                           lib,
                           top_n = k,
                           add_library_metadata = "sample_name")  %>%
          group_by(object_id, material_class) %>%
          summarise(prob = (n() / k + mean(match_val)) / 2) %>%
          group_by(object_id) %>%
          top_n(prob, n = 1) %>%
          ungroup() %>%
          rename(max_cor_val = prob)
      }
      if (k_weighting == "multiple") {
        cors <- match_spec(proc_map,
                           lib,
                           top_n = k,
                           add_library_metadata = "sample_name")  %>%
          group_by(object_id, material_class) %>%
          summarise(prob = n() / k * mean(match_val)) %>%
          group_by(object_id) %>%
          top_n(prob, n = 1) %>%
          ungroup() %>%
          rename(max_cor_val = prob)
      }
      
      proc_map$metadata <- left_join(proc_map$metadata, cors, by = c("feature_id" = "object_id"))
      
    } else{
      #Run standard top 1 k classification
      cors <- cor_spec(proc_map, lib)
      
      max_cors <- max_cor_named(cors)
      
      proc_map$metadata$max_cor_val <- max_cors
      
      proc_map$metadata$max_cor_name <- names(max_cors)
      
      proc_map$metadata <- left_join(
        proc_map$metadata,
        lib$metadata %>%
          dplyr::select(sample_name, material_class),
        by = c("max_cor_name" = "sample_name")
      )
      
    }
    #If using all cell id set the thresholds and identify the features.
    if (particle_id_strategy == "all cell id") {
      map$metadata <- left_join(
        map$metadata,
        dplyr::select(proc_map$metadata, x, y, max_cor_val, material_class),
        by = c("x", "y")
      ) %>%
        mutate(threshold = snr > sn_threshold &
                 max_cor_val > cor_threshold) %>%
        mutate(material_class = ifelse(threshold, material_class, "background"))
      
      map <- def_features(
        x = map,
        features = map$metadata$material_class,
        close = T,
        close_kernel = c(4, 4),
        close_shape = "box"
      )
      
      map$metadata <- map$metadata %>%
        mutate(material_class = ifelse(
          !grepl("background", feature_id) &
            area > area_threshold,
          gsub("_.*", "", feature_id),
          NA
        ))
      
      proc_map <- filter_spec(map, !is.na(map$metadata$feature_id)) %>%
        collapse_spec(fun = collapse_function)
      
    }
    
    #If using particle cell voting conduct the voting process
    if (particle_id_strategy == "particle cell vote") {
      ts <- transpose(proc_map$spectra)
      ts$id <- proc_map$metadata$feature_id
      proc_map$spectra <- ts[, lapply(.SD, median, na.rm = T), by = "id"] |>
        transpose(make.names = "id")
      
      # Step 1: Subset proc_map$metadata and then select unique rows based on certain columns
      subset_metadata <- proc_map$metadata[, .(
        feature_id,
        area,
        perimeter,
        feret_max,
        feret_min,
        centroid_y,
        centroid_x,
        first_x,
        first_y,
        rand_x,
        rand_y
      )]
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
    if (label_unknown) {
      proc_map$metadata[[material_class]] <- ifelse(proc_map$metadata$max_cor_val < cor_threshold,
                                                    "unknown",
                                                    proc_map$metadata[[material_class]])
    }
    
    #Remove any nonplastic particles from further analysis
    if (remove_nonplastic) {
      if (sum(!proc_map$metadata[[material_class]] %in% c("mineral", "organic matter", "other material")) == 0)
        next
      proc_map <- filter_spec(
        proc_map,
        !proc_map$metadata[[material_class]] %in% c("mineral", "organic matter", "other material")
      )
    }
    
    #Remove any nonplastic particles from further analysis
    if (remove_unknown) {
      if (sum(!proc_map$metadata[[material_class]] %in% c("unknown")) == 0)
        next
      proc_map <- filter_spec(proc_map, !proc_map$metadata[[material_class]] %in% c("unknown"))
    }
    
    if (particle_id_strategy %in% c("particle cell vote", "collapse")) {
      map$metadata <- left_join(
        map$metadata,
        dplyr::select(proc_map$metadata, feature_id, material_class),
        by = "feature_id"
      )
    }
    
    #Plot overlay on mosaic
    if (!is.null(img)) {
      map_dim <- c(length(unique(map$metadata$x)), length(unique(map$metadata$y)))
      xscale = (tr[1] - bl[1]) / map_dim[1]
      yscale = (bl[2] - tr[2]) / map_dim[2]
      x_vals = as.integer(map$metadata$x * xscale + bl[1])
      y_vals = as.integer(bl[2] - map$metadata$y * yscale)
      particles  = map$metadata$threshold
      
      # Convert the image to a raster
      mosaic <- jpeg::readJPEG(img)
      image_raster <- as.raster(mosaic)
      
      # Create a matrix of coordinates for indexing
      coords <- cbind(y_vals[particles], x_vals[particles])
      
      image_raster[coords] <- "#00FF00"
      
      png(
        paste0(wd, "/overlay_", file_names, ".png"),
        height = nrow(image_raster),
        width = ncol(image_raster)
      )
      
      plot(image_raster)
      dev.off()
    }
    
    #Plot sample level figures for particle types
    if ("particle_image" %in% types & "sample" %in% by) {
      plot <- particle_image(
        proc_map$metadata,
        map$metadata,
        material_class,
        pixel_length,
        origin = c(originx, originy),
        title = "All Particles"
      )
      ggsave(
        filename = paste0("particle_image_", file_names, ".png"),
        path = wd,
        plot = plot,
        width = width,
        height = height,
        units = units
      )
    }
    
    #Plot material-sample level figures for particle types
    if ("particle_image" %in% types & "material" %in% by) {
      for (material in unique(proc_map$metadata[[material_class]])) {
        map_subset <- map$metadata %>%
          mutate(!!material_class := ifelse(.data[[material_class]] == material, .data[[material_class]], NA))
        proc_map_subset <- proc_map$metadata %>%
          dplyr::filter(.data[[material_class]] == material)
        plot <- particle_image(
          proc_map_subset,
          map_subset,
          material_class,
          pixel_length,
          origin = c(originx, originy),
          title = paste0(material, " Particles")
        )
        ggsave(
          filename = paste0(material, "_particle_image_", file_names, ".png"),
          plot = plot,
          path = wd,
          width = width,
          height = height,
          units = units
        )
      }
    }
    
    #Plot some example spectra for each material type in each sample
    if ("median_spec_plot" %in% types & "sample" %in% by) {
      metadata_to_use <- proc_map$metadata %>%
        group_by(material_class) %>%
        sample_n(10, replace = T) %>%
        distinct() %>%
        ungroup()
      
      spectra_sample <- proc_map$spectra %>%
        mutate(wavenumber = proc_map$wavenumber) %>%
        pivot_longer(
          cols = -wavenumber,
          names_to = "feature_id",
          values_to = "intensity"
        ) %>%
        inner_join(metadata_to_use %>% dplyr::select(feature_id, material_class),
                   by = "feature_id")
      
      spectra_plot <- ggplot(spectra_sample) +
        geom_line(aes(x = wavenumber, y = intensity, group = feature_id)) +
        facet_grid(rows = vars(.data[[material_class]])) +
        theme_bw(base_size = 5) +
        labs(
          title = paste0("Example spectra in ", file_names[file]),
          x = "Wavenumbers (1/cm)",
          y = "Absorbance (Min-Max Norm)"
        ) +
        scale_x_reverse() +
        theme(plot.title = element_text(size = 5))
      
      ggsave(
        filename = paste0("median_spec_plot_", file_names, ".png"),
        plot = spectra_plot,
        path = wd,
        width = width,
        height = height,
        units = units
      )
    }
    
    #Create a unique figure for each material-sample spectra examples.
    if ("median_spec_plot" %in% types & "material" %in% by) {
      for (material in unique(proc_map$metadata[[material_class]])) {
        metadata_to_use <- proc_map$metadata %>%
          filter(.data[[material_class]] == material) %>%
          sample_n(10, replace = T) %>%
          distinct() %>%
          ungroup()
        
        spectra_sample <- proc_map$spectra %>%
          mutate(wavenumber = proc_map$wavenumber) %>%
          pivot_longer(
            cols = -wavenumber,
            names_to = "feature_id",
            values_to = "intensity"
          ) %>%
          inner_join(metadata_to_use %>% dplyr::select(feature_id, material_class),
                     by = "feature_id")
        
        spectra_plot <- ggplot(spectra_sample) +
          geom_line(aes(
            x = wavenumber,
            y = intensity,
            group = feature_id
          )) +
          facet_grid(rows = vars(feature_id)) +
          theme_bw(base_size = 5) +
          labs(
            title = paste0("Example Spectra ", material, " in ", file_names),
            x = "Wavenumbers (1/cm)",
            y = "Absorbance (Min-Max Norm)"
          ) +
          scale_x_reverse() +
          theme(plot.title = element_text(size = 5))
        
        ggsave(
          filename = paste0(material, "_median_spec_plot_", file_names, ".png"),
          plot = spectra_plot,
          path = wd,
          width = width,
          height = height,
          units = units
        )
      }
    }
    
    #Create a dataset for the output.
    if ("particle_details" %in% types) {
      particle_info <- proc_map$metadata %>%
        rename(particle_id = feature_id) %>%
        mutate(
          area = area * pixel_length ^ 2,
          perimeter = perimeter * pixel_length,
          feret_min = feret_min * pixel_length,
          feret_max = feret_max * pixel_length,
          centroid_x = centroid_x * pixel_length + originx,
          centroid_y = centroid_y * pixel_length + originy,
          first_x = first_x * pixel_length + originx,
          #dont work with vote option yet.
          first_y = first_y * pixel_length + originy,
          #rand_x = rand_x * pixel_length + originx,
          #rand_y = rand_y * pixel_length + originy,
          bad_spectra = max_cor_val < cor_threshold,
          acc_analy_conf = ifelse(
            max_cor_val > 0.6,
            "confident",
            ifelse(max_cor_val < 0.3, "undetermined", "possible")
          ),
          sample_id = file_names
        ) %>%
        mutate(
          aspect_ratio = feret_max / feret_min,
          circularity = (perimeter ^ 2) / (4 * pi * area)
        )
      
      selected_columns <- c(
        "particle_id",
        "sample_id",
        "max_cor_val",
        "bad_spectra",
        "material_class",
        "area",
        "perimeter",
        "feret_max",
        "feret_min",
        "aspect_ratio",
        "circularity",
        "centroid_x",
        "centroid_y",
        "first_x",
        "first_y",
        "acc_analy_conf"#,
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
        rename(
          max_length_um = feret_max,
          min_length_um = feret_min,
          perimeter_um = perimeter,
          area_um2 = area,
          centroid_y = centroid_y,
          centroid_x = centroid_x
        )
      fwrite(particle_info,
             paste0(wd, "/particle_details_", file_names, ".csv"))
      
    }
    
    #Create summary tables for each sample
    if ("particle_summary" %in% types) {
      particle_summary <- proc_map$metadata %>%
        rename(particle_id = feature_id) %>%
        group_by(material_class) %>%
        summarise(count = n()) %>%
        mutate(sample_id = file_names)
      
      fwrite(particle_summary,
             paste0(wd, "/particle_summary_", file_names, ".csv"))
      
    }
    
    #Save the raw spectra
    if ("spectra_raw" %in% types) {
      saveRDS(map, file = paste0(wd, "/particles_raw_", file_names, ".rds"))
    }
    
    #Save the processed spectra
    if ("spectra_processed" %in% types) {
      saveRDS(proc_map, file = paste0(wd, "/particles_", file_names, ".rds"))
    }
    
    #Save the median spectra as a csv
    if ("median_spec" %in% types) {
      spectra_sample <- proc_map$spectra %>%
        mutate(wavenumber = proc_map$wavenumber) %>%
        pivot_longer(
          cols = -wavenumber,
          names_to = "feature_id",
          values_to = "intensity"
        ) %>%
        inner_join(metadata_to_use %>% dplyr::select(feature_id, material_class),
                   by = "feature_id") %>%
        mutate(sample_id = file_names)
      
      fwrite(spectra_sample,
             paste0(wd, "/median_spec_", file_names, ".csv"))
    }
    #Print the final run time.
    time_diff = Sys.time() - time_start
    if ("time" %in% types) {
      saveRDS(time_diff, file = paste0(wd, "/time_", file_names, ".rds"))
    }
    
    #After all samples are finished, concatenate the summary, details, ann spectra csv tables.
    if ("particle_summary" %in% types & "all" %in% by) {
      if (file.exists(paste0(wd, "/particle_summary_all.csv")))
        file.remove(paste0(wd, "/particle_summary_all.csv"))
      
      fwrite(rbindlist(lapply(
        list.files(
          path = wd,
          pattern = "^particle_summary.*\\.csv$",
          full.names = T
        ),
        fread
      ), fill = T),
      paste0(wd, "/particle_summary_all.csv"))
    }
    if ("particle_details" %in% types & "all" %in% by) {
      if (file.exists(paste0(wd, "/particle_details_all.csv")))
        file.remove(paste0(wd, "/particle_details_all.csv"))
      
      fwrite(rbindlist(lapply(
        list.files(
          path = wd,
          pattern = "^particle_details.*\\.csv$",
          full.names = T
        ),
        fread
      ), fill = T),
      paste0(wd, "/particle_details_all.csv"))
    }
    
    if ("median_spec" %in% types & "all" %in% by) {
      fwrite(rbindlist(lapply(
        list.files(
          path = wd,
          pattern = "^median_spec.*\\.csv$",
          full.names = T
        ),
        fread
      ), fill = T),
      paste0(wd, "/median_spec_all.csv"))
    }
    # Export files ----
    # Find id of Project folder
    # folder_id <- shared_drive_find("Project") |>
    #   drive_ls("Customer Projects") |>
    #   drive_ls() |>
    #   filter(name == project_name) |>
    #   drive_ls() |>
    #   filter(name == "Export_Files") |>
    #   pull(id)
    # 
    # Create new folder
    
    folder_id <- "1M8c8RXsk5ayeINRsP36FAWHB-7I2Zej7"
    
    new_folder <- drive_mkdir("Spectral_Results",
                              path = as_id(folder_id),
                              overwrite = TRUE
                              )
    # Upload files
    # Send data up to drive
    wd <-file.path("data/JHENG/Results")
    data_upload <- list.files(wd)

    for (file in data_upload){

      drive_upload(media = file.path(wd, file),
                   path = as_id(as.character(new_folder$id))
      )

    }
    
  }
}

# Run analysis ----
analyze_features(project_name = project_name,
                 lib = lib,
                 #img = img,
                 #bottom_left = list(c(171, 472)),
                 #top_right = list(c(632, 10)),
                 #origins = NULL, 
                 spectral_smooth = T, 
                 #sigma1 = c(0.001,2),
                 #sigma2 = c(0.001,2),
                 #spatial_smooth = ,
                 close = F,
                 adj_map_baseline = F, 
                 sn_threshold = 0.01, #0.01 the lowest signal considered "particles"
                 cor_threshold = 0.66, #the lowest level you would consider something known
                 area_threshold = 1, #to remove artifacts from background. 
                 label_unknown = F, #whether to scrub information when below correlation threshold.
                 remove_nonplastic = F, #whether to remove all non plastic ided particles
                 remove_unknown = F,
                 pixel_length = 25, #25 conversion from pixels to length, set up for microns
                 metric = "sig_times_noise", #technique for the signal
                 abs = F, #F This determines whether to take the absolute value of the metric.
                 particle_id_strategy = "collapse", #This describes how the matching happens. 
                 vote_count = 10,
                 collapse_function = median,
                 k = 1,
                 k_weighting = "mean", #or multiple
                 wd = local_store_results, #will put results here. 
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
