get_lib <- function(type = c("derivative", 
                             "nobaseline", 
                             "raw", 
                             "medoid",
                             "model"),
                    path = "system",
                    ...) {
  
  lp <- ifelse(path == "system",
               system.file("extdata", package = "OpenSpecy"),
               path)
  
  # List all .rds files in the directory
  rds_files <- list.files(lp, pattern = "\\.rds$", full.names = TRUE)
  
  # Check if any .rds files exist
  if (length(rds_files) == 0) {
    
    message("Fetching Open Specy reference libraries from OSF ...")
    
    if("derivative" %in% type){
      message("Fetching derivative library...")
      download.file("https://osf.io/download/2qbkt/", destfile = file.path(lp, "derivative.rds"), mode = "wb")
    }
    
    if("nobaseline" %in% type){
      message("Fetching nobaseline library...")
      download.file("https://osf.io/download/jy7zk/", destfile = file.path(lp, "nobaseline.rds"), mode = "wb")
    }
    
    if("medoid" %in% type){
      message("Fetching mediod library...")
      download.file("https://osf.io/download/yzscg/", destfile = file.path(lp, "mediod.rds"), mode = "wb")
    }
    
    if("model" %in% type){
      message("Fetching model library...")
      download.file("https://osf.io/download/v2yr3/", destfile = file.path(lp, "model.rds"), mode = "wb")
    }
    
    if("raw" %in% type){
      message("Fetching raw library...")
      download.file("https://osf.io/download/kzv3n/", destfile = file.path(lp, "raw.rds"), mode = "wb")
    }
    
    message("Use 'load_lib()' to load the library")
    
    
  } else {
    print("Already Installed")
  }
  
}
