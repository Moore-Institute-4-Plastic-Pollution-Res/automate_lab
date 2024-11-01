validate_inputs <- function(params) {
  with(params, {
    if (!all(particle_id_strategy %in% c("collapse", "all cell id", "particle cell vote")))
      stop("Incorrect particle_id_strategy")
    if (!all(is.character(files)))
      stop("Incorrect files")
    if (!all(is.character(img) | is.null(img)))
      stop("Incorrect img")
    if (!all(is.list(bottom_left) | is.null(bottom_left)))
      stop("Incorrect bottom_left")
    if (!all(is.list(top_right) | is.null(top_right)))
      stop("Incorrect top_right")
    if (!all(is_OpenSpecy(lib) | is.list(lib)))
      stop("Incorrect lib")
    if (!all(is.logical(adj_map_baseline)))
      stop("Incorrect adj_map_baseline")
    if (!all(is.character(material_class)))
      stop("Incorrect material_class")
    if (!all(is.list(origins)))
      stop("Incorrect origins")
    if (!all(is.logical(spectral_smooth)))
      stop("Incorrect spectral_smooth")
    if (!all(is.numeric(sigma1)))
      stop("Incorrect sigma1")
    if (!all(is.logical(spatial_smooth)))
      stop("Incorrect spatial_smooth")
    if (!all(is.numeric(sigma2)))
      stop("Incorrect sigma2")
    if (!all(is.logical(close)))
      stop("Incorrect close")
    if (!all(is.numeric(close_kernel)))
      stop("Incorrect close_kernel")
    if (!all(is.numeric(sn_threshold)))
      stop("Incorrect sn_threshold")
    if (!all(is.numeric(cor_threshold)))
      stop("Incorrect cor_threshold")
    if (!all(is.numeric(area_threshold)))
      stop("Incorrect area_threshold")
    if (!all(is.logical(label_unknown)))
      stop("Incorrect label_unknown")
    if (!all(is.logical(remove_nonplastic)))
      stop("Incorrect remove_nonplastic")
    if (!all(is.logical(remove_unknown)))
      stop("Incorrect remove_unknown")
    if (!all(is.numeric(pixel_length)))
      stop("Incorrect pixel_length")
    if (!all(is.character(metric)))
      stop("Incorrect metric")
    if (!all(is.logical(abs)))
      stop("Incorrect abs")
    if (!all(is.numeric(k)))
      stop("Incorrect k")
    if (!all(is.character(k_weighting)))
      stop("Incorrect k_weighting")
    if (!all(is.character(wd)))
      stop("Incorrect wd")
    if (!all(is.character(types)))
      stop("Incorrect types")
    if (!all(is.character(by)))
      stop("Incorrect by")
    if (!all(is.numeric(width)))
      stop("Incorrect width")
    if (!all(is.numeric(height)))
      stop("Incorrect height")
    if (!all(is.character(units)))
      stop("Incorrect units")
  })
}
