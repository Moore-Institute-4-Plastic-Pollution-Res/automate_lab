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
