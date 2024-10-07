source("Application.R")
library(ggspatial)

################ France data
France <- terra::vect('data/regbiofr/region_biogeographique.shp')

France <- France[France$CODE == "ALP" | France$CODE == "CON" | France$CODE == "ATL" | France$CODE == "MED" ]

######### Grid

gridPCM <-  SpatialPosition::CreateGrid(w = France, resolution = 2000, returnclass = "sf")#2000 = 2km x 2km

############## Prepare new data ##############
new <- c()

new$lon <- gridPCM$COORDX
new$lat <- gridPCM$COORDY
new <- as.data.frame(new)
new_rast <- terra::rast(new, type ="xyz", crs = "EPSG:2154", digits = 6)

new_rast$EMEP_air <- add_rast(df_rast = new_rast, 
                              rast_to_add = cd_air_2016_rast)

new_data <- as.data.frame(new_rast$EMEP_air, xy=TRUE) %>%
    st_as_sf(coords = c("x", "y"), crs = 2154, remove = FALSE) %>%
    st_transform("EPSG:2154")

#Spatial+1
true_EMEP_new <- new_data$EMEP_air

res_EMEP <- bru(
    components = ~ -1 + u(geometry, model = spde),
    like(
        formula = EMEP_air ~ .,
        family = "gaussian",
        data = new_data
    ),
    options = list(
        control.compute = list(waic = TRUE, cpo = FALSE),
        control.inla = list(int.strategy = "eb"),
        verbose = FALSE
    )
)

new_data$fitted_EMEP_air_bru <- true_EMEP_new - res_EMEP$summary.fitted.values[0:139199, ]$mean 

mod_EMEP_lm <- lm(EMEP_air ~ -1 + x + y, data = new_data)
new_data$fitted_EMEP_air_lm <- true_EMEP_new - mod_EMEP_lm$fitted.values

mod_EMEP_gam <-mgcv::gam(EMEP_air ~ s(x,y,k=60), data=new_data)
new_data$fitted_EMEP_air_gam <- true_EMEP_new - mod_EMEP_gam$fitted.values


# Spatial+2
# dist <- sp::spDists(st_coordinates(new_data[, c("x","y")]), longlat=TRUE)

# Sigma <- inla.matern.cov(nu = 1, kappa = sqrt(8) / log_range_u, dist)
# Q <- solve(Sigma)
# r <- eigen(Q) # spectral decomposition
# eigen.vect <- r$vectors # eigen vectors


# # find the decomposition in eigen vector
# coef <- solve(eigen.vect, new_data$EMEP_air)

# new_data$EMEP_air.eigen341 <- as.vector(eigen.vect[, 0:310]%*%coef[0:310])


############## Plot maps of predictions ##############

make_prediction_map <- function(pred_data, title, min = 0, max = 0.55, scale_title="prediction Cd (µg/g)\n mean") {
    if(class(pred_data)[1] == "bru_prediction") {
        pred_data$var <- pred_data$sd^2
        data_plot <- gg(pred_data, aes(fill = mean), geom = "tile")
    } else {
        data_plot <- tidyterra::geom_spatraster(data = pred_data, maxcell = 5e+07)
    }
    ggplot() +
        data_plot +
        scale_fill_distiller(scale_title,
            palette = "Spectral",
            na.value = "Transparent",
            limits = c(min, max)
        ) +
        annotation_scale(
            location = "bl", line_width = .5, height = unit(0.1, "cm"),
            tick_height = unit(0.1, "cm")
        ) +
        annotation_north_arrow(location = "tl", style = north_arrow_nautical) +
        coord_sf(crs = 2154, datum = sf::st_crs(2154)) +
        theme_classic() + # themes
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            plot.title = element_text(size = 16),
            legend.title = element_text(
                size = 14,
                color = "#4e4d47"
            ),
            legend.key.size = unit(1, "cm"),
            legend.text = element_text(
                size = 12,
                margin = margin(t = 3)
            ),
            plot.subtitle = element_text(
                size = 15,
                family = "Avantgarde",
                color = "#4e4d47"
            )
        )+
        ggtitle(title)
}


pred_null <- predict(res_null, newdata = new_data, formula = ~ exp(intercept + EMEP_air))
map_null <- make_prediction_map(pred_null, "Null", max=0.19)
# ggsave(map_null,
#        filename = "plot/map_sd_null.png", 
#        device='png', 
#        height=7,
#        width=7,
#        units="in",
#        dpi=300,
#        bg = "white"
# )

pred_spatial <- predict(res_spatial, newdata = new_data, formula = ~ exp(intercept + EMEP_air + u))
map_spatial <- make_prediction_map(pred_spatial, "Spatial" ,max=0.01)
# ggsave(map_spatial,
#        filename = "plot/map_sd_spatial_001.png", 
#        device='png', 
#        height=7,
#        width=7,
#        units="in",
#        dpi=300,
#        bg = "white"
# )

pred_RSR <- predict(res_RSR_with_i, newdata = new_data, formula = ~ exp(intercept + X_Effect + u  - (sum(u) * sum(EMEP_air^2) + EMEP_air * (length(EMEP_air) * sum(EMEP_air * u) - sum(u) * sum(EMEP_air)) - sum(EMEP_air * u) * sum(EMEP_air)) / (length(as.vector(EMEP_air)) * sum(EMEP_air^2) - sum(EMEP_air)^2)))
map_RSR <- make_prediction_map(pred_RSR, "RSR with intercept", max=0.19)
# ggsave(map_RSR,
#        filename = "plot/map_sd_RSR.png", 
#        device='png', 
#        height=7,
#        width=7,
#        units="in",
#        dpi=300,
#        bg = "white"
# )


pred_spatial_plus1_bru <- predict(res_spatial_plus1, newdata = new_data, formula = ~ exp(intercept + fitted_EMEP_air_bru + u))
map_spatial_plus1_bru_max45 <- make_prediction_map(pred_spatial_plus1_bru, "Spatial+ 1 GF", max=0.45)
map_spatial_plus1_bru <- make_prediction_map(pred_spatial_plus1_bru, "Spatial+ 1 GF", max=0.19)
# ggsave(map_spatial_plus1_bru,
#        filename = "plot/map_sd_spp1_GF.png", 
#        device='png', 
#        height=7,
#        width=7,
#        units="in",
#        dpi=300,
#        bg = "white"
# )

pred_spatial_plus1_lm <- predict(res_spatial_plus1_lm, newdata = new_data, formula = ~ exp(intercept + fitted_EMEP_air_lm + u))
map_spatial_plus1_lm_max45 <- make_prediction_map(pred_spatial_plus1_lm, "Spatial+ 1 lm", max= 0.45)
map_spatial_plus1_lm <- make_prediction_map(pred_spatial_plus1_lm, "Spatial+ 1 lm", max=0.01)
# ggsave(map_spatial_plus1_lm,
#        filename = "plot/map_sd_spp1_lm_001.png", 
#        device='png', 
#        height=7,
#        width=7,
#        units="in",
#        dpi=300,
#        bg = "white"
# )

pred_spatial_plus1_gam <- predict(res_spatial_plus1_gam_bru, newdata = new_data, formula = ~ exp(intercept + fitted_EMEP_air_gam + u))
map_spatial_plus1_gam_max45 <- make_prediction_map(pred_spatial_plus1_gam, "Spatial+ 1 gam", max=0.45)
map_spatial_plus1_gam <- make_prediction_map(pred_spatial_plus1_gam, "Spatial+ 1 gam", max= 0.01)
# ggsave(map_spatial_plus1_gam,
#        filename = "plot/map_sd_spp1_gam_001.png", 
#        device='png', 
#        height=7,
#        width=7,
#        units="in",
#        dpi=300,
#        bg = "white"
# )

###### all models
multiplot(map_null, map_RSR, map_spatial, map_spatial_plus1_bru, cols=3)
dev.print(png,
    file = "plot/prediction_maps_sd.png",
    width = 4000,
    height = 3000,
    units = "px"
)