# libraries
library(maps)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra) # raster plotting
library(tidyr)
library(scales)
library(dplyr)
library(fmesher)
library(INLA)
library(inlabru)
library(INLAutils)
# Note: the 'splancs' package also needs to be installed,
# but doesn't need to be loaded

# Bibliothèques
source(file.path("code", "function.R", fsep = "/"))
# set option
select <- dplyr::select
options(scipen = 99999)
options(max.print = 99999)
options(stringsAsFactors = FALSE)

#  Load the datas
rep <- "data"
load(file.path(rep, "df_FRANCE", fsep = "/"))

df_FRANCE_subset <- df_FRANCE %>% select(Cd_mousse, EMEP_air, lon, lat)

df_FRANCE_count <- df_FRANCE_subset %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    st_transform("EPSG:2154")
df_FRANCE_count$site_idx <- row.names(df_FRANCE_count)


# # Charger les packages
# library(spdep)
# library(sf)
# library(dplyr)

# # Charger les données spatiales
# df_FRANCE_count <- df_FRANCE_count %>%
#     st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
#     st_transform("EPSG:2154")

# # Calculer la matrice des voisins basée sur une distance
# coords <- st_coordinates(df_FRANCE_count)
# nb <- dnearneigh(coords, 0, 100000)  # Ajustez la distance selon vos besoins

# # Calculer les poids spatiaux
# listw <- nb2listw(nb, style = "W")

# # Calculer l'indice de Moran pour Cd_mousse
# moran_test <- moran.test(df_FRANCE_count$Cd_mousse, listw)

# # Afficher les résultats
# print(moran_test)



FRANCE <- terra::vect("data/regbiofr/region_biogeographique.shp")
France <- FRANCE[FRANCE$CODE == "ALP" | FRANCE$CODE == "CON" | FRANCE$CODE == "ATL" | FRANCE$CODE == "MED"]

load(file.path(rep, "df_france_cd_air_2016", fsep = "/"))
cd_air_2016 <- df_france_cd_air_2016 %>% select(-c("i", "j"))
cd_air_2016_rast <- terra::rast(cd_air_2016,
                                type = "xyz",
                                crs = "EPSG:4326",
                                digits = 6) %>% 
    terra::project("EPSG:2154")
cd_air_2016_rast <- terra::crop(cd_air_2016_rast, France, mask = TRUE)
r2 <- cd_air_2016_rast
terra::res(r2) <- 2000
cd_air_2016_rast <- terra::resample(cd_air_2016_rast, r2)
cd_air_2016_rast <- terra::crop(cd_air_2016_rast, France, mask = TRUE)

# France base with gographical areas
states <- st_as_sf(France) %>%
    filter(CODE %in% c(
        "ALP", "ATL", "CON", "MED"
    ))

# 445 concentrations measures on France map with geographical areas
ggplot() +
    geom_sf(
        data = df_FRANCE_count,
        aes(col = Cd_mousse)
    ) +
    geom_sf(data = states, fill = NA) +
    coord_sf(datum = NA) +
    scale_color_distiller(palette = "Spectral") +
    theme_bw()

# make a set of distinct study sites for mapping
site_map <- df_FRANCE_count %>%
    select(Cd_mousse, lon, lat) %>%
    distinct() %>%
    select(Cd_mousse, lon, lat)

# make a mesh restricted to France geography

st_coordinates(df_FRANCE_count)
coords <- st_coordinates(df_FRANCE_count)
bdry <- inla.sp2segment(st_as_sf(as.polygons(cd_air_2016_rast)))
bdry$loc <- inla.mesh.map(bdry$loc)
new_mesh <- inla.mesh.2d(
                         loc = coords,
                         boundary = bdry,
                         max.edge = c(100000, 170000), # km inside and outside
                         cutoff = 2000,
                         offset = c(100000, 170000),
                         crs = fm_crs(as.matrix(df_FRANCE_count)))

# Mesh plot with measure sites
ggplot() +
    gg(data = new_mesh) +
    geom_sf(data = site_map, col = "darkgreen", size = 1) +
    theme_bw() +
    labs(x = "", y = "")


# make spde
spde <- inla.spde2.pcmatern(
    mesh = new_mesh,
    prior.range = c(60000, 0.05), # P(range < 60000) = 0.05
    prior.sigma = c(5, 0.05) # P(sigma > 5) = 0.05,
)


################################# Models and methods#######################################

############# Null model

res_null <- bru(
    components = log(Cd_mousse) ~ -1 + intercept(1) + EMEP_air,
    family = "gaussian",
    data = df_FRANCE_count,
    options = list(
        control.compute = list(waic = TRUE, cpo = FALSE),
        control.inla = list(int.strategy = "eb"),
        verbose = FALSE
    )
)
summary(res_null)
############# Standard spatial model

res_spatial <- bru(
    components = log(Cd_mousse) ~ -1 + intercept(1) + EMEP_air + u(geometry, model = spde),
    family = "gaussian",
    data = df_FRANCE_count,
    options = list(
        control.compute = list(waic = TRUE, cpo = FALSE),
        control.inla = list(int.strategy = "eb"),
        verbose = FALSE
    )
)
summary(res_spatial)

############# Spatial+ v1 model with EMEP fitted with bru
spde_X <- spde <- inla.spde2.pcmatern(
    mesh = new_mesh,
    prior.range = c(130000, 0.05), # P(range < 130000) = 0.05
    prior.sigma = c(13, 0.05) # P(sigma > 13) = 0.05,
)

res_EMEP <- bru(
    components = ~ -1 + u(geometry, model = spde_X),
    like(
        formula = EMEP_air ~ .,
        family = "gaussian",
        data = df_FRANCE_count
    ),
    options = list(
        control.compute = list(waic = TRUE, cpo = FALSE),
        control.inla = list(int.strategy = "eb"),
        verbose = FALSE
    )
)

true_EMEP <- df_FRANCE_count$EMEP_air

df_FRANCE_count$fitted_EMEP_air_bru <- true_EMEP - 
    res_EMEP$summary.fitted.values[0:445, ]$mean

res_spatial_plus1 <- bru(
    components = ~ -1 + intercept(1) +
        fitted_EMEP_air_bru + u(geometry, model = spde),
    like(
        formula = log(Cd_mousse) ~ .,
        family = "gaussian",
        data = df_FRANCE_count
    ),
    options = list(
        control.compute = list(waic = TRUE, cpo = FALSE),
        control.inla = list(int.strategy = "eb"),
        verbose = FALSE
    )
)
summary(res_spatial_plus1)

############# Spatial+ v1 model with EMEP fitted with lm
mod_EMEP_lm <- lm(EMEP_air ~ -1 + lon + lat, data = df_FRANCE_count)
# summary(mod_EMEP_lm$fitted.values)

true_EMEP <- df_FRANCE_count$EMEP_air

df_FRANCE_count$fitted_EMEP_air_lm <- true_EMEP - mod_EMEP_lm$fitted.values


res_spatial_plus1_lm <- bru(
    components = ~ -1 + intercept(1) +
        fitted_EMEP_air_lm + u(geometry, model = spde),
    like(
        formula = log(Cd_mousse) ~ .,
        family = "gaussian",
        data = df_FRANCE_count
    ),
    options = list(
        control.compute = list(waic = TRUE, cpo = FALSE),
        control.inla = list(int.strategy = "eb"),
        verbose = FALSE
    )
)
summary(res_spatial_plus1_lm)

############# Spatial+ v1 model with EMEP fitted with gam

library(mgcv)


k_sp=60

mod_EMEP_gam <-gam(EMEP_air ~ -1 + s(lon,lat, k=60), data=df_FRANCE_count)
# b<-mod_EMEP_gam
# plot(b,pages=1,residuals=TRUE)
# gam.check(b)
# summary(mod_EMEP_gam)

df_FRANCE_count$fitted_EMEP_air_gam <- true_EMEP - mod_EMEP_gam$fitted.values

res_spatial_plus1_gam_bru <- bru(
    components = ~ -1 + intercept(1) +
        fitted_EMEP_air_gam + u(geometry, model = spde),
    like(
        formula = log(Cd_mousse) ~ .,
        family = "gaussian",
        data = df_FRANCE_count
    ),
    options = list(
        control.compute = list(waic = TRUE, cpo = FALSE),
        control.inla = list(int.strategy = "eb"),
        verbose = FALSE
    )
)
summary(res_spatial_plus1_gam_bru)

res_spatial_plus1_gam_gam <- gam(
                                 log(Cd_mousse) ~ fitted_EMEP_air_gam +
                                     s(lon, lat, k = 160),
                                 data = df_FRANCE_count)
summary(res_spatial_plus1_gam_gam)

###### Spatial+ 2.0 model
hyperparam <- res_spatial$summary.hyperpar
log_range_u <- hyperparam["Range for u", "mean"]

dist <- sp::spDists(coordinates(df_FRANCE_subset[, c("lon", "lat")]), longlat = TRUE)

Sigma <- inla.matern.cov(nu = 1, kappa = sqrt(8) / log_range_u, dist)
# Q <- inla.spde2.precision(spde = spde, theta = c(log_range_u, log_sigma_u))
Q <- solve(Sigma)
r <- eigen(Q) # spectral decomposition
eigen.vect <- r$vectors # eigen vectors

# find the decomposition in eigen vector
coef <- solve(eigen.vect, df_FRANCE_count$EMEP_air)

fit_with_eigen_vectors <- function(nbr_eigen_vectors = 445,
                                   data = df_FRANCE_count,
                                   decomposition=coef) {
    col_name <- paste("EMEP_air.eigen", nbr_eigen_vectors, sep = "")
    if(nbr_eigen_vectors == 1){
        data[, col_name] <-  as.vector(eigen.vect[, 0:nbr_eigen_vectors] * decomposition[0:nbr_eigen_vectors])
    }else{
        data[, col_name] <-  as.vector(eigen.vect[, 0:nbr_eigen_vectors]%*%decomposition[0:nbr_eigen_vectors])
    }
    res <- bru(
        components = as.formula(
          paste0("~ -1 + intercept(1) + ",
                 col_name, " + u(geometry, model = spde)", sep = "")
        ),
        like(
            formula = log(Cd_mousse) ~ .,
            family = "gaussian",
            data = data
        ),
        options = list(
            control.compute = list(waic = TRUE, cpo = FALSE),
            control.inla = list(int.strategy = "eb"),
            verbose = FALSE
        )
    )
    return(res)
}

res_spatial_plus2 <- fit_with_eigen_vectors(341)
summary(res_spatial_plus2)

#####RSR
spde_RSR <- inla.spde2.pcmatern(
    mesh = new_mesh,
    prior.range = c(18000000, 0.9999), # P(range < 18000000) = 0.9999
    prior.sigma = c(2.5, 0.0001) # P(sigma > 2.5) = 0.0001,
)

res_RSR_with_i <- bru( components = ~ -1 + intercept(1) + X_Effect(EMEP_air) + u(geometry, model = spde_RRS),
    like(
        formula = log(Cd_mousse) ~ 0 + intercept + X_Effect + u  - (sum(u) * sum(EMEP_air^2) + EMEP_air * (length(EMEP_air) * sum(EMEP_air * u) - sum(u) * sum(EMEP_air)) - sum(EMEP_air * u) * sum(EMEP_air)) / (length(EMEP_air) * sum(EMEP_air^2) - sum(EMEP_air)^2), 
        family = "gaussian",
        data = df_FRANCE_count
    ),
    options = list(
            control.compute = list(waic = TRUE, cpo = FALSE),
            control.inla = list(int.strategy = "eb"),
            verbose = FALSE
        )
)
summary(res_RSR_with_i)
