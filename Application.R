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
library(mgcv)

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

FRANCE <- terra::vect("data/regbiofr/region_biogeographique.shp")
France <- FRANCE[FRANCE$CODE == "ALP" | FRANCE$CODE == "CON" | FRANCE$CODE == "ATL" | FRANCE$CODE == "MED"]

load(file.path(rep, "df_france_cd_air_2016", fsep = "/"))
cd_air_2016 <- df_france_cd_air_2016 %>% select(-c("i", "j"))
cd_air_2016_rast <- terra::rast(cd_air_2016,
    type = "xyz",
    crs = "EPSG:4326",
    digits = 6
) %>%
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
    crs = fm_crs(as.matrix(df_FRANCE_count))
)

# Mesh plot with measure sites
ggplot() +
    gg(data = new_mesh) +
    geom_sf(data = site_map, col = "darkgreen", size = 1) +
    theme_bw() +
    labs(x = "", y = "")

# make spde
spde <- inla.spde2.pcmatern(
    # alpha = 1.000001,
    mesh = new_mesh,
    prior.range = c(60000, 0.05), # P(range < 60000) = 0.05
    prior.sigma = c(5, 0.05) # P(sigma > 5) = 0.05,
)

################################# Models implementation #######################################

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

############# Spatial model
### With inlabru
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

### With mgcv
res_spatial_gam <- gam(log(Cd_mousse) ~ EMEP_air + s(lon, lat, k = 250, fx = FALSE), data = df_FRANCE_count, method = "GCV.Cp")
summary(res_spatial_gam)

############# RSR model
### With inlabru
res_RSR_formula <- bru(
    components = ~ -1 + intercept(1) + X_Effect(EMEP_air) + u(geometry, model = spde),
    like(
        formula = log(Cd_mousse) ~ 0 + intercept + X_Effect + u - (sum(u) * sum(EMEP_air^2) + EMEP_air * (length(EMEP_air) * sum(EMEP_air * u) - sum(u) * sum(EMEP_air)) - sum(EMEP_air * u) * sum(EMEP_air)) / (length(EMEP_air) * sum(EMEP_air^2) - sum(EMEP_air)^2),
        family = "gaussian",
        data = df_FRANCE_count
    ),
    options = list(
        control.compute = list(waic = TRUE, cpo = FALSE),
        control.inla = list(int.strategy = "eb"),
        verbose = FALSE
    )
)
summary(res_RSR_formula)

############# RSR model with extraconstr
A.rsr <- inla.spde.make.A(new_mesh, loc = df_FRANCE_count$geometry)
e.rsr <- rep(0, nrow(A.rsr))
res_RSR_extra <- bru(
    components = ~ -1 + intercept(1) + X_Effect(EMEP_air) + u(geometry, model = spde, extraconstr = list(A = A.rsr, e = e.rsr)),
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

summary(res_RSR_extra)

############# Spatial+ model
### With inlabru
spde_X <- inla.spde2.pcmatern(
    mesh = new_mesh,
    prior.range = c(130000, 0.05), # P(range < 130000) = 0.05
    prior.sigma = c(13, 0.05) # P(sigma > 13) = 0.05,
)

res_EMEP <- bru(
    components = ~ -1 + intercept(1) + u(geometry, model = spde_X),
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

### With mgcv
f_X_hat <- gam(EMEP_air ~ s(lon, lat, k = 250), data = df_FRANCE_count, method = "GCV.Cp")
df_FRANCE_count$fitted_EMEP_air_bru <- df_FRANCE_count$EMEP_air - f_X_hat$fitted.values
res_gam <- gam(log(Cd_mousse) ~ fitted_EMEP_air_bru + s(lon, lat, k = 250, fx = FALSE), data = df_FRANCE_count, method = "GCV.Cp")
summary(res_gam)

############# Spatial+ 2.0 model
hyperparam <- res_spatial$summary.hyperpar
log_range_u <- hyperparam["Range for u", "mean"]

dist <- sp::spDists(coordinates(df_FRANCE_subset[, c("lon", "lat")]), longlat = TRUE)

Sigma <- inla.matern.cov(nu = 1, kappa = sqrt(8) / log_range_u, dist)
Q <- solve(Sigma)
r <- eigen(Q) # spectral decomposition
eigen.vect <- r$vectors # eigen vectors

# find the decomposition in eigen vector
coef <- solve(eigen.vect, df_FRANCE_count$EMEP_air)

fit_with_eigen_vectors <- function(nbr_eigen_vectors = 445,
                                   data = df_FRANCE_count,
                                   decomposition = coef) {
    col_name <- paste("EMEP_air.eigen", nbr_eigen_vectors, sep = "")
    if (nbr_eigen_vectors == 1) {
        data[, col_name] <- as.vector(eigen.vect[, 0:nbr_eigen_vectors] * decomposition[0:nbr_eigen_vectors])
    } else {
        data[, col_name] <- as.vector(eigen.vect[, 0:nbr_eigen_vectors] %*% decomposition[0:nbr_eigen_vectors])
    }
    res <- bru(
        components = as.formula(
            paste0("~ -1 + intercept(1) + ",
                col_name, " + u(geometry, model = spde)",
                sep = ""
            )
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

res_spatial_plus2 <- fit_with_eigen_vectors(350)
summary(res_spatial_plus2)

############# gSEM model
### With inlabru
f_X_hat_res <- bru(EMEP_air ~ -1 + intercept(1) + u(geometry, model = spde), data = df_FRANCE_count, family = "gaussian")
df_FRANCE_count$r_X <- df_FRANCE_count$EMEP_air - f_X_hat_res$summary.fitted.values[1:445, ]$mean

f_Y_hat_res <- bru(log(Cd_mousse) ~ -1 + intercept(1) + u(geometry, model = spde), data = df_FRANCE_count, family = "gaussian")
df_FRANCE_count$r_Y <- log(df_FRANCE_count$Cd_mousse) - f_Y_hat_res$summary.fitted.values[1:445, ]$mean

res_gSEM <- bru(r_Y ~ -1 + intercept(1) + r_X, family = "gaussian", data = df_FRANCE_count)

### With mgcv
f_X_hat <- gam(EMEP_air ~ s(lon, lat, k = 250, fx = FALSE), data = df_FRANCE_count, method = "GCV.Cp")
df_FRANCE_count$r_X <- df_FRANCE_count$EMEP_air - f_X_hat$fitted.values
f_Y_hat <- gam(log(Cd_mousse) ~ s(lon, lat, k = 250, fx = FALSE), data = df_FRANCE_count, method = "GCV.Cp")
df_FRANCE_count$r_Y <- log(df_FRANCE_count$Cd_mousse) - f_Y_hat$fitted.values
mod <- lm(r_Y ~ r_X, data = df_FRANCE_count)
summary(mod)

# Boxplots as a function of the number
make_boxplot_sensiblity <- function(max = 445, min = 1, precision = 10) {
    N <- round(seq(min, max, length.out = precision))
    sensibility_df <- data.frame(null = matrix(NA, nrow = 3, ncol = 1))

    for (n in N) {
        print(n)
        model <- fit_with_eigen_vectors(n)
        col_name <- n # paste("eigen", n, sep = "")
        sensibility_df[, col_name] <- as.vector(
            unlist(
                model$summary.fixed[2, 3:5]
            )
        )
    }

    sensibility_df <- subset(sensibility_df, select = -null)
    return(sensibility_df)
}

sensibility_bx <- make_boxplot_sensiblity(max = 445, min = 1, precision = 445)

# png(file = "bxp_sensi_eigen_data.png", width = 1000, height = 500)
boxplot(sensibility_bx,
    names = names(sensibility_bx),
    xlab = "k",
    ylab = "Beta Estimates"
    #  = paste("Distribution of Beta Estimates Across Different Numbers of Eigenvectors in X Decomposition with Spatial+ v2")
)
# dev.off()
