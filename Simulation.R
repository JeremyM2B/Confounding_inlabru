## Libraries
library(INLA)
library(inlabru)
library(dplyr)
library(ggplot2)
library(fields)
library(furrr)
library(knitr)
library(mgcv)

## Fixed parameters
n <- 500
n_sim <- 50
beta <- 3

## Create mesh
fake.locations <- matrix(c(0, 0, 10, 10, 0, 10, 10, 0), nrow = 4, byrow = T)
crs1 <- inla.CRS("EPSG:4326")
mesh.sim <- inla.mesh.2d(loc = fake.locations, max.edge = c(0.61, 1), crs = crs1)

ggplot() +
  gg(mesh.sim)

## SPDE
# Simulate spatial field
seed <- 7

spde_xs <- inla.spde2.pcmatern(mesh.sim, alpha = 2, prior.range = c(.5, .5), prior.sigma = c(.5, .5))

Qu_xs <- inla.spde.precision(spde_xs, theta = c(2, .4))

u_xs <- inla.qsample(n = 1, Q = Qu_xs, seed = seed)

spde_s <- inla.spde2.pcmatern(mesh.sim, alpha = 2, prior.range = c(.5, .5), prior.sigma = c(.5, .5))

Qu_s <- inla.spde.precision(spde_s, theta = c(0, 1))

u_s <- inla.qsample(n = 1, Q = Qu_s, seed = seed)

## Plot fields
local.plot.field(u_xs[, 1] - mean(u_xs[, 1]), mesh.sim)
local.plot.field(u_s[, 1] - mean(u_s[, 1]), mesh.sim)

## Simulate at measurement locations
set.seed(123)

u_s_xs <- matrix(data = NA, nrow = n, ncol = n_sim)
loc.data <- c()
for (i in 1:n_sim) {
  set.seed(i)
  loc.data[[i]] <- matrix(runif(2 * n), n) * 10
  A <- inla.spde.make.A(mesh = mesh.sim, loc = loc.data[[i]])
  u_s_xs[, i] <- drop(A %*% u_xs[, 1])
}

u_s_s <- matrix(data = NA, nrow = n, ncol = n_sim)
for (i in 1:n_sim) {
  set.seed(i + n_sim)
  B <- inla.spde.make.A(mesh = mesh.sim, loc = loc.data[[i]])
  u_s_s[, i] <- drop(B %*% u_s[, 1])
}

quilt.plot(
  x = loc.data[[1]][, 1], y = loc.data[[1]][, 2], z = u_s_xs[, 1], nx = 80, ny = 80,
  col = plasma(101), main = "Field projected to data locations",
  zlim = range(u_s_xs)
)

quilt.plot(
  x = loc.data[[1]][, 1], y = loc.data[[1]][, 2], z = u_s_s[, 1], nx = 80, ny = 80,
  col = plasma(101), main = "Field projected to data locations",
  zlim = range(u_s_s)
)

test <- c()
test$lon <- loc.data[[1]][, 1]
test$lat <- loc.data[[1]][, 2]
testdf <- as.data.frame(test)
testdf$u_s <- u_s_s[, 1]
testdf$u_xs <- u_s_xs[, 1]
testdf$dif <- -u_s_xs[, 1] - u_s_s[, 1]

testdfsf <- testdf %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  sf::st_transform("EPSG:4326")

plot(gstat::variogram(0.35 * (u_xs - mean(u_xs)) ~ 1, testdfsf))
plot(gstat::variogram(0.275 * (u_s - mean(u_s)) ~ 1, testdfsf))

## Generate random effects
generate_epsilon <- function() {
  rnorm(n, mean = 0, sd = 1)
}

epsilon_list <- lapply(1:n_sim, function(x) generate_epsilon())

epsilon <- as(do.call(cbind, epsilon_list), "dgeMatrix")
plot(density(epsilon_list[[1]]))

## Covariate x
generate_epsilon_x <- function() {
  rnorm(n, mean = 0, sd = 0.1)
}
epsilon_list_x <- lapply(1:n_sim, function(x) generate_epsilon_x())

epsilon_x <- as(do.call(cbind, epsilon_list_x), "dgeMatrix")

x <- 0.35 * (u_s_xs - mean(u_s_xs)) + epsilon_x

## Variable y
ux <- (u_s_xs - mean(u_s_xs))
us <- 0.275 * (u_s_s - mean(u_s_s))
Z <- us - ux
y <- beta * x + Z + epsilon

## Simulation dataframes
df_list <- lapply(1:n_sim, function(i) {
  data.frame(
    y = as.matrix(y)[, i],
    locx = loc.data[[i]][, 1],
    locy = loc.data[[i]][, 2],
    x = as.matrix(x)[, i],
    Z = as.matrix(Z)[, i]
  )
})

df_sf_list <- lapply(df_list, function(df) {
  df %>%
    sf::st_as_sf(coords = c("locx", "locy"), crs = 4326, remove = FALSE) # %>%
  # sf::st_transform("EPSG:2154")
})

test <- lapply(df_sf_list, function(sf_obj) {
  coords <- sf::st_coordinates(sf_obj)
  as.data.frame(coords)
})

# Variogramms
plot(gstat::variogram(x ~ 1, df_list[[1]], locations = ~ locx + locy))
plot(gstat::variogram(y ~ 1, df_list[[1]], locations = ~ locx + locy))
plot(gstat::variogram(Z ~ 1, df_list[[1]], locations = ~ locx + locy))

rm(A, B, crs1, epsilon, epsilon_list, fake.locations, Qu_s, Qu_xs, spde_s, spde_xs, test, testdf, testdfsf, u_s, u_s_s, u_s_xs, u_xs, us, ux, x, y, Z)

############# Apply models to simulated data #############
CI <- function(model, nbr = n_sim, beta = 3) {
  interval <- lapply(1:nbr, function(i) model[[i]]$summary.fixed[1, c("0.025quant", "0.975quant")])
  valid_intervals <- interval[!sapply(interval, is.null)]
  is_in_interval <- lapply(valid_intervals, function(x) {
    beta >= x["0.025quant"] && beta <= x["0.975quant"]
  })
  CI <- (sum(unlist(is_in_interval)) / length(valid_intervals)) * 100
  return(CI)
}
# load simulated data that contains a dataframe df_sf_list
load("data/simu_data.Rdata")

## SPDE
spde_spat <- inla.spde2.pcmatern(
  mesh = mesh.sim,
  prior.range = c(0.05, 0.05), # P(range < 0.05) = 0.05
  prior.sigma = c(3, 0.05) # P(sigma > 3) = 0.05
)

### NULL model ###
plan(multisession, workers = 4)
model_function_NULL <- function(i) {
  bru(y ~ 0 + x, family = "gaussian", data = df_sf_list[[i]])
}

model_NULL <- furrr::future_map(1:n_sim, model_function_NULL, .options = furrr_options(seed = 123))

save(model_NULL, file = "result_simu/result_model_Null.Rda")
# beta model NULL
bru_beta_NULL <- list()
bru_beta_NULL <- sapply(1:n_sim, function(i) model_NULL[[i]]$summary.fixed$mean)
bru_beta_NULL <- as.data.frame(bru_beta_NULL) %>% rename(value = "bru_beta_NULL")
bru_beta_NULL$Group <- rep(c("beta NULL"), each = n_sim)

# beta sd model NULL
bru_beta_sd_NULL <- list()
bru_beta_sd_NULL <- sapply(1:n_sim, function(i) model_NULL[[i]]$summary.fixed$sd)
bru_beta_sd_NULL <- as.data.frame(bru_beta_sd_NULL) %>% rename(value = "bru_beta_sd_NULL")
bru_beta_sd_NULL$Group <- rep(c("beta sd NULL"), each = n_sim)

# WAIC
WAIC_model_NULL <- c()
WAIC_model_NULL <- lapply(1:n_sim, function(i) model_NULL[[i]]$waic$waic) %>%
  unlist() %>%
  as.data.frame()

names(WAIC_model_NULL) <- "value"
WAIC_model_NULL$Group <- rep(c("WAIC NULL"), each = n_sim)

# DIC
DIC_model_NULL <- c()
DIC_model_NULL <- lapply(1:n_sim, function(i) model_NULL[[i]]$dic$dic) %>%
  unlist() %>%
  as.data.frame()

names(DIC_model_NULL) <- "value"
DIC_model_NULL$Group <- rep(c("DIC NULL"), each = n_sim)

# rm(model_function_NULL, model_NULL)

### Spatial model ###
plan(multisession, workers = 4)

model_function_spatial <- function(i) {
  bru(y ~ 0 + x + u_spa(geometry, model = spde_spat), family = "gaussian", data = df_sf_list[[i]])
}

model_function_spatial_gam <- function(i) {
  gam(y ~ -1 + x + s(locx, locy, k = 300, fx = TRUE), family = gaussian(), data = df_sf_list[[i]])
}

# Apply the function to each element in parallel
model_spatial <- furrr::future_map(1:n_sim, model_function_spatial_gam, .options = furrr_options(seed = 123))

summary(model_spatial[[1]])$se[1]
# save(model_spatial, file = "result_simu/result_model_Spatial.Rda")

# beta model spatial
bru_beta_spatial <- list()
bru_beta_spatial <- sapply(1:n_sim, function(i) model_spatial[[i]]$summary.fixed$mean)
# bru_beta_spatial <- sapply(1:n_sim, function(i) model_spatial[[i]]$coefficients["x"])
# bru_beta_spatial <- sapply(1:n_sim, function(i) summary(model_spatial[[i]])$se[1])
bru_beta_spatial <- as.data.frame(bru_beta_spatial) %>% rename(value = "bru_beta_spatial")
bru_beta_spatial$Group <- rep(c("beta spatial"), each = n_sim)

# beta sd model spatial
bru_beta_sd_spatial <- list()
bru_beta_sd_spatial <- sapply(1:n_sim, function(i) model_spatial[[i]]$summary.fixed$sd)
bru_beta_sd_spatial <- as.data.frame(bru_beta_sd_spatial) %>% rename(value = "bru_beta_sd_spatial")
bru_beta_sd_spatial$Group <- rep(c("beta sd spatial"), each = n_sim)

# WAIC
WAIC_model_spatial <- c()
WAIC_model_spatial <- lapply(1:n_sim, function(i) model_spatial[[i]]$waic$waic) %>%
  unlist() %>%
  as.data.frame()

names(WAIC_model_spatial) <- "value"
WAIC_model_spatial$Group <- rep(c("WAIC spatial"), each = n_sim)

# DIC
DIC_model_spatial <- c()
DIC_model_spatial <- lapply(1:n_sim, function(i) model_spatial[[i]]$dic$dic) %>%
  unlist() %>%
  as.data.frame()

names(DIC_model_spatial) <- "value"
DIC_model_spatial$Group <- rep(c("DIC spatial"), each = n_sim)

CI(model_spatial)
### RSR model ###
spde_RSR <- inla.spde2.pcmatern(
  mesh = mesh.sim,
  prior.range = c(15, .9999), # P(range < 15) = 0.9999
  prior.sigma = c(1.5, .0001) # P(sigma > 1.5) = 0.0001
)

plan(multisession, workers = 4)

model_function_RSR_formula <- function(i) {
  bru(
    components = ~ 0 + X_Effect(x) + u(geometry, model = spde_RSR),
    like(
      formula = y ~ 0 + X_Effect + u - x * sum(x * u) / sum(x^2),
      family = "gaussian",
      data = df_sf_list[[i]]
    ),
    options = list(
      control.compute = list(waic = TRUE, cpo = FALSE),
      control.inla = list(int.strategy = "eb"),
      verbose = FALSE
    )
  )
}

# start_time <- Sys.time()

plan(multisession, workers = 4)
model_RSR_formula <- furrr::future_map(
  1:n_sim,
  model_function_RSR_formula,
  .options = furrr_options(seed = TRUE)
)

# end_time <- Sys.time()
# time_taken <- end_time - start_time
# print(paste("Time taken:", time_taken))

# save(model_RSR_formula, file = "result_simu/result_model_RSR_formula_spde_spat.Rda")

# beta model RSR formula
bru_beta_RSR_formula <- list()
bru_beta_RSR_formula <- sapply(1:n_sim, function(i) model_RSR_formula[[i]]$summary.fixed$mean)
mean(unlist(bru_beta_RSR_formula))
bru_beta_RSR_formula <- as.data.frame(bru_beta_RSR_formula) %>% rename(value = "bru_beta_RSR_formula")
bru_beta_RSR_formula$Group <- rep(c("beta RSR formula"), each = n_sim)

# beta sd model RSR
bru_beta_sd_RSR_formula <- list()
bru_beta_sd_RSR_formula <- sapply(1:n_sim, function(i) model_RSR_formula[[i]]$summary.fixed$sd)
bru_beta_sd_RSR_formula <- as.data.frame(bru_beta_sd_RSR_formula) %>% rename(value = "bru_beta_sd_RSR_formula")
bru_beta_sd_RSR_formula$Group <- rep(c("beta sd RSR"), each = n_sim)

# WAIC
WAIC_model_RSR_formula <- c()
WAIC_model_RSR_formula <- lapply(1:n_sim, function(i) model_RSR_formula[[i]]$waic$waic) %>%
  unlist() %>%
  as.data.frame()

names(WAIC_model_RSR_formula) <- "value"
WAIC_model_RSR_formula$Group <- rep(c("WAIC RSR formula"), each = n_sim)

# DIC
DIC_model_RSR_formula <- c()
DIC_model_RSR_formula <- lapply(1:n_sim, function(i) model_RSR_formula[[i]]$dic$dic) %>%
  unlist() %>%
  as.data.frame()

names(DIC_model_RSR_formula) <- "value"
DIC_model_RSR_formula$Group <- rep(c("DIC RSR formula"), each = n_sim)

CI(model_RSR_formula)
# rm(model_function_RSR_formula, model_RSR_formula)

### RSR model with extraconstr ###

n_nodes <- mesh.sim$n
model_function_RSR_extra <- function(i) {
  print(i)
  A.rsr <- inla.spde.make.A(mesh.sim, loc = df_sf_list[[i]]$geometry)
  e.rsr <- rep(0, nrow(A.rsr))

  bru(
    components = ~ 0 + X_Effect(x) + u(geometry, model = spde_RSR, extraconstr = list(A = A.rsr, e = e.rsr)),
    like(
      formula = y ~ .,
      family = "gaussian",
      data = df_sf_list[[i]]
    ),
    options = list(
      control.compute = list(waic = TRUE, cpo = FALSE),
      control.inla = list(int.strategy = "eb"),
      verbose = FALSE
    )
  )
}

# start_time <- Sys.time()
plan(multisession, workers = 4)
model_RSR_extra <- furrr::future_map(
  1:n_sim,
  model_function_RSR_extra,
  .options = furrr_options(seed = 123)
)
# end_time <- Sys.time()
# time_taken <- end_time - start_time
# print(paste("Time taken:", time_taken))

# save(model_RSR_extra, file = "result_simu/result_model_RSR_extraconstr_spde_RSR.Rda")

# beta model RSR_extra
bru_beta_RSR_extra <- list()
bru_beta_RSR_extra <- sapply(1:n_sim, function(i) model_RSR_extra[[i]]$summary.fixed$mean)
bru_beta_RSR_extra <- as.data.frame(bru_beta_RSR_extra) %>% rename(value = "bru_beta_RSR_extra")
bru_beta_RSR_extra$Group <- rep(c("beta RSR_extra"), each = n_sim)

# beta sd model RSR_extra
bru_beta_sd_RSR_extra <- list()
bru_beta_sd_RSR_extra <- sapply(1:n_sim, function(i) model_RSR_extra[[i]]$dic$dic)
bru_beta_sd_RSR_extra <- as.data.frame(bru_beta_sd_RSR_extra) %>% rename(value = "bru_beta_sd_RSR_extra")
bru_beta_sd_RSR_extra$Group <- rep(c("beta sd RSR_extra"), each = n_sim)

# WAIC
WAIC_model_RSR_extra <- c()
WAIC_model_RSR_extra <- lapply(1:n_sim, function(i) model_RSR_extra[[i]]$waic$waic) %>%
  unlist() %>%
  as.data.frame()
names(WAIC_model_RSR_extra) <- "value"
WAIC_model_RSR_extra$Group <- rep(c("WAIC RSR_extra"), each = n_sim)

# DIC
DIC_model_RSR_extra <- c()
DIC_model_RSR_extra <- lapply(1:n_sim, function(i) model_RSR_extra[[i]]$dic$dic) %>%
  unlist() %>%
  as.data.frame()

names(DIC_model_RSR_extra) <- "value"
DIC_model_RSR_extra$Group <- rep(c("DIC RSR_extra"), each = n_sim)

CI(model_RSR_extra)
# rm(model_function_RSR_extra, model_RSR_extra)

### Spatial+ model ###
spde_X <- inla.spde2.pcmatern(
  mesh = mesh.sim,
  prior.range = c(0.01, 0.01),
  prior.sigma = c(0.15, 0.01)
)

plan(multisession, workers = 4)

x_fit <- lapply(df_sf_list, function(df) bru(x ~ 0 + u(geometry, model = spde_X), data = df))

x_true <- lapply(df_sf_list, function(i) i$x)

df_sf_list_new <- mapply(function(df, true, fit) {
  df$fitted_x_bru <- true - fit$summary.fitted.values[1:n, ]$mean
  return(df)
}, df_sf_list, x_true, x_fit) %>% as.data.frame()

model_function_spatial_plus <- function(i) {
  print(i)
  bru(
    y ~ 0 + fitted_x_bru +
      u(geometry, model = spde_spat),
    family = "gaussian", data = df_sf_list_new[[i]]
  )
}

model_spatial_plus <- furrr::future_map(1:n_sim, model_function_spatial_plus, .options = furrr_options(seed = 123))
save(model_spatial_plus, file = "result_simu/result_model_spatial_plus.Rda")

# beta model spatial plus
bru_beta_spatial_plus <- list()
bru_beta_spatial_plus <- sapply(1:n_sim, function(i) model_spatial_plus[[i]]$summary.fixed$mean)
bru_beta_spatial_plus <- as.data.frame(bru_beta_spatial_plus) %>% rename(value = "bru_beta_spatial_plus")
bru_beta_spatial_plus$Group <- rep(c("beta spatial plus"), each = n_sim)

# beta sd model spatial plus
bru_beta_sd_spatial_plus <- list()
bru_beta_sd_spatial_plus <- sapply(1:n_sim, function(i) model_spatial_plus[[i]]$summary.fixed$sd)
bru_beta_sd_spatial_plus <- as.data.frame(bru_beta_sd_spatial_plus) %>% rename(value = "bru_beta_sd_spatial_plus")
bru_beta_sd_spatial_plus$Group <- rep(c("beta sd spatial plus"), each = n_sim)

# WAIC
WAIC_model_spatial_plus <- c()
WAIC_model_spatial_plus <- lapply(1:n_sim, function(i) model_spatial_plus[[i]]$waic$waic) %>%
  unlist() %>%
  as.data.frame()

names(WAIC_model_spatial_plus) <- "value"
WAIC_model_spatial_plus$Group <- rep(c("WAIC spatial plus"), each = n_sim)

# DIC
DIC_model_spatial_plus <- c()
DIC_model_spatial_plus <- lapply(1:n_sim, function(i) model_spatial_plus[[i]]$dic$dic) %>%
  unlist() %>%
  as.data.frame()

names(DIC_model_spatial_plus) <- "value"
DIC_model_spatial_plus$Group <- rep(c("DIC spatial plus"), each = n_sim)

# rm(model_x_fit, x_fit, x_true, fitted_x_bru, df_sf_list_new, model_function_spatial_plus, model_spatial_plus)

### Spatial+ 2.0 ###

dist <- lapply(df_sf_list, function(df) sp::spDists(coordinates(df[, c("locx", "locy")])))

Q <- lapply(1:n_sim, function(i) inla.matern.cov(nu = 1 / 2, kappa = 3.1, dist[[i]]))
Q_inv <- lapply(Q, function(i) solve(i))
r <- lapply(Q_inv, function(i) eigen(i)) # spectral decomposition
eigen.vect <- lapply(r, function(i) i$vectors) # eigen vectors

# find the decomposition in eigen vector
coef <- lapply(1:n_sim, function(i) solve(eigen.vect[[i]], df_sf_list[[i]]$x))

fit_with_eigen_vectors <- function(nbr_eigen_vectors = 500,
                                   data = df_sf_list[[1]],
                                   decomposition = coef, e_vect = eigen.vect) {
  col_name <- paste("x.eigen", nbr_eigen_vectors, sep = "")
  if (nbr_eigen_vectors == 1) {
    data[, col_name] <- as.vector(e_vect[, 0:nbr_eigen_vectors] * decomposition[0:nbr_eigen_vectors])
  } else {
    data[, col_name] <- as.vector(e_vect[, 0:nbr_eigen_vectors] %*% decomposition[0:nbr_eigen_vectors])
  }
  res <- bru(
    components = as.formula(
      paste0("~ -1 + ",
        col_name, " + u(geometry, model = spde_spat)",
        sep = ""
      )
    ),
    like(
      formula = y ~ .,
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

model_function_spatial_plus2 <- function(i) {
  dist <- sp::spDists(sp::coordinates(df_sf_list[[i]]))

  Q <- inla.matern.cov(nu = 1 / 2, kappa = 3.1, dist)
  Q_inv <- solve(Q)
  r <- eigen(Q_inv) # spectral decomposition
  eigen.vect <- r$vectors

  # find the decomposition in eigen vector
  coef <- solve(eigen.vect, df_sf_list[[i]]$x)
  fit_with_eigen_vectors(nbr_eigen_vectors = 400, data = df_sf_list[[i]], decomposition = coef, e_vect = eigen.vect)
}
plan(multisession, workers = 4)
model_spatial_plus2 <- furrr::future_map(1:n_sim, model_function_spatial_plus2, .options = furrr_options(seed = 123))


model_spatial_plus2 <- lapply(1:n_sim, function(i) fit_with_eigen_vectors(nbr_eigen_vectors = 400, data = df_sf_list[[i]], decomposition = coef[[i]], e_vect = eigen.vect[[i]]))

# beta model spatial plus V2
bru_beta_spatial_plus2 <- list()
bru_beta_spatial_plus2 <- sapply(1:n_sim, function(i) model_spatial_plus2[[i]]$summary.fixed$mean)
bru_beta_spatial_plus2 <- as.data.frame(bru_beta_spatial_plus2) %>% rename(value = "bru_beta_spatial_plus2")
bru_beta_spatial_plus2$Group <- rep(c("beta spatial plus V2"), each = n_sim)

# beta sd model spatial plus V2
bru_beta_sd_spatial_plus2 <- list()
bru_beta_sd_spatial_plus2 <- sapply(1:n_sim, function(i) model_spatial_plus2[[i]]$summary.fixed$sd)
bru_beta_sd_spatial_plus2 <- as.data.frame(bru_beta_sd_spatial_plus2) %>% rename(value = "bru_beta_sd_spatial_plus2")
bru_beta_sd_spatial_plus2$Group <- rep(c("beta sd spatial plus V2"), each = n_sim)

# WAIC
WAIC_model_spatial_plus2 <- c()
WAIC_model_spatial_plus2 <- lapply(1:n_sim, function(i) model_spatial_plus2[[i]]$waic$waic) %>%
  unlist() %>%
  as.data.frame()

names(WAIC_model_spatial_plus2) <- "value"
WAIC_model_spatial_plus2$Group <- rep(c("WAIC spatial plus V2"), each = n_sim)

# DIC
DIC_model_spatial_plus2 <- c()
DIC_model_spatial_plus2 <- lapply(1:n_sim, function(i) model_spatial_plus2[[i]]$dic$dic) %>%
  unlist() %>%
  as.data.frame()

names(DIC_model_spatial_plus2) <- "value"
DIC_model_spatial_plus2$Group <- rep(c("DIC spatial plus V2"), each = n_sim)

rm(hyperparam, log_range_u, log_sigma_u, dist, Q, r, coef, model_spatial_plus2)

### gSEM ###
spde_spat <- inla.spde2.pcmatern(
  mesh = mesh.sim,
  prior.range = c(0.05, 0.05), # P(range < 0.05) = 0.05
  prior.sigma = c(3, 0.05) # P(sigma > 3) = 0.05
)

model_function_gSEM_bru <- function(i) {
  f_X_hat <- bru(x ~ u(geometry, model = spde_spat), data = df_sf_list[[i]], family = "gaussian")
  df_sf_list[[i]]$r_X <- df_sf_list[[i]]$x - f_X_hat$summary.fitted.values[1:n, ]$mean
  f_Y_hat <- bru(y ~ u(geometry, model = spde_spat), data = df_sf_list[[i]], family = "gaussian")
  df_sf_list[[i]]$r_Y <- df_sf_list[[i]]$y - f_Y_hat$summary.fitted.values[1:n, ]$mean
  lm(r_Y ~ 0 + r_X, data = df_sf_list[[i]])
}

# start_time <- Sys.time()
# plan(multisession, workers = 4)
model_gSEM_bru <- furrr::future_map(1:n_sim, model_function_gSEM_bru, .options = furrr_options(seed = 123))
# end_time <- Sys.time()

# time_taken <- end_time - start_time
# print(paste("Time taken:", time_taken))

model_gSEM[[1]]$summary.hyperpar$mean
bru_beta_gSEM <- sapply(1:n_sim, function(i) model_gSEM[[i]]$summary.fixed$mean)
bru_beta_gSEM <- as.data.frame(bru_beta_gSEM) %>% rename(value = "bru_beta_gSEM")
bru_beta_gSEM$Group <- rep(c("beta gSEM"), each = n_sim)
boxplot(bru_beta_gSEM$value)

bru_beta_gSEM_sd <- sapply(1:n_sim, function(i) model_gSEM[[i]]$summary.fixed$sd)
WAIC_model_gSEM <- sapply(1:n_sim, function(i) model_gSEM[[i]]$dic$dic)
DIC_model_gSEM <- sapply(1:n_sim, function(i) model_gSEM[[i]]$dic$dic)
range_u <- sapply(1:n_sim, function(i) model_gSEM[[i]]$summary$dic)
model_gSEM[[1]]$summary.
mean(WAIC_model_gSEM)
median(bru_beta_gSEM)
sd(bru_beta_gSEM)

CI(model_gSEM)
## Plot beta Standard error
combined_sd <- rbind(bru_beta_sd_NULL, bru_beta_sd_spatial, bru_beta_sd_RSR, bru_beta_sd_spatial_plus, bru_beta_sd_spatial_plus2)

ggplot(combined_sd, aes(x = "", y = value, fill = Group)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  labs(x = "", y = "Value", title = "Boxplot posterior beta sd") +
  theme_minimal()


## Plot WAIC
combined_WAIC <- rbind(WAIC_model_NULL, WAIC_model_spatial, WAIC_model_RSR_formula, WAIC_model_spatial_plus, WAIC_model_spatial_plus2)

ggplot(combined_WAIC, aes(x = "", y = value, fill = Group)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  labs(x = "", y = "Value", title = "Boxplot WAIC") +
  theme_minimal()


## Plot DIC
combined_DIC <- rbind(DIC_model_NULL, DIC_model_spatial, DIC_model_RSR_formula, DIC_model_spatial_plus, DIC_model_spatial_plus2)

ggplot(combined_DIC, aes(x = "", y = value, fill = Group)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  labs(x = "", y = "Value", title = "Boxplot DIC") +
  theme_minimal()


############################################## GAM #####################################################
model_fx <- TRUE

## NULL
model_GAM_function_NULL <- function(i) {
  mod <- lm(y ~ -1 + x, df_sf_list[[i]])
  bru_beta_hat <- mod$coefficients[1]
  return(mod)
}
plan(multisession, workers = 4)
model_GAM_NULL <- furrr::future_map(1:n_sim, model_GAM_function_NULL, .options = furrr_options(seed = TRUE))
r_X_values <- purrr::map_dbl(model_GAM_NULL, "x")
bru_beta_GAM_NULL <- r_X_values
bru_beta_GAM_NULL <- as.data.frame(bru_beta_GAM_NULL) %>% rename(value = "bru_beta_GAM_NULL")
bru_beta_GAM_NULL$Group <- rep(c("GAM NULL"), each = n_sim)


## spatial
model_GAM_function_spatial <- function(i) {
  mod <- gam(y ~ -1 + x + s(locx, locy, k = 300, fx = model_fx, sp = 1), data = df_sf_list[[i]], method = "GCV.Cp")
  bru_beta_hat <- mod$coefficients[1]
  return(mod)
}
plan(multisession, workers = 4)
model_GAM_spatial <- furrr::future_map(1:n_sim, model_GAM_function_spatial, .options = furrr_options(seed = TRUE))
r_X_values <- purrr::map_dbl(model_GAM_spatial, "x")
bru_beta_GAM_spatial <- r_X_values
bru_beta_GAM_spatial <- as.data.frame(bru_beta_GAM_spatial) %>% rename(value = "bru_beta_GAM_spatial")
bru_beta_GAM_spatial$Group <- rep(c("GAM spatial"), each = n_sim)


## RSR
model_GAM_function_RSR <- function(i) {
  mod_list <- gam(y ~ x + s(locx, locy, k = 300, fx = model_fx, sp = 1), data = df_sf_list[[i]], fit = FALSE)
  B_sp <- mod_list$X[, -2]
  x <- df_sf_list[[i]]$x
  P <- 1 / sum(x^2) * x %*% t(x)
  B_sp_tilde <- (diag(x = 1, nrow = n, ncol = n) - P) %*% B_sp
  mod_list$X[, -2] <- B_sp_tilde
  mod <- gam(G = mod_list, method = "GCV.Cp")
  bru_beta_hat <- mod$coefficients[2]
  return(mod)
}
plan(multisession, workers = 4)
model_GAM_RSR <- furrr::future_map(1:n_sim, model_GAM_function_RSR, .options = furrr_options(seed = TRUE))
r_X_values <- purrr::map_dbl(model_GAM_RSR, "x")
bru_beta_GAM_RSR <- r_X_values
bru_beta_GAM_RSR <- as.data.frame(bru_beta_GAM_RSR) %>% rename(value = "bru_beta_GAM_RSR")
bru_beta_GAM_RSR$Group <- rep(c("GAM RSR"), each = n_sim)

## gSEM
model_GAM_function_gSEM <- function(i) {
  f_X_hat <- gam(x ~ -1 + s(locx, locy, k = 300, fx = model_fx, sp = 1), data = df_sf_list[[i]], method = "GCV.Cp")$fitted.values
  r_X <- df_sf_list[[i]]$x - f_X_hat
  f_Y_hat <- gam(y ~ -1 + s(locx, locy, k = 300, fx = model_fx, sp = 1), data = df_sf_list[[i]], method = "GCV.Cp")$fitted.values
  r_Y <- df_sf_list[[i]]$y - f_Y_hat
  mod <- lm(r_Y ~ -1 + r_X)
  bru_beta_hat <- mod$coefficients[1]
  return(mod)
}
plan(multisession, workers = 4)
model_GAM_gSEM <- furrr::future_map(1:n_sim, model_GAM_function_gSEM, .options = furrr_options(seed = TRUE))
r_X_values <- purrr::map_dbl(model_GAM_gSEM, "r_X")
bru_beta_GAM_gSEM <- r_X_values
bru_beta_GAM_gSEM <- as.data.frame(bru_beta_GAM_gSEM) %>% rename(value = "bru_beta_GAM_gSEM")
bru_beta_GAM_gSEM$Group <- rep(c("GAM gSEM"), each = n_sim)


## GAM spatial+
model_function_GAM_spatialplus <- function(i) {
  f_X_hat <- gam(x ~ -1 + s(locx, locy, k = 300, fx = model_fx, sp = 1), data = df_sf_list[[i]], method = "GCV.Cp")$fitted.values
  r_X <- df_sf_list[[i]]$x - f_X_hat
  mod <- gam(y ~ -1 + r_X + s(locx, locy, k = 300, fx = model_fx, sp = 1), data = df_sf_list[[i]], method = "GCV.Cp")
  bru_beta_hat <- mod$coefficients[1]
  return(mod)
}
plan(multisession, workers = 4)
model_GAM_spatialplus <- furrr::future_map(1:n_sim, model_function_GAM_spatialplus, .options = furrr_options(seed = TRUE))
r_X_values <- purrr::map_dbl(model_GAM_spatialplus, "r_X")
bru_beta_GAM_spatialplus <- r_X_values
bru_beta_GAM_spatialplus <- as.data.frame(bru_beta_GAM_spatialplus) %>% rename(value = "bru_beta_GAM_spatialplus")
bru_beta_GAM_spatialplus$Group <- rep(c("GAM spatial+"), each = n_sim)


## Plot beta bru + GAM

combined <- rbind(bru_beta_NULL, bru_beta_spatial, bru_beta_RSR, bru_beta_RSR_extra, bru_beta_gSEM, bru_beta_spatial_plus, bru_beta_spatial_plus2, bru_beta_GAM_spatial, bru_beta_GAM_RSR, bru_beta_GAM_gSEM, bru_beta_GAM_spatialplus)
desired_order <- c("beta Null", "beta Spatial", "beta RSR extra", "beta RSR formula", "beta gSEM", "beta spatial plus", "beta spatial plus V2", "GAM spatial", "GAM RSR", "GAM gSEM", "GAM spatial+")

combined$Group <- factor(combined$Group, levels = desired_order)

ggplot(combined, aes(x = Group, y = value, fill = Group)) +
  geom_blank(data = subset(combined_GAM, Group == "beta SVC")) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  # scale_fill_discrete(breaks=c('beta NULL', 'beta RSR i', 'beta spatial', "beta spatial plus", "beta spatial plus V2", "beta SVC plus", "beta SVC")) +
  labs(x = "", y = "", title = expression(hat(beta))) +
  geom_hline(yintercept = beta, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 7.5, linetype = "dashed", color = "black") +
  # guides(fill = guide_legend(override.aes = list(alpha = c(1, 1, 1, 1, 1, 1, 0)))) +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

load("result_simu/result_model_RSR_extraconstr_spde_RSR.Rda")
load("result_simu/result_model_RSR_extraconstr_spde_spat.Rda")
load("result_simu/result_model_RSR_formula_spde_RSR.Rda")

model_RSR_extra_RSR <- model_RSR_extra
model_RSR_extra_spat <- model_RSR_extra
model_RSR_formula_RSR <- model_RSR

get_beta <- function(model, name = "beta") {
  df_beta <- data.frame(matrix(ncol = 2, nrow = n_sim))
  df_beta$X1 <- sapply(1:n_sim, function(i) model[[i]]$summary.fixed$mean)
  df_beta$X2 <- rep(c(name), each = n_sim)
  colnames(df_beta) <- c("value", "Group")
  return(df_beta)
}

beta_RSR_extra_RSR <- get_beta(model_RSR_extra_RSR, "beta RSR extra RSR")
beta_RSR_extra_spat <- get_beta(model_RSR_extra_spat, "beta RSR extra spat")
beta_RSR_formula_RSR <- get_beta(model_RSR_formula_RSR, "beta RSR formula RSR")

combined <- rbind(beta_RSR_extra_RSR, beta_RSR_extra_spat, beta_RSR_formula_RSR)

ggplot(combined, aes(x = "", y = value, fill = Group)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  labs(x = "", y = "Value", title = "Boxplot posterior beta") +
  geom_hline(yintercept = beta, linetype = "dashed", color = "red") +
  theme_minimal()

sd(beta_RSR_extra_RSR$value)
sd(beta_RSR_extra_spat$value)
sd(beta_RSR_formula_RSR$value)
