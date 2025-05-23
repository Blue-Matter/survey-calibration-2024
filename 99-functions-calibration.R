

fit_ghl <- function(dat, formula = num ~ 1 + vessel, ...) {

  require(sdmTMB)

  mesh <- sdmTMB::make_mesh(
    dat,
    c("X", "Y"),
    type = "cutoff",
    cutoff = 30
  )
  #plot(mesh)

  all_time <- seq(min(dat$year), max(dat$year))
  extra_time <- setdiff(all_time, unique(dat$year))

  fit <- sdmTMB::sdmTMB(
    formula,
    dat,
    mesh = mesh,
    time = "year",
    family = nbinom2(),
    spatiotemporal = "rw",
    extra_time = extra_time,
    anisotropy = TRUE,
    ...
  )
  return(fit)

}



make_index <- function(fit, grid, make_index = TRUE) {
  year <- seq(min(fit$data$year), max(fit$data$year))

  require(sdmTMB)

  newdata <- sdmTMB::replicate_df(grid, "year", year)

  pred <- predict(
    fit,
    newdata = newdata,
    return_tmb_object = TRUE
  )

  if (make_index) {
    index <- get_index(pred, bias_correct = TRUE, area = newdata$area)
    return(index)
  } else {
    return(pred)
  }
}


