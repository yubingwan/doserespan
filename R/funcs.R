fit_4pl_model <- function(df, log50_0 = -6, h_0 = -1, top_0 = 42, bottom_0 = 0, constr_top = NULL, constr_bottom = NULL, ...) {

  if(!is.null(constr_top) & is.null(constr_bottom)) {
    top <- constr_top
    f1 <- response ~ bottom + (top - bottom)/(1 + 10^((logEC50 - log10(dose)) * h))
    m1 <- nls(f1, data = df, start = list(logEC50 = log50_0, h = h_0, bottom = bottom_0), ...)
    coeffs <- list(logEC50 = coef(m1)[1], h = coef(m1)[2], top = top, bottom = coef(m1)[3])
  }

  if(is.null(constr_top) & !is.null(constr_bottom)) {
    bottom <- constr_bottom
    f2 <- response ~ bottom + (top - bottom)/(1 + 10^((logEC50 - log10(dose)) * h))
    m1 <- nls(f2, data = df, start = list(logEC50 = log50_0, h = h_0, top = top_0), ...)
    coeffs <- list(logEC50 = coef(m1)[1], h = coef(m1)[2], top = coef(m1)[3], bottom = bottom)
  }

  if(!is.null(constr_top) & !is.null(constr_bottom)) {
    top <- constr_top
    bottom <- constr_bottom
    f3 <- response ~ bottom + (top - bottom)/(1 + 10^((logEC50 - log10(dose)) * h))
    m1 <- nls(f3, data = df, start = list(logEC50 = log50_0, h = h_0), ...)
    coeffs <- list(logEC50 = coef(m1)[1], h = coef(m1)[2], top = top, bottom  = bottom)
  }

  if(is.null(constr_top) & is.null(constr_bottom)) {
    f4 <- response ~ bottom + (top - bottom)/(1 + 10^((logEC50 - log10(dose)) * h))
    m1 <- nls(f4, data = df, start = list(logEC50 = log50_0, h = h_0, top = top_0, bottom = bottom_0), ...)
    coeffs <- list(logEC50 = coef(m1)[1], h = coef(m1)[2], top = coef(m1)[3], bottom  = coef(m1)[4])
  }

  coeffs
}

plot_dose_response_curve <- function(c0, df, ...) {
  plot(x = df$dose, y = df$response, log = "x", xlab = "Dose (M)", ylab = "Response", ...)
  X1 <- seq(-9, -2, 0.1) # X-data for the curve
  y <- c0$bottom + (c0$top - c0$bottom)/(1 + 10^((c0$logEC50 - X1) * c0$h)) # Y data
  lines(10^X1, y, col="blue", lwd=2, ...) # plot the curve
}

simulate_dose_response_data <- function(
    n = 100, type = c("inhibition", "stimulation"),
    ec50 = 50, h = 1, top = 100, bottom = 0, noiseLevel = 0.05) {

  # Ensure valid type is selected
  type <- match.arg(type)

  # Generate dose values (log-scale can be more realistic for certain scenarios)
  dose <- seq(from = bottom, to = top, length.out = n)

  # Generate response values
  if (type == "inhibition") {
    response <- top / (1 + (dose/ec50)^h)
  } else { # stimulation
    response <- bottom + (top - bottom) * (dose^h / (ec50^h + dose^h))
  }

  # Add random noise
  response <- response + rnorm(n, mean = 0, sd = noiseLevel * top)

  # Return a data frame
  return(data.frame(dose = dose, response = response))
}

## Analysis by "drc" R package - (under development)
dose_response_analysis <- function(data, response_type = "inhibition", constrained = FALSE) {
  # Choose the model based on response type and constraints
  if (response_type == "inhibition") {
    model <- if (constrained) LL.4() else LL.3()
  } else if (response_type == "stimulation") {
    model <- if (constrained) LL.5() else LL.2()
  } else {
    stop("Invalid response type")
  }

  # Fit the model
  fit <- drm(response ~ dose, data = data, fct = model)

  return(fit)
}

## 2nd approach generating simulation data
# simulate_dose_response_data <- function(n = 10, response_type = "inhibition", variability = 0.1) {
#   set.seed(123) # for reproducibility
#
#   dose <- seq(0, 10, length.out = n)
#   if (response_type == "inhibition") {
#     response <- 100 / (1 + exp(5 - 0.5 * dose)) + rnorm(n, sd = variability * 100)
#   } else if (response_type == "stimulation") {
#     response <- 100 * (1 - exp(-0.5*dose)) + rnorm(n, sd = variability * 100)
#   } else {
#     stop("Invalid response type")
#   }
#
#   data.frame(dose = dose, response = response)
# }
