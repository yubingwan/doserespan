if (!require("drc", character.only = TRUE)) {
    install.packages("drc", dependencies = TRUE)
    library("drc", character.only = TRUE)
}

if (!require("broom", character.only = TRUE)) {
    install.packages("broom", dependencies = TRUE)
    library("broom", character.only = TRUE)
}

if (!require("reticulate", character.only = TRUE)) {
    install.packages("reticulate", dependencies = TRUE)
    library("reticulate", character.only = TRUE)
}

theme_set(theme_bw())

fit_4pl_model <- function(
    df, ec50 = 1e-6, h = -1, top = 42, bottom = 0,
    constr_top = NULL, constr_bottom = NULL, constr_h = NULL, ...) {

  if(!is.null(constr_top) & is.null(constr_bottom)) {
    top <- constr_top
    start_list <- list(ec50 = ec50, h = h, bottom = bottom)
    f1 <- response ~ bottom + (top - bottom)/(1 + 10^((log10(ec50) - log10(dose)) * h))
    m1 <- nls(f1, data = df, start = start_list, ...)
    coeffs <- list(ec50 = coef(m1)[1], h = coef(m1)[2], top = top, bottom = coef(m1)[3])
  }

  if(is.null(constr_top) & !is.null(constr_bottom)) {
    bottom <- constr_bottom
    f2 <- response ~ bottom + (top - bottom)/(1 + 10^((log10(ec50) - log10(dose)) * h))
    m1 <- nls(f2, data = df, start = list(ec50 = ec50, h = h, top = top), ...)
    coeffs <- list(ec50 = coef(m1)[1], h = coef(m1)[2], top = coef(m1)[3], bottom = bottom)
  }

  if(!is.null(constr_top) & !is.null(constr_bottom)) {
    top <- constr_top
    bottom <- constr_bottom
    f3 <- response ~ bottom + (top - bottom)/(1 + 10^((log10(ec50) - log10(dose)) * h))
    m1 <- nls(f3, data = df, start = list(ec50 = ec50, h = h), ...)
    coeffs <- list(ec50 = coef(m1)[1], h = coef(m1)[2], top = top, bottom  = bottom)
  }

  if(is.null(constr_top) & is.null(constr_bottom)) {
    f4 <- response ~ bottom + (top - bottom)/(1 + 10^((log10(ec50) - log10(dose)) * h))
    m1 <- nls(f4, data = df, start = list(ec50 = ec50, h = h, top = top, bottom = bottom), ...)
    coeffs <- list(ec50 = coef(m1)[1], h = coef(m1)[2], top = coef(m1)[3], bottom  = coef(m1)[4])
  }

  if(!is.null(constr_h)) {
    h <- constr_h
    f5 <- response ~ bottom + (top - bottom)/(1 + 10^((log10(ec50) - log10(dose)) * h))
    m1 <- nls(f5, data = df, start = list(ec50 = ec50, top = top, bottom = bottom), ...)
    coeffs <- list(ec50 = coef(m1)[1], h = h, top = coef(m1)[2], bottom  = coef(m1)[3])
  }

  coeffs
}

plot_dose_response_curve <- function(c0, df, ...) {
  plot(x = df$dose, y = df$response, log = "x", xlab = "Dose (M)", ylab = "Response", ...)
  x <- 10^seq(-9, -2, 0.1)
  y <- c0$bottom + (c0$top - c0$bottom)/(1 + 10^((log10(c0$ec50) - log10(x)) * c0$h))
  lines(x, y, col="blue", lwd=2, ...)
}

plot_dose_response_with_fit <- function(
    data, coef, title = "Dose-Response Curve", color = "blue", width = 8, height = 6) {

  data <- data %>%
    mutate(fitted = coef$bottom + (coef$top - coef$bottom)/(1 + 10^((log10(coef$ec50) - log10(dose)) * coef$h)))

  p <- ggplot(data, aes(x = dose)) +
    geom_point(aes(y = response), color = color, size = 3, show.legend = TRUE) +
    geom_line(aes(y = fitted), color = color, show.legend = TRUE) +
    ggtitle(title) +
    xlab("Dose") +
    ylab("Response") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("Observed Data" = color, "Fitted Curve" = color))

  # Print the plot
  print(p)
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
    signal <- top / (1 + (dose/ec50)^h)
  } else { # stimulation
    signal <- bottom + (top - bottom) * (dose^h / (ec50^h + dose^h))
  }

  # Add random noise
  response <- signal + rnorm(n, mean = 0, sd = noiseLevel * top)

  # Return a data frame
  return(data.frame(dose = dose, response = response, signal = signal))
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
