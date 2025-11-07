# Symuliatsii Monte-Karlo dlia otsinky efektyvnosti PMM2
# Chastyna 1: Modeliuvannia ta porivniannia metodiv na riznykh rozpodilakh pokhybok

# Perevirka naiavnosti ta pidkliuchennia neobkhidnykh pakunkiv
required_pkgs <- c("EstemPMM", "ggplot2", "gridExtra", "dplyr", "parallel")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Dlia tsoho demo vstanovit pakunky: ",
       paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

library(EstemPMM)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(parallel)

#############################################################
# Dopomizhni funktsii dlia symuliatsii Monte-Karlo
#############################################################

# Funktsiia dlia heneratsii danykh z riznymy rozpodilamy pomylok
generate_data <- function(n, distribution, a0, a1, ...) {
  x <- rnorm(n, mean = 0, sd = 1)

  # Heneratsiia pomylok z zadanym rozpodilom
  errors <- switch(distribution,
                   "gaussian" = rnorm(n, mean = 0, sd = 1),
                   "t" = rt(n, df = 3),
                   "gamma" = {
                     # Parametryzovana hamma dlia nulovoi serednoi
                     alpha <- 2
                     beta <- 1/sqrt(alpha)
                     rgamma(n, shape = alpha, scale = beta) - alpha*beta
                   },
                   "exponential" = {
                     # Eksponentsiina zi zsuvom dlia nulovoi serednoi
                     lambda <- 1
                     rexp(n, rate = lambda) - 1/lambda
                   },
                   "chi-squared" = {
                     # Khi-kvadrat zi zsuvom dlia nulovoi serednoi
                     df <- 3
                     rchisq(n, df = df) - df
                   },
                   "lognormal" = {
                     # Lohnormalna zi zsuvom dlia nulovoi serednoi
                     sigma <- 0.5
                     exp(rnorm(n, mean = -sigma^2/2, sd = sigma)) - 1
                   })

  # Obchyslennia znachen y
  y <- a0 + a1 * x + errors

  return(data.frame(x = x, y = y, errors = errors))
}

# Funktsiia dlia porivniannia metodiv PMM2 ta OLS
compare_methods <- function(data, true_a0, true_a1) {
  # Pidhonka OLS
  ols_fit <- lm(y ~ x, data = data)

  # Pidhonka PMM2
  pmm2_fit <- lm_pmm2(y ~ x, data = data, verbose = FALSE)

  # Obchyslennia zalyshkiv
  ols_resid <- residuals(ols_fit)
  pmm2_resid <- pmm2_fit@residuals

  # Obchyslennia MSE
  ols_mse <- mean(ols_resid^2)
  pmm2_mse <- mean(pmm2_resid^2)

  # Obchyslennia AIC
  ols_aic <- AIC(ols_fit)
  pmm2_aic <- AIC(pmm2_fit)

  # Obchyslennia zmishchennia otsinok
  ols_bias_a0 <- coef(ols_fit)[1] - true_a0
  ols_bias_a1 <- coef(ols_fit)[2] - true_a1

  pmm2_bias_a0 <- pmm2_fit@coefficients[1] - true_a0
  pmm2_bias_a1 <- pmm2_fit@coefficients[2] - true_a1

  # Momenty rozpodilu pomylok
  moments <- EstemPMM::compute_moments(data$errors)

  return(list(
    ols_coef = coef(ols_fit),
    pmm2_coef = pmm2_fit@coefficients,
    ols_mse = ols_mse,
    pmm2_mse = pmm2_mse,
    ols_aic = ols_aic,
    pmm2_aic = pmm2_aic,
    ols_bias = c(a0 = ols_bias_a0, a1 = ols_bias_a1),
    pmm2_bias = c(a0 = pmm2_bias_a0, a1 = pmm2_bias_a1),
    mse_ratio = pmm2_mse / ols_mse,
    moments = moments
  ))
}

# Funktsiia dlia provedennia symuliatsii Monte-Karlo
monte_carlo <- function(n_sim, n_samples, distribution, true_a0, true_a1, parallel = FALSE) {

  run_sim <- function(i) {
    data <- generate_data(n_samples, distribution, true_a0, true_a1)
    results <- compare_methods(data, true_a0, true_a1)

    # Zberehty tilky vazhlyvi rezultaty dlia ekonomii pam'iati
    return(list(
      ols_coef = results$ols_coef,
      pmm2_coef = results$pmm2_coef,
      ols_mse = results$ols_mse,
      pmm2_mse = results$pmm2_mse,
      mse_ratio = results$mse_ratio
    ))
  }

  # Vykorystannia paralelnykh obchyslen, iakshcho vkazano
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- parallel::detectCores() - 1
    results <- parallel::mclapply(1:n_sim, function(i) run_sim(i), mc.cores = n_cores)
  } else {
    # Initsializatsiia prohres-baru
    pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

    results <- vector("list", n_sim)
    for (i in 1:n_sim) {
      results[[i]] <- run_sim(i)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  # Obchyslennia statystyk po rezultatam
  ols_a0 <- sapply(results, function(x) x$ols_coef[1])
  ols_a1 <- sapply(results, function(x) x$ols_coef[2])

  pmm2_a0 <- sapply(results, function(x) x$pmm2_coef[1])
  pmm2_a1 <- sapply(results, function(x) x$pmm2_coef[2])

  ols_mse <- sapply(results, function(x) x$ols_mse)
  pmm2_mse <- sapply(results, function(x) x$pmm2_mse)

  mse_ratio <- sapply(results, function(x) x$mse_ratio)

  # Heneratsiia odnorazovoho naboru danykh dlia obchyslennia momentiv
  data <- generate_data(10000, distribution, true_a0, true_a1)
  moments <- EstemPMM::compute_moments(data$errors)

  return(list(
    ols_a0_mean = mean(ols_a0),
    ols_a0_sd = sd(ols_a0),
    ols_a1_mean = mean(ols_a1),
    ols_a1_sd = sd(ols_a1),

    pmm2_a0_mean = mean(pmm2_a0),
    pmm2_a0_sd = sd(pmm2_a0),
    pmm2_a1_mean = mean(pmm2_a1),
    pmm2_a1_sd = sd(pmm2_a1),

    ols_mse_mean = mean(ols_mse),
    pmm2_mse_mean = mean(pmm2_mse),

    mse_ratio_mean = mean(mse_ratio),

    theoretical_g = moments$g,
    c3 = moments$c3,
    c4 = moments$c4,

    ols_a0 = ols_a0,
    ols_a1 = ols_a1,
    pmm2_a0 = pmm2_a0,
    pmm2_a1 = pmm2_a1
  ))
}

# Funktsiia dlia vizualizatsii rezultativ Monte-Karlo
plot_monte_carlo_results <- function(mc_results, distribution_name, true_a0, true_a1) {
  # Rozpodil otsinok a0
  p1 <- ggplot() +
    geom_density(aes(x = mc_results$ols_a0, fill = "OLS"), alpha = 0.5) +
    geom_density(aes(x = mc_results$pmm2_a0, fill = "PMM2"), alpha = 0.5) +
    geom_vline(xintercept = true_a0, linetype = "dashed") +
    labs(title = paste("Rozpodil otsinok a0 -", distribution_name),
         x = "a0", y = "Hustyna") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                      name = "Metod") +
    theme_minimal()

  # Rozpodil otsinok a1
  p2 <- ggplot() +
    geom_density(aes(x = mc_results$ols_a1, fill = "OLS"), alpha = 0.5) +
    geom_density(aes(x = mc_results$pmm2_a1, fill = "PMM2"), alpha = 0.5) +
    geom_vline(xintercept = true_a1, linetype = "dashed") +
    labs(title = paste("Rozpodil otsinok a1 -", distribution_name),
         x = "a1", y = "Hustyna") +
    scale_fill_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                      name = "Metod") +
    theme_minimal()

  # QQ-hrafiky
  p3 <- ggplot() +
    geom_qq(aes(sample = mc_results$ols_a0 - true_a0, color = "OLS")) +
    geom_qq(aes(sample = mc_results$pmm2_a0 - true_a0, color = "PMM2")) +
    geom_qq_line(aes(sample = mc_results$ols_a0 - true_a0)) +
    labs(title = paste("QQ-hrafik dlia a0 -", distribution_name),
         x = "Teoretychni kvantyli", y = "Vybirkovi kvantyli") +
    scale_color_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                       name = "Metod") +
    theme_minimal()

  p4 <- ggplot() +
    geom_qq(aes(sample = mc_results$ols_a1 - true_a1, color = "OLS")) +
    geom_qq(aes(sample = mc_results$pmm2_a1 - true_a1, color = "PMM2")) +
    geom_qq_line(aes(sample = mc_results$ols_a1 - true_a1)) +
    labs(title = paste("QQ-hrafik dlia a1 -", distribution_name),
         x = "Teoretychni kvantyli", y = "Vybirkovi kvantyli") +
    scale_color_manual(values = c("OLS" = "blue", "PMM2" = "red"),
                       name = "Metod") +
    theme_minimal()

  # Vidobrazhennia statystyk
  stats_text <- paste(
    paste("Asymetriia (c3):", round(mc_results$c3, 4)),
    paste("Ekstses (c4):", round(mc_results$c4, 4)),
    paste("Teoretychnyi koefitsiient g:", round(mc_results$theoretical_g, 4)),
    paste("Faktychnyi koefitsiient (MSE):", round(mc_results$mse_ratio_mean, 4)),
    paste("Serednie znachennia a0 (OLS):", round(mc_results$ols_a0_mean, 4),
          "+/-", round(mc_results$ols_a0_sd, 4)),
    paste("Serednie znachennia a0 (PMM2):", round(mc_results$pmm2_a0_mean, 4),
          "+/-", round(mc_results$pmm2_a0_sd, 4)),
    paste("Serednie znachennia a1 (OLS):", round(mc_results$ols_a1_mean, 4),
          "+/-", round(mc_results$ols_a1_sd, 4)),
    paste("Serednie znachennia a1 (PMM2):", round(mc_results$pmm2_a1_mean, 4),
          "+/-", round(mc_results$pmm2_a1_sd, 4)),
    sep = "\n"
  )

  p5 <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = stats_text, hjust = 0) +
    theme_void() +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = paste("Statystyka -", distribution_name))

  # Ob'iednannia hrafikiv
  grid.arrange(p1, p2, p3, p4, p5, ncol = 2)
}

#############################################################
# Nalashtuvannia ta zapusk symuliatsii
#############################################################

set.seed(42)

# Parametry symuliatsii
n_sim <- 1000       # Kilkist symuliatsii Monte-Karlo
n_samples <- 100    # Rozmir vybirky v kozhnii symuliatsii
true_a0 <- 2        # Spravzhnie znachennia a0
true_a1 <- 1.5      # Spravzhnie znachennia a1

# Dlia shvydkoho demonstratsiinoho zapusku mozhna vykorystaty:
# n_sim <- 100      # Mensha kilkist symuliatsii dlia shvydshoho zapusku
# n_samples <- 50   # Menshyi rozmir vybirky

# Spysok rozpodiliv dlia testuvannia
distributions <- c("gaussian", "t", "gamma", "exponential", "chi-squared", "lognormal")
distribution_names <- c("Normalnyi", "Stiudenta (df=3)", "Hamma (a=2)",
                        "Eksponentsiinyi", "Khi-kvadrat (df=3)", "Lohnormalnyi")

# Funktsiia dlia vykonannia vsikh symuliatsii
run_all_simulations <- function() {
  results <- list()

  for(i in 1:length(distributions)) {
    cat("\nVykonuietsia symuliatsiia dlia rozpodilu:", distribution_names[i], "\n")

    mc_results <- monte_carlo(n_sim, n_samples,
                              distributions[i],
                              true_a0, true_a1,
                              parallel = TRUE)

    results[[distributions[i]]] <- mc_results

    # Vizualizatsiia rezultativ
    plot_monte_carlo_results(mc_results, distribution_names[i], true_a0, true_a1)

    # Vyvid rezultativ u konsol
    cat("\nRezultaty dlia rozpodilu:", distribution_names[i], "\n")
    cat("Asymetriia (c3):", round(mc_results$c3, 4), "\n")
    cat("Ekstses (c4):", round(mc_results$c4, 4), "\n")
    cat("Teoretychnyi koefitsiient g:", round(mc_results$theoretical_g, 4), "\n")
    cat("Faktychnyi koefitsiient (MSE):", round(mc_results$mse_ratio_mean, 4), "\n")
    cat("Pokrashchennia efektyvnosti:",
        round((1 - mc_results$mse_ratio_mean) * 100, 2), "%\n")

    cat("a0 (OLS):", round(mc_results$ols_a0_mean, 4),
        "+/-", round(mc_results$ols_a0_sd, 4), "\n")
    cat("a0 (PMM2):", round(mc_results$pmm2_a0_mean, 4),
        "+/-", round(mc_results$pmm2_a0_sd, 4), "\n")

    cat("a1 (OLS):", round(mc_results$ols_a1_mean, 4),
        "+/-", round(mc_results$ols_a1_sd, 4), "\n")
    cat("a1 (PMM2):", round(mc_results$pmm2_a1_mean, 4),
        "+/-", round(mc_results$pmm2_a1_sd, 4), "\n")
  }

  # Stvorennia pidsumkovoho porivniannia vsikh rozpodiliv
  summary_table <- data.frame(
    Distribution = distribution_names,
    Skewness = sapply(results, function(x) round(x$c3, 4)),
    Kurtosis = sapply(results, function(x) round(x$c4, 4)),
    Theoretical_g = sapply(results, function(x) round(x$theoretical_g, 4)),
    Actual_g = sapply(results, function(x) round(x$mse_ratio_mean, 4)),
    Improvement = sapply(results, function(x)
      round((1 - x$mse_ratio_mean) * 100, 2))
  )

  cat("\n\nPidsumok vsikh symuliatsii:\n")
  print(summary_table)

  # Vizualizatsiia porivniannia efektyvnosti
  p <- ggplot(summary_table, aes(x = reorder(Distribution, -Improvement))) +
    geom_bar(aes(y = Improvement, fill = "Faktychne"), stat = "identity",
             alpha = 0.7, position = position_dodge()) +
    geom_point(aes(y = (1 - Theoretical_g) * 100, color = "Teoretychne"),
               size = 3) +
    labs(title = "Porivniannia efektyvnosti PMM2 vidnosno OLS",
         x = "Rozpodil pomylok",
         y = "Pokrashchennia efektyvnosti (%)") +
    scale_fill_manual(values = c("Faktychne" = "steelblue"),
                      name = "Pokrashchennia") +
    scale_color_manual(values = c("Teoretychne" = "red"),
                       name = "Pokrashchennia") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p)

  return(list(results = results, summary = summary_table))
}

# Shchob zapustyty vsi symuliatsii, vyklychte vruchnu:
# run_all_simulations()
