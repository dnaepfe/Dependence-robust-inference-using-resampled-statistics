# if needed:
#install.packages("clubSandwich") 

library(clubSandwich)

#X = Data matrix, mu = Null hypothesis, Rn = numbers of resamples
Mean_stat <- function(X,mu,Rn){
  n <- length(X)
  
  sigma_hat <- var(X) * (n-1)/n
  
  indices <- sample(1:n,Rn,replace=TRUE)
  W <- (X[indices] - mu) / sqrt(sigma_hat)
  
  T_M_tilde <- sum(W) / sqrt(Rn)
  T_M <- T_M_tilde^2
  
  return(T_M)
}

#X = Data, mu = Null hypothesis, Rn = numbers of resamples  
U_stat <- function(X,mu,Rn){
  n <- length(X)
  
  sigma_hat <- var(X) * (n-1)/n
  
  i <- sample.int(n, Rn, replace = TRUE)
  j <- sample.int(n, Rn, replace = TRUE)
  
  duplicates <- which(i == j)
  while (length(duplicates)) {
    j[duplicates] <- sample.int(n, length(duplicates), replace = TRUE)
    duplicates <- duplicates[i[duplicates] == j[duplicates]]
  }
  
  W <- (X[i] - mu) * (X[j] - mu) / sigma_hat
  T_U <- sum(W) / sqrt(Rn)
  
  return(T_U)
}

alpha <- 0.05
iterations <- 2000
beta_null <- c(3,1.5)  #(intercept,slope)
beta_alt <- c(3.5,2.5)   #(intercept,slope)
Sigma <- matrix(c(0.9, 0.4,
                  0.4, 0.9), 2, 2)


#parameters
n_vals <- c(100, 500, 1000)
clusters <- c(10, 50, 100)

grid <- expand.grid(
  num_clusters = clusters,
  n           = n_vals
)

summary <- data.frame(
  n            = grid$n,
  num_clusters = grid$num_clusters,
  Type1_b0_t = NA, Type1_b1_t = NA, Type1_b0_m = NA, Type1_b1_m = NA,
  Type1_b0_u = NA, Type1_b1_u = NA, Type1_b0_s = NA, Type1_b1_s = NA,
  Power_b0_t = NA, Power_b1_t = NA, Power_b0_m = NA, Power_b1_m = NA,
  Power_b0_u = NA, Power_b1_u = NA, Power_b0_s = NA, Power_b1_s = NA
)


#simulation
for (s in seq_along(n_vals)){
  n <- n_vals[s]
  
  cat("Simulation for n=", n, "\n")
  
  Rn_M <- ceiling(sqrt(n))
  Rn_U <- ceiling((n / 2)^(4 / 3))
  
  for (m in seq_along(clusters)){
    num_clusters <- clusters[m]
    cluster_size <- n/num_clusters
    cluster_id    <- rep(1:num_clusters, times = cluster_size)
    
    cat("Simulation for n_c=", num_clusters, "\n")
    
    results <- data.frame(Type1_b0_t = logical(iterations), Type1_b1_t = logical(iterations),
                          Power_b0_t = logical(iterations), Power_b1_t = logical(iterations),
                          Type1_b0_m = logical(iterations), Type1_b1_m = logical(iterations),
                          Power_b0_m = logical(iterations), Power_b1_m = logical(iterations),
                          Type1_b0_u = logical(iterations), Type1_b1_u = logical(iterations),
                          Power_b0_u = logical(iterations), Power_b1_u = logical(iterations),
                          Type1_b0_s = logical(iterations), Type1_b1_s = logical(iterations),
                          Power_b0_s = logical(iterations), Power_b1_s = logical(iterations))
    
    for (i in 1:iterations){
      if (i %% 500 == 0) cat("  Iteration", i, "of", iterations, "\n")
      
      #generate clustered data
      re_null <- MASS::mvrnorm(num_clusters, mu = c(0, 0), Sigma = Sigma)
      x_null <- rnorm(n)
      y_null  <- (beta_null[1] + re_null[cluster_id, 1]) +
        (beta_null[2] + re_null[cluster_id, 2]) * x_null +
        rnorm(n, 0, 1)
      fit_null <- lm(y_null~x_null)
      s_null <- summary(fit_null)
      
      re_alt <- MASS::mvrnorm(num_clusters, mu = c(0, 0), Sigma = Sigma)
      x_alt <- rnorm(n)
      y_alt <- (beta_alt[1] + re_alt[cluster_id, 1]) +
        (beta_alt[2] + re_alt[cluster_id, 2]) * x_alt +
        rnorm(n, 0, 1)
      fit_alt <- lm(y_alt~x_alt)
      s_alt <- summary(fit_alt)
      
      #size
      
      #t-test
      beta_hat_null <- c(coef(fit_null)[1], coef(fit_null)[2])
      se_null <- c(s_null$coefficients[1,"Std. Error"], s_null$coefficients[2,"Std. Error"])
      df_null = fit_null$df.residual
      t_null <- c((beta_hat_null[1] - beta_null[1]) / se_null[1],
                  (beta_hat_null[2] - beta_null[2]) / se_null[2])
      
      results$Type1_b0_t[i] <- abs(t_null[1]) > qt(0.975, df = df_null)
      results$Type1_b1_t[i] <- abs(t_null[2]) > qt(0.975, df = df_null)
      
      #mean- and u-type
      X_null_mat <- cbind(1,x_null)
      W_mat_null <- n * solve(t(X_null_mat) %*% X_null_mat) %*% t(X_null_mat)
      data_null <- cbind(W_mat_null[1,]*y_null - beta_null[1], 
                         W_mat_null[2,]*y_null - beta_null[2])
      
      tm_null <- c(Mean_stat(data_null[,1], mu=0, Rn=Rn_M), Mean_stat(data_null[,2], mu=0, Rn=Rn_M))
      tu_null <- c(U_stat(data_null[,1], mu=0, Rn=Rn_U),U_stat(data_null[,2], mu=0, Rn=Rn_U))
      
      results$Type1_b0_m[i] <- tm_null[1] > qchisq(1 - alpha, df = 1)
      results$Type1_b1_m[i] <- tm_null[2] > qchisq(1 - alpha, df = 1)
      results$Type1_b0_u[i] <- tu_null[1] > qnorm(1 - alpha)
      results$Type1_b1_u[i] <- tu_null[2] > qnorm(1 - alpha)
      
      #sandwich
      y_null_shift <- y_null - beta_null[1] - beta_null[2] * x_null
      fit_null_shift <- lm(y_null_shift ~ x_null)
      
      ps_null <- coef_test(fit_null_shift,
                           vcov    = vcovCR(fit_null_shift, cluster = cluster_id, type = "CR2"),
                           test    = "Satterthwaite")$p_Satt
      
      results$Type1_b0_s[i] <- ps_null[1] < alpha
      results$Type1_b1_s[i] <- ps_null[2] < alpha
      
      #power
      
      #t-test
      beta_hat_alt <- c(coef(fit_alt)[1], coef(fit_alt)[2])
      se_alt <- c(s_alt$coefficients[1,"Std. Error"], s_alt$coefficients[2,"Std. Error"])
      df_alt = fit_alt$df.residual
      t_alt <- c((beta_hat_alt[1] - beta_null[1]) / se_alt[1],
                 (beta_hat_alt[2] - beta_null[2]) / se_alt[2])
      
      results$Power_b0_t[i] <- abs(t_alt[1]) > qt(0.975, df = df_null)
      results$Power_b1_t[i] <- abs(t_alt[2]) > qt(0.975, df = df_null)
      
      #mean- and u-type
      X_alt_mat <- cbind(1,x_alt)
      W_mat_alt <- n * solve(t(X_alt_mat) %*% X_alt_mat) %*% t(X_alt_mat)
      data_alt <- cbind(W_mat_alt[1,]*y_alt - beta_null[1], 
                        W_mat_alt[2,]*y_alt - beta_null[2])
      
      tm_alt <- c(Mean_stat(data_alt[,1], mu=0, Rn=Rn_M), Mean_stat(data_alt[,2], mu=0, Rn=Rn_M))
      tu_alt <- c(U_stat(data_alt[,1], mu=0, Rn=Rn_U),U_stat(data_alt[,2], mu=0, Rn=Rn_U))
      
      results$Power_b0_m[i] <- tm_alt[1] > qchisq(1 - alpha, df = 1)
      results$Power_b1_m[i] <- tm_alt[2] > qchisq(1 - alpha, df = 1)
      results$Power_b0_u[i] <- tu_alt[1] > qnorm(1 - alpha)
      results$Power_b1_u[i] <- tu_alt[2] > qnorm(1 - alpha)
      
      #sandwich
      y_alt_shift <- y_alt - beta_null[1] - beta_null[2] * x_alt
      fit_alt_shift <- lm(y_alt_shift ~ x_alt)
      
      ps_alt <- coef_test(fit_alt_shift,
                          vcov    = vcovCR(fit_alt_shift, cluster = cluster_id, type = "CR2"),
                          test    = "Satterthwaite")$p_Satt
      
      results$Power_b0_s[i] <- ps_alt[1] < alpha
      results$Power_b1_s[i] <- ps_alt[2] < alpha
    }
    row_idx <- (s - 1) * length(clusters) + m
    
    summary$Type1_b0_t[row_idx] <- mean(results$Type1_b0_t)
    summary$Type1_b1_t[row_idx] <- mean(results$Type1_b1_t)
    summary$Type1_b0_m[row_idx] <- mean(results$Type1_b0_m)
    summary$Type1_b1_m[row_idx] <- mean(results$Type1_b1_m)
    summary$Type1_b0_u[row_idx] <- mean(results$Type1_b0_u)
    summary$Type1_b1_u[row_idx] <- mean(results$Type1_b1_u)
    summary$Type1_b0_s[row_idx] <- mean(results$Type1_b0_s)
    summary$Type1_b1_s[row_idx] <- mean(results$Type1_b1_s)
    
    summary$Power_b0_t[row_idx] <- mean(results$Power_b0_t)
    summary$Power_b1_t[row_idx] <- mean(results$Power_b1_t)
    summary$Power_b0_m[row_idx] <- mean(results$Power_b0_m)
    summary$Power_b1_m[row_idx] <- mean(results$Power_b1_m)
    summary$Power_b0_u[row_idx] <- mean(results$Power_b0_u)
    summary$Power_b1_u[row_idx] <- mean(results$Power_b1_u)
    summary$Power_b0_s[row_idx] <- mean(results$Power_b0_s)
    summary$Power_b1_s[row_idx] <- mean(results$Power_b1_s)
  }
}
#plot
library(ggplot2)
library(tidyr)
library(dplyr)


long <- summary |>
  pivot_longer(
    cols = -c(n, num_clusters),
    names_to  = "metric",
    values_to = "rate"
  ) |>
  mutate(
    type      = if_else(startsWith(metric, "Type1"), "Type I Error", "Power"),
    param     = if_else(grepl("b0", metric), "β_0 (intercept)", "β_1 (slope)"),
    test      = case_when(
      endsWith(metric, "_t") ~ "t-test",
      endsWith(metric, "_m") ~ "Mean",
      endsWith(metric, "_u") ~ "U-stat",
      endsWith(metric, "_s") ~ "Sandwich"
    ),
    n_label      = factor(paste0("n = ", n), levels = c("n = 100", "n = 500", "n = 1000")),
    num_clusters = factor(num_clusters)
  )


test_colors <- c(
  "t-test"   = "#E69F00",
  "Mean"     = "#2A5A8C",
  "U-stat"   = "#2E7D5C",
  "Sandwich" = "#B03E3E"
)

test_shapes <- c(
  "Mean"     = 19,
  "U-stat"   = 17,
  "t-test"   = 15,
  "Sandwich" = 18
)


base_theme <- list(
  scale_colour_manual(values = test_colors, name = "Test"),
  scale_shape_manual(values = test_shapes, name = "Test"),
  scale_linetype_manual(
    values = c("t-test" = "dotted", "Mean" = "dotted",
               "U-stat" = "dotted", "Sandwich" = "dotted"),
    name = "Test"
  ),
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)),
  labs(x = "Number of clusters", y = "Rejection rate"),
  theme_bw(base_size = 11),
  theme(
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(size = 9, face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5)
  )
)


p_size <- long |>
  filter(type == "Type I Error") |>
  ggplot(aes(x = num_clusters, y = rate,
             colour = test, group = test, linetype = test, shape = test)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.5) +
  geom_hline(aes(yintercept = 0.05),
             linetype = "dashed", colour = "#A01010", linewidth = 0.4) +
  facet_grid(rows = vars(param), cols = vars(n_label)) +
  labs(title = "Type I Error (Size)") +
  base_theme


p_power <- long |>
  filter(type == "Power") |>
  ggplot(aes(x = num_clusters, y = rate,
             colour = test, group = test, linetype = test, shape = test)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.5) +
  facet_grid(rows = vars(param), cols = vars(n_label)) +
  labs(title = "Power") +
  base_theme

ggsave("size_plot.png",  p_size,  width = 8, height = 6, dpi = 150)
ggsave("power_plot.png", p_power, width = 8, height = 6, dpi = 150)

p_size
p_power