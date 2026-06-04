#leung vs club sandwich

library(clubSandwich)

#sanity check
iid = FALSE

#Mean-type statistic
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

#U-type statistic
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

#Simulation

#parameters
scenarios <- c(20,50,100) #num of clusters

alpha <- 0.05
iterations <- 2000
rho <- 1

mu_0 <- 0
mu_alt <- 0.5


summary <- data.frame(
  num_clusters = scenarios,
  Type1_Mean = NA, Type1_U = NA, Type1_t = NA, Type1_S = NA,
  Power_Mean = NA, Power_U = NA, Power_t = NA, Power_S = NA
)

for (s in seq_along(scenarios)) {
  num_clusters <- scenarios[s]
  cat("\n Scenario: num_clusters =", num_clusters, " \n")
  
  cluster_sizes <- sample(4:8, num_clusters, replace = TRUE)
  cluster_id    <- rep(1:num_clusters, times = cluster_sizes)
  n             <- length(cluster_id)
  
  Rn_M <- ceiling(sqrt(n))
  Rn_U <- ceiling((n / 2)^(4 / 3))
  
  res <- data.frame(
    Type1_Mean = logical(iterations), Type1_U = logical(iterations),
    Type1_t    = logical(iterations), Type1_S  = logical(iterations),
    Power_Mean = logical(iterations), Power_U  = logical(iterations),
    Power_t    = logical(iterations), Power_S  = logical(iterations)
  )
  
  for (i in 1:iterations) {
    if (i %% 500 == 0) cat("  Iteration", i, "of", iterations, "\n")
    
    # Generate data
    if (iid) {
      X_null <- rnorm(n, mu_0, sqrt(rho))
      X_alt  <- rnorm(n, mu_alt, sqrt(rho))
    } else {
      X_null <- rnorm(num_clusters, mu_0,   sqrt(rho))[cluster_id] + rnorm(n)
      X_alt  <- rnorm(num_clusters, mu_alt, sqrt(rho))[cluster_id] + rnorm(n)
    }
    
    #type 1 error
    tm_null <- Mean_stat(X_null, mu = mu_0, Rn = Rn_M)
    tu_null <- U_stat(X_null,    mu = mu_0, Rn = Rn_U)
    t_null  <- t.test(X_null, mu = mu_0)$p.value
    
    fit_null <- lm(X_null ~ 1)
    s_null   <- coef_test(fit_null,
                          vcov    = vcovCR(fit_null, cluster = cluster_id, type = "CR2"),
                          test    = "Satterthwaite")$p_Satt
    
    res$Type1_Mean[i] <- tm_null > qchisq(1 - alpha, df = 1)
    res$Type1_U[i]    <- tu_null > qnorm(1 - alpha)
    res$Type1_t[i]    <- t_null  < alpha
    res$Type1_S[i]    <- s_null  < alpha
    
    #power
    tm_alt <- Mean_stat(X_alt, mu = mu_0, Rn = Rn_M)
    tu_alt <- U_stat(X_alt,    mu = mu_0, Rn = Rn_U)
    t_alt  <- t.test(X_alt, mu = mu_0)$p.value
    
    fit_alt <- lm(X_alt ~ 1)
    s_alt   <- coef_test(fit_alt,
                         vcov    = vcovCR(fit_alt, cluster = cluster_id, type = "CR2"),
                         test    = "Satterthwaite")$p_Satt
    
    res$Power_Mean[i] <- tm_alt > qchisq(1 - alpha, df = 1)
    res$Power_U[i]    <- tu_alt > qnorm(1 - alpha)
    res$Power_t[i]    <- t_alt  < alpha
    res$Power_S[i]    <- s_alt  < alpha
  }
  
  summary[s, "Type1_Mean"] <- mean(res$Type1_Mean)
  summary[s, "Type1_U"]    <- mean(res$Type1_U)
  summary[s, "Type1_t"]    <- mean(res$Type1_t)
  summary[s, "Type1_S"]    <- mean(res$Type1_S)
  summary[s, "Power_Mean"] <- mean(res$Power_Mean)
  summary[s, "Power_U"]    <- mean(res$Power_U)
  summary[s, "Power_t"]    <- mean(res$Power_t)
  summary[s, "Power_S"]    <- mean(res$Power_S)
}
    
  
#plot
scenario_cols  <- c("#2A5A8C", "#2E7D5C", "#B03E3E")
scenario_pch   <- c(19, 17, 15, 18)          # one shape per statistic
stat_labels    <- c("Mean-type", "U-type", "t-test", "clubSandwich")

all_type1 <- as.matrix(summary[, c("Type1_Mean","Type1_U","Type1_t","Type1_S")])
all_power <- as.matrix(summary[, c("Power_Mean","Power_U","Power_t","Power_S")])

par(mar = c(5.5, 4.5, 2, 1), mgp = c(3, 0.7, 0), tcl = -0.3)

plot(NA,
     xlab = "Type I error (rejection rate under H0)",
     ylab = "Power (rejection rate under H1)",
     xlim = c(0, max(all_type1) * 1.35),
     ylim = c(0, 1),
     las  = 1,
     bty  = "o")

grid(nx = NULL, ny = NULL, col = "grey85", lty = 1, lwd = 0.7)

abline(v = 0.05, col = "#A01010", lty = 2, lwd = 1.2)

for (s in seq_along(scenarios)) {
  col    <- scenario_cols[s]
  type1s <- all_type1[s, ]
  powers <- all_power[s, ]
  
  ord <- order(type1s)
  lines(type1s[ord], powers[ord], col = col, lty = 2, lwd = 0.9)
  
  points(type1s, powers,
         pch = scenario_pch,
         col = col, bg = col,
         cex = 1.5, lwd = 1.8)
}

# Legend: scenarios (color)
legend("topright",
       legend = paste0("n_c = ", scenarios),
       col    = scenario_cols,
       lty    = 2, lwd    = 0.9,
       pch    = 19, pt.cex = 1.3,
       title  = expression(bold("No. clusters")),
       bty    = "n", cex = 0.82)

# Legend: statistics (shape)
legend("bottomright",
       legend = stat_labels,
       pch    = scenario_pch,
       col    = "grey30",
       pt.cex = 1.3, pt.lwd = 1.8,
       title  = expression(bold("Statistic")),
       bty    = "n", cex = 0.82)



