#R_n control 

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
num_clusters <- 200
cluster_size <- 5
n <- num_clusters * cluster_size
cluster_id <- rep(1:num_clusters, each = cluster_size)

alpha <- 0.05
rho <- 1
iterations <- 2000

mu_0 <- 0
mu_alt <- 0.5

Rn_M_base <- ceiling(sqrt(n))
Rn_U_base <- ceiling((n / 2)^(4 / 3))

epsilons <- c(0.1,0.3,0.5,0.7,1,2,3,5,7,10)


results <- data.frame(epsilon = epsilons,
                      Power_Mean = NA_real_,
                      Power_U    = NA_real_,
                      Size_Mean = NA_real_,
                      Size_U    = NA_real_)

for (k in seq_along(epsilons)){
  epsilon <- epsilons[k]
  Rn_M <- pmax(2, floor(epsilon*Rn_M_base))
  Rn_U <- pmax(2, floor(epsilon*Rn_U_base))
  
  rej_Mean_power <- logical(iterations)
  rej_U_power    <- logical(iterations)
  rej_Mean_size  <- logical(iterations)
  rej_U_size     <- logical(iterations)
  
  for (i in 1:iterations){
    
    #generate clusterd data
    alpha_f_null <- rnorm(num_clusters,mu_0,sqrt(rho))
    eps_null <- rnorm(n)
    X_null <- alpha_f_null[cluster_id] + eps_null
    
    alpha_f_alt <- rnorm(num_clusters,mu_alt,sqrt(rho))
    eps_alt <- rnorm(n)
    X_alt <- alpha_f_alt[cluster_id] + eps_alt
    
    if (iid){
      X_null <- rnorm(n,mu_0,sqrt(rho))
      X_alt <- rnorm(n,mu_alt,sqrt(rho))
    }
    
    #Type-1 error
    tm_null <- Mean_stat(X_null, mu=mu_0, Rn=Rn_M)
    tu_null <- U_stat(X_null, mu=mu_0, Rn=Rn_U)
    
    rej_Mean_size[i] <- tm_null > qchisq(1 - alpha, df = 1)
    rej_U_size[i] <- tu_null > qnorm(1 - alpha)
    
    #power
    tm_alt <- Mean_stat(X_alt, mu=mu_0, Rn=Rn_M)
    tu_alt <- U_stat(X_alt, mu=mu_0, Rn=Rn_U)
    
    rej_Mean_power[i] <- tm_alt > qchisq(1 - alpha, df = 1)
    rej_U_power[i] <- tu_alt > qnorm(1 - alpha)
  }
  results$Power_Mean[k] <- mean(rej_Mean_power)
  results$Power_U[k] <- mean(rej_U_power)
  results$Size_Mean[k] <- mean(rej_Mean_size)
  results$Size_U[k] <- mean(rej_U_size)
}

#plot
par(mar = c(5, 4, 4, 2) + 0.1)

plot(NULL,
     xlim = range(epsilons),
     ylim = c(0, 1),
     xlab = expression(paste("R"[n], " multiplier  ", epsilon)),
     ylab = "Rejection rate",
     log  = "x")   

grid(nx = NULL, ny = NULL, col = "grey85", lty = 1, lwd = 0.7)

# reference lines
abline(h = alpha, lty = 2, col = "#A01010")
abline(v = 1,     lty = 3, col = "#A01010")
text(1.05, 0.97, expression(epsilon == 1), col = "grey40", cex = 0.8, adj = 0)

# size curves (solid)
lines(epsilons, results$Size_Mean, col = "#2A5A8C", lwd = 2, lty = 1)
lines(epsilons, results$Size_U,   col = "#B03E3E",    lwd = 2, lty = 1)

# power curves (dashed)
lines(epsilons, results$Power_Mean, col = "#2A5A8C", lwd = 2, lty = 3)
lines(epsilons, results$Power_U,    col = "#B03E3E",    lwd = 2, lty = 3)

# points to make individual epsilon values visible
points(epsilons, results$Size_Mean,  col = "#2A5A8C", pch = 16)
points(epsilons, results$Size_U,     col = "#B03E3E",    pch = 16)
points(epsilons, results$Power_Mean, col = "#2A5A8C", pch = 1)
points(epsilons, results$Power_U,    col = "#B03E3E",    pch = 1)

legend("topleft",
       legend = c(expression(T[M] ~ size),
                  expression(T[U] ~ size),
                  expression(T[M] ~ power),
                  expression(T[U] ~ power),
                  "nominal size",
                  expression(epsilon == 1)),
       col    = c("#2A5A8C", "#B03E3E", "#2A5A8C", "#B03E3E", "#A01010", "#A01010"),
       lty    = c(1, 1, 2, 2, 2, 3),
       lwd    = c(2, 2, 2, 2, 1, 1),
       pch    = c(16, 16, 1, 1, NA, NA),
       bty    = "n",
       cex    = 0.85)
