

# PRELIMINARY FUNCTIONS -------------------------------------------------------

# Function to read in all required packages in one go:
loadPackages <- function(x) {
  for(i in x) {
    if(!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

# Load the packages
loadPackages(c("tidyverse", "parallel", "foreach", "doParallel", 
               "Rfast", "data.table"))

# FUNCTIONS TO CREATE SAMPLE MATRICES -----------------------------------------

scrambled_sobol <- function(matrices, A, B, order, cluster) {
  first <- 1:ncol(A)
  N <- nrow(A)
  if(order == "first") {
    loop <- first
  } else if(order == "second" | order == "third") {
    second <- c(first, utils::combn(1:ncol(A), 2, simplify = FALSE))
    loop <- second
  } else if(order == "third") {
    third <- c(second, utils::combn(1:ncol(A), 3, simplify = FALSE))
    loop <- third
  } else {
    stop("order should be either first, second or third")
  }
  if(is.null(cluster) == FALSE) {
    loop <- cluster
  }
  AB.mat <- "AB" %in% matrices
  BA.mat <- "BA" %in% matrices
  if(AB.mat == TRUE) {
    X <- rbind(A, B)
    for(i in loop) {
      AB <- A
      AB[, i] <- B[, i]
      X <- rbind(X, AB)
    }
    AB <- X[(2 * N + 1):nrow(X), ]
  } else if(AB.mat == FALSE) {
    AB <- NULL
  }
  if(BA.mat == TRUE) {
    W <- rbind(A, B)
    for(i in loop) {
      BA <- B
      BA[, i] <- A[, i]
      W <- rbind(W, BA)
    }
    BA <- W[(2 * N + 1) : nrow(W), ]
  } else if(BA.mat == FALSE) {
    BA <- NULL
  }
  return(rbind(AB, BA))
}

sobol_matrices <- function(matrices = c("A", "B", "AB"), 
                           N, params, order = "first",
                           cluster = NULL) {
  if(length(matrices) >= 4 & !order == "first" ) {
    stop("higher orders should be computed with an A, B and AB or BA matrices")
  }
  k <- length(params)
  df <- randtoolbox::sobol(n = N, dim = k * 2)
  A <- df[, 1:k]
  B <- df[, (k + 1) : (k * 2)]
  out <- scrambled_sobol(matrices = matrices, 
                         A = A, B = B, order = order, 
                         cluster = cluster)
  A.mat <- "A" %in% matrices
  B.mat <- "B" %in% matrices
  if(A.mat == FALSE) {
    A <- NULL
  }
  if(B.mat == FALSE) {
    B <- NULL
  }
  final <- rbind(A, B, out)
  colnames(final) <- params
  return(final) 
}


# COMPUTATION OF SOBOL' INDICES -----------------------------------------------

sobol_boot <- function(d, i, N, params, R, first, total, order) {
  k <- length(params)
  m <- d[i, ]
  if(first == "monod" & !total == "monod" | !first == "monod" & total == "monod") {
    stop("monod should be used simultaneously in first and total indices 
         \n with an A, AB and BA matrices")
  }
  if(first == "azzini" & !total == "azzini" | !first == "azzini" & total == "azzini") {
    stop("azzini should be used simultaneously in first and total indices 
         \n with an A, B, AB and BA matrices")
  }
  if(first == "jansen" | first == "saltelli" & 
     total == "jansen" | total == "sobol" | total == "homma") {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, -c(1, 2)]
    f0 <- 1 / N * sum(m[, 1])
    VY <- 1 / N * sum((m[, 1] - f0) ^ 2)
  } else if(first == "monod" & total == "monod") {
    Y_A <- m[, 1]
    Y_AB <- m[, 2:(2 + k - 1)]
    Y_BA <- m[, (ncol(m) - k + 1):ncol(m)]
  } else if(first == "azzini" | total == "azzini") {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, 3:(3 + k - 1)]
    Y_BA <- m[, (ncol(m) - k + 1):ncol(m)]
  }
  if(order == "second" | order == "third") {
    if(first == "azzini" | first == "monod") {
      stop("high-order indices require to use either jansen or saltelli as 
         \n first-order indices")
    }
  }
  if(first == "saltelli") {
    Vi <- 1 / N * Rfast::colsums(Y_B * (Y_AB - Y_A))
  } else if(first == "jansen") {
    Vi <- VY - 1 / (2 * N) * Rfast::colsums((Y_B - Y_AB) ^ 2)
  } else if(first == "azzini") {
    Vi <- Rfast::colsums((Y_A - Y_AB) * (Y_B - Y_BA))
  } else if(first == "monod") {
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA - 1 / (2 * N) * Rfast::colsums(Y_A + Y_BA))
  } else {
    stop("Si should  be either jansen, saltelli, azzini or monod")
  } 
  if(first == "saltelli" | first == "jansen") {
    Si <- Vi[1:k] / VY
  } 
  if(first == "azzini") {
    Si <- abs(Vi / (1 / 2 * Rfast::colsums((Y_A - Y_B) ^ 2 + (Y_AB - Y_BA) ^ 2))) 
  } 
  if (first == "monod") { # Monod
    Si <- Vi / (1 / (2 * N) * Rfast::colsums(Y_A ^ 2 + Y_BA ^ 2) - 
                  1 / (2 * N) * Rfast::colsums(Y_A + Y_BA))
  }
  if(total == "jansen") {
    STi <- (1 / (2 * N) * Rfast::colsums((Y_A - Y_AB) ^ 2)) / VY
  } else if(total == "homma") {
    STi <- (VY - (1 / N) * Rfast::colsums(Y_A * Y_AB) + f0 ^ 2) / VY
  } else if(total == "sobol") {
    STi <- ((1 / N) * Rfast::colsums(Y_A * (Y_A - Y_AB))) / VY
  } else if(total == "azzini") {
    STi <- 1 - abs(Rfast::colsums((Y_A - Y_BA) * (Y_B - Y_AB)) / 
                     (1 / 2 * Rfast::colsums((Y_A - Y_B) ^ 2 + (Y_AB - Y_BA) ^ 2)))
  } else if(total == "monod") {
    STi <- 1 - (1 / N * Rfast::colsums(Y_A * Y_AB - 1 / (2 * N) * Rfast::colsums(Y_A + Y_AB)) / 
                  ((1 / (2 * N) * Rfast::colsums(Y_A ^ 2 + Y_AB ^ 2) - 
                      1 / (2 * N) * Rfast::colsums(Y_A + Y_AB))))
  } else {
    stop("STi should  be either jansen, sobol, azzini, homma or monod")
  }
  STi <- STi[1:k]
  if(order == "second" | order == "third") {
    com2 <- utils::combn(1:k, 2, simplify = FALSE)
    mat2 <- t(mapply(c, Vi[(k + 1):(k + length(com2))], lapply(com2, function(x) Vi[x])))
    Vij <- apply(mat2, 1, function(x) Reduce("-", x))
    Sij <- Vij / VY
  } else {
    Sij <- NULL
  }
  if(order == "third") {
    tmp <- do.call(rbind, com2)
    Vij.vec <- as.numeric(paste(tmp[, 1], tmp[, 2], sep = ""))
    Vij.named <- Vij
    names(Vij.named) <- Vij.vec
    com3 <- utils::combn(1:k, 3, simplify = FALSE)
    Vi.only <- do.call(rbind, lapply(com3, function(x) Vi[x])) # Extract Vi, Vj, Vk
    Vijk.only <- tail(Vi, length(com3)) # Extract Vijk
    tmp3 <- do.call(rbind, com3)
    first.pairwise <- lapply(paste(tmp3[, 1], tmp3[, 2], sep = ""), function(x) Vij.named[x])
    second.pairwise <- lapply(paste(tmp3[, 1], tmp3[, 3], sep = ""), function(x) Vij.named[x]) 
    third.pairwise <- lapply(paste(tmp3[, 2], tmp3[, 3], sep = ""), function(x) Vij.named[x])
    Vij.only <- t(mapply(cbind, first.pairwise, second.pairwise, third.pairwise))
    mat3 <- cbind(Vijk.only, Vij.only, Vi.only)
    Vijk <- apply(mat3, 1, function(x) Reduce("-", x))
    Sijk <- Vijk / VY
  } else {
    Sijk <- NULL
  }
  return(c(Si, STi, Sij, Sijk))
}

# BOOTSTRAP OF SOBOL' INDICES -------------------------------------------------

# Extract bootstrap confidence intervals
bootstats <- function(b, conf = conf, type = type) {
  p <- length(b$t0)
  lab <- c("original", "bias", "std.error", "low.ci", "high.ci")
  tmp <- as.data.frame(matrix(nrow = p,
                              ncol = length(lab),
                              dimnames = list(NULL, lab)))
  for (i in 1:p) {
    # original estimation, bias, standard deviation
    tmp[i, "original"] <- b$t0[i]
    tmp[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    tmp[i, "std.error"] <- stats::sd(b$t[, i])
    # confidence interval
    if (type == "norm") {
      ci <- boot::boot.ci(b, index = i, type = "norm", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$norm[2]
        tmp[i, "high.ci"] <- ci$norm[3]
      }
    } else if (type == "basic") {
      ci <- boot::boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$basic[4]
        tmp[i, "high.ci"] <- ci$basic[5]
      }
    } else if (type == "percent") {
      ci <- boot::boot.ci(b, index = i, type = "perc", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$percent[4]
        tmp[i, "high.ci"] <- ci$percent[5]
      }
    } else if (type == "bca") {
      ci <- boot::boot.ci(b, index = i, conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$bca[4]
        tmp[i, "high.ci"] <- ci$bca[5]
      }
    }
  }
  return(data.table::data.table(tmp))
}

# Compute bootstrap Sobol' indices
sobol_indices <- function(Y, N, params, R, first = "saltelli", total = "jansen", 
                          order = "first", parallel = "no", 
                          ncpus = 1, conf = 0.95, type = "norm") {
  k <- length(params)
  d <- matrix(Y, nrow = N)
  out.boot <- boot::boot(data = d, statistic = sobol_boot, R = R, N = N, params = params, 
                         first = first, total = total, order = order, 
                         parallel = parallel, ncpus = ncpus)
  out.ci <- bootstats(out.boot, conf = conf, type = type)
  if(order == "first") {
    parameters <- c(rep(params, times = 2))
    indices <- c(rep(c("Si", "Ti"), each = k))
  } else if(order == "second" | order == "third") {
    unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x) 
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    indices <- c(rep(c("Si", "Ti"), each = length(params)), 
                 rep("Sij", times = length(vector.second)))
  } else if(order == "third") {
    vector.third <- unlist(lapply(utils::combn(params, 3, simplify = FALSE), function(x) 
      paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.third)
    indices <- c(rep(c("Si", "Ti"), each = k), 
                 rep("Sij", times = length(vector.second)), 
                 rep("Sijk", times = length(vector.third)))
  }
  out.ci[, sensitivity:= cbind(indices)][, parameters:= cbind(parameters)]
  return(out.ci)
}

# Create function for custom plot themes
theme_AP <- function() {
  theme_bw() +
    theme(legend.position = "top", 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA))
}

# Plot indices
plot_sobol <- function(indices, dummy = NULL, order = "first") {
  sensitivity <- low.ci <- high.ci <- parameters <- original <- NULL
  if(order == "first") {
    if(is.null(dummy) == FALSE) {
      plot.dummy <- geom_rect(data = dummy,
                              aes(ymin = 0,
                                  ymax = high.ci,
                                  xmin = -Inf,
                                  xmax = Inf,
                                  fill = sensitivity),
                              alpha = 0.2,
                              inherit.aes = FALSE)
    } else {
      plot.dummy <- NULL
    }
    p <- indices[sensitivity == "Si" | sensitivity == "Ti"]
    gg <- ggplot2::ggplot(p, aes(parameters, original,
                                 fill = sensitivity)) +
      geom_bar(stat = "identity",
               position = position_dodge(0.6),
               color = "black") +
      plot.dummy +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci),
                    position = position_dodge(0.6)) +
      scale_fill_discrete(name = "Sobol' indices",
                          labels = c(expression(S[italic(i)]),
                                     expression(T[italic(i)]))) +
      labs(x = "",
           y = "Sobol' index") +
      theme_AP()
  } else if(!order == "first") {
    if(order == "second") {
      plot.type <- "Sij"
    } else if(order == "third") {
      plot.type <- "Sijk"
    } else {
      stop("Order should be either first, second or third")
    }
    p <- indices[sensitivity == plot.type]
    gg <- ggplot2::ggplot(p, aes(stats::reorder(parameters, original),
                                 original)) +
      geom_point() +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci)) +
      theme_bw() +
      labs(x = "",
           y = "Sobol' index") +
      geom_hline(yintercept = 0,
                 lty = 2,
                 color = "red") +
      theme_AP()
  }
  return(gg)
}

# CREATE METAFUNCTION ---------------------------------------------

function_list <- list(
  Linear = function(x) x,
  Quadratic = function(x) x ^ 2,
  Cubic = function(x) x ^ 3,
  Exponential = function(x) exp(1) ^ x / (exp(1) - 1),
  Periodic = function(x) sin(2 * pi * x) / 2,
  Discontinuous = function(x) ifelse(x > 0.5, 1, 0),
  No.effect = function(x) x * 0,
  Non.monotonic = function(x) 4 * (x - 0.5) ^ 2,
  Inverse = function(x) (10 - 1 / 1.1) ^ -1 * (x + 0.1) ^ - 1
)

metafunction <- function(mat) {
  all_functions <- names(function_list)
  tests <- sample(all_functions, 6)
  # 3 Coefficients from a Gaussian mixture
  components <- sample(1:3, prob = c(0.3, 0.5, 0.2), size = N, replace = TRUE)
  mus <- c(0, 10, 3)
  sds <- sqrt(c(1, 1, 0.1))
  coefficients <- sample(rnorm(N) * sds[components] + mus[components], 3)
  # Compute metafunction
  output <- vector()
  for(i in 1:nrow(mat)) {
    # First order
    output[[i]] <- sum(coefficients[[1]] * function_list[[tests[[1]]]](mat[i, ])) + 
      # Second order
      sum(coefficients[[2]] * function_list[[tests[[2]]]](mat[i, seq_len(2)]) * 
            function_list[[tests[[3]]]](mat[i, seq_len(2)])) + 
      # Third order
      sum(coefficients[[3]] * function_list[[tests[[4]]]](mat[i, seq_len(3)]) * 
            function_list[[tests[[5]]]](mat[i, seq_len(3)]) * 
            function_list[[tests[[6]]]](mat[i, seq_len(3)]))
  }
  final <- list(coefficients, tests, output)
  names(final) <- c("coefficients", "functions", "output")
  return(final)
}

# PLOT METAFUNCTION -----------------------------------------------------------

ggplot(data.frame(x = runif(100)), aes(x)) +
  map(1:length(function_list), function(nn) {
    stat_function(fun = function_list[[nn]], 
                  geom = "line", 
                  aes_(color = factor(names(function_list[nn])), 
                       linetype = factor(names(function_list[nn]))))
  }) + 
  labs(color= "Function", linetype = "Function") +
  theme_AP()

# CREATE SAMPLE MATRIX --------------------------------------------------------

N <- 500
R <- 100 # Number of bootstrap replicas for Sobol' indices
params <- c("k", "N") # Dimension and initial sample size
matrices <- c("A", "B", "AB", "BA")
mat <- randtoolbox::sobol(N, length(params))
mat[, 1] <- floor(qunif(mat[, 1], 3, 100))
mat[, 2] <- floor(qunif(mat[, 2], 100, 1000))
colnames(mat) <- params

# CREATE MODEL ----------------------------------------------------------------

model_battle <- function(k, N) {
  mt <- sobol_matrices(N = N, params = paste("X", 1:k, sep = ""), 
                       matrices = matrices)
  out <- metafunction(mt)
  return(out)
}

# RUN MODEL -------------------------------------------------------------------

# Set number of cores at 75%
n_cores <- floor(detectCores() * 0.75)

# Define parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Compute
Y <- foreach(i=1:nrow(mat), 
             .packages = "Rfast") %dopar%
  {
    model_battle(N = mat[[i, "N"]],
                 k = mat[[i, "k"]])
  }
# Stop parallel cluster
stopCluster(cl)

# # EXTRACT OUTPUT ------------------------------------------------------------

estimators_AB <- c("jansen", "sobol", "homma")

# SI:
## jansen, saltelli: A, B, AB

# STI:
## jansen, sobol, homma: A, B, AB
## monod/janon: A, AB, BA
## azzini/rosati: A, B, AB, BA.

# Extract output
A_B_AB <- A_AB_BA <- list()
for(i in 1:nrow(mat)) {
  # A, B and AB matrices
  A_B_AB[[i]] <- Y[[i]]$output[1:(mat[i, "N"] * (mat[i, "k"] + 2))]
  # A, AB and BA matrices
  A_AB_BA[[i]] <- Y[[i]]$output[c(1:mat[i, "N"], 
                                  (2 * mat[i, "N"] + 1):length(Y[[i]]$output))]
}

# COMPUTE SOBOL' INDICES ------------------------------------------------------

ind.jansen.sobol.homma <- ind.monod <- ind.azzini <- list()
for(i in 1:nrow(mat)) {
  ind.jansen.sobol.homma[[i]] <- lapply(estimators_AB, function(x) 
    sobol_indices(A_B_AB[[i]], 
                  N = mat[i, "N"], 
                  params = paste("X", 1:mat[i, "k"], sep = ""), 
                  total = x,
                  R = R, 
                  parallel = "multicore", 
                  ncpus = n_cores))
  ind.azzini[[i]] <- sobol_indices(Y[[i]]$output, 
                                   N = mat[i, "N"], 
                                   params = paste("X", 1:mat[i, "k"], sep = ""), 
                                   first = "azzini", 
                                   total = "azzini", 
                                   R = R, 
                                   parallel = "multicore", 
                                   ncpus = n_cores)
}





N <- 1000
order = "first"
A <- sobol_matrices(N = N, 
                    params = paste("X", 1:3, sep = ""), 
                    matrices = c("A", "B", "AB"), 
                    order = order)
out <- metafunction(A)
ind <- sobol_indices(Y = out$output, N = N, params = paste("X", 1:3, sep = ""), R = 100, 
                     order=order )
plot_sobol(ind, order = order)


out <- sensobol::ishigami_Mapply(A)
ind <- sobol_indices(Y = out, N = N, params = paste("X", 1:3, sep = ""), R = 100, 
                     first = "monod", total = "monod", order=order )
plot_sobol(ind, order = order)






































plot_sobol(ind[[39]])


# CHECK THE K = 2 ISSUE
###########################################


liu <- function(X1, X2) {
  X1 / X2
}
liu_Mapply <- function(X) {
  X[, 1] <- qchisq(X[, 1], df = 10)
  X[, 2] <- qchisq(X[, 2], df = 13.978)
  return(mapply(liu, X[, 1], X[, 2]))
}


N <- 1000
params <- paste("X", 1:3, sep = "")
R <- 100
order = "second"
A <- sobol_matrices(N = N, params = params, order = order, matrices = c("A", "B", "AB"))
Y <- sensobol::ishigami_Mapply(A)
ind <- sobol_indices(Y = Y, N = N, params = params, R = R, order = order)


# Compute bootstrap Sobol' indices
sobol_indices <- function(Y, N, params, R, first = "saltelli", total = "jansen", 
                          order = "first", parallel = "no", 
                          ncpus = 1, conf = 0.95, type = "norm") {
  k <- length(params)
  d <- matrix(Y, nrow = N)
  out.boot <- boot::boot(data = d, statistic = sobol_boot, R = R, N = N, params = params, 
                         first = first, total = total, order = order, 
                         parallel = parallel, ncpus = ncpus)
  return(out.boot)
}

N <- 10000
params <- paste("X", 1:3, sep = "")
R <- 100
order = "second"
A <- sobol_matrices(N = N, params = params, order = order, matrices = c("A", "B", "AB"))
Y <- sensobol::ishigami_Mapply(A)
ind <- sobol_indices(Y = Y, N = N, params = params, R = R, order = order)
plot_sobol(ind, order = "second")



sobol_indices <- function(Y, N, params, R, first = "saltelli", total = "jansen", 
                          order = "first", parallel = "no", 
                          ncpus = 1, conf = 0.95, type = "norm") {
  k <- length(params)
  d <- matrix(Y, nrow = N)
  out.boot <- boot::boot(data = d, statistic = sobol_boot, R = R, N = N, params = params, 
                         first = first, total = total, order = order, 
                         parallel = parallel, ncpus = ncpus)
  out.ci <- bootstats(out.boot, conf = conf, type = type)
  if(order == "first") {
    parameters <- c(rep(params, times = 2))
    indices <- c(rep(c("Si", "Ti"), each = k))
  } else if(order == "second" | order == "third") {
    vector.second <- paste0(unlist(utils::combn(params, 2, simplify = FALSE)), collapse = ".")
    parameters <- c(c(rep(params, times = 2)), vector.second)
    indices <- c(rep(c("Si", "Ti"), each = length(params)), 
                 rep("Sij", times = length(vector.high)))
  } else if(order == "third") {
    vector.third <- paste0(unlist(utils::combn(params, 3, simplify = FALSE)), collapse = ".")
    parameters <- c(parameters, vector.third)
    indices <- c(rep(c("Si", "Ti"), each = k), 
                 rep("Sij", times = length(vector.second)), 
                 rep("Sijk", times = length(vector.third)))
  }
  out.ci[, sensitivity:= cbind(indices)][, parameters:= cbind(parameters)]
  return(out.ci)
}



sensobol::sobol_matrices(n = 5, k = 2, second = TRUE)


sobol










paste0(unlist(utils::combn(params, 2, simplify = FALSE)), 
       collapse = ".")

vector.high <- paste0(unlist(utils::combn(params, 2, simplify = FALSE)), 
                      collapse = ".")










vector.high <- lapply(2:3, function(x) {
  perm <- utils::combn(params, x, simplify = FALSE)
  out <- unlist(lapply(perm, function(x) paste0(x, collapse = "."))) 
})
indices <- c(rep(c("Si", "Ti"), each = k), 
             rep("Sij", times = length(vector.high[[1]])), 
             rep("Sijk", times = length(vector.high[[2]])))
parameters <- c(rep(params, times = 2), vector.high[[1]], vector.high[[2]])
if(order == "first") {
  out.ci[, sensitivity:= cbind(indices[1:(k * 2)])][, parameters:= cbind(parameters[1:(k * 2)])]
} else if(order == "second") {
  out.ci[, sensitivity:= cbind(c(indices[1:(k * 2)], indices[indices == "Sij"]))][
    , parameters:= cbind(head(parameters, -length(vector.high[[2]])))
    ]
} else if(order == "third") {
  out.ci[, sensitivity:= cbind(c(indices[1:(k * 2)], indices[indices == "Sij"], indices[indices == "Sijk"]))][
    , parameters:= cbind(parameters)
    ]
}
return(out.ci)
}
















sobol_indices <- function(Y, N, params, R, first = "saltelli", total = "jansen", 
                          order = "first", parallel = "no", 
                          ncpus = 1, conf = 0.95, type = "norm") {
  k <- length(params)
  d <- matrix(Y, nrow = N)
  out.boot <- boot::boot(data = d, statistic = sobol_boot, R = R, N = N, params = params, 
                         first = first, total = total, order = order, 
                         parallel = parallel, ncpus = ncpus)
  out.ci <- bootstats(out.boot, conf = conf, type = type)
  if(length(k) == 2) {
    vector.high <- paste0(unlist(utils::combn(params, 2, simplify = FALSE)), 
                          collapse = ".")
    parameters <- c(rep(params, times = 2), vector.high)
    indices <- c(rep(c("Si", "Ti"), each = k), 
                 rep("Sij", times = length(vector.high)))
  } else if(length(k) > 2) {
    vector.high <- lapply(2:3, function(x) {
      perm <- utils::combn(params, x, simplify = FALSE)
      out <- unlist(lapply(perm, function(x) paste0(x, collapse = "."))) 
    })
    parameters <- c(rep(params, times = 2), vector.high[[1]], vector.high[[2]])
    indices <- c(rep(c("Si", "Ti"), each = k), 
                 rep("Sij", times = length(vector.high[[1]])), 
                 rep("Sijk", times = length(vector.high[[2]])))
  }
  if(order == "first") {
    out.ci[, sensitivity:= indices][, parameters:= cbind(parameters)]
  } else if(order == "second") {
    out.ci[, sensitivity:= cbind(c(indices[1:(k * 2)], indices[indices == "Sij"]))][
      , parameters:= cbind(head(parameters, -length(vector.high[[2]])))
      ]
  } else if(order == "third") {
    out.ci[, sensitivity:= cbind(c(indices[1:(k * 2)], indices[indices == "Sij"], indices[indices == "Sijk"]))][
      , parameters:= cbind(parameters)
      ]
  }
  return(out.ci)
}






indices <- c(rep(c("Si", "Ti"), each = k), 
             rep("Sij", times = length(vector.high[[1]])), 
             rep("Sijk", times = length(vector.high[[2]])))





paste0(unlist(utils::combn(params, 2, simplify = FALSE)), collapse = ".")

vector.high <- lapply(2:3, function(x) {
  perm <- utils::combn(params, x, simplify = FALSE)
  out <- unlist(lapply(perm, function(x) paste0(x, collapse = "."))) 
})
indices <- c(rep(c("Si", "Ti"), each = k), 
             rep("Sij", times = length(vector.high[[1]])), 
             rep("Sijk", times = length(vector.high[[2]])))
vector.high <- lapply(2:3, function(x) {
  perm <- utils::combn(params, x, simplify = FALSE)
  out <- unlist(lapply(perm, function(x) paste0(x, collapse = "."))) 
})
parameters <- c(rep(params, times = 2), vector.high[[1]], vector.high[[2]])
if(order == "first") {
  out.ci[, sensitivity:= cbind(indices[1:(k * 2)])][, parameters:= cbind(parameters[1:(k * 2)])]
} else if(order == "second") {
  out.ci[, sensitivity:= cbind(c(indices[1:(k * 2)], indices[indices == "Sij"]))][
    , parameters:= cbind(head(parameters, -length(vector.high[[2]])))
    ]
} else if(order == "third") {
  out.ci[, sensitivity:= cbind(c(indices[1:(k * 2)], indices[indices == "Sij"], indices[indices == "Sijk"]))][
    , parameters:= cbind(parameters)
    ]
}
return(out.ci)
}


























vector.high <- lapply(2:3, function(x) {
  perm <- utils::combn(c("X1", "X2"), x, simplify = FALSE)
  out <- unlist(lapply(perm, function(x) paste0(x, collapse = "."))) 
})

































da <- sobol_indices(Y = A_B_AB[[500]], 
                    params = paste("X", 1:45, sep = ""), 
                    N = 158, 
                    R = 2)



do.call(rbind, lapply(Y, "[", "coefficients"))
do.call(rbind, lapply(Y, "[", "functions"))


Y[["coefficients"]]


do.call(rbind, lapply(Y, function(x) x$coefficients))
data.table(do.call(rbind, lapply(Y, function(x) x$functions)))

































# CREATE SAMPLE MATRIX --------------------------------------------------------

N <- 10
params <- c("k", "N") # Dimension and initial sample size
estimators <- c("jansen", "sobol", "homma", "janon", "azzini")
mat <- randtoolbox::sobol(N, length(params))
dt <- replicate(length(estimators), mat, simplify = FALSE) %>%
  do.call(rbind, .) %>%
  data.table() 

dt[, estimators:= cbind(rep(estimators, each = N))]
setnames(dt, c("V1", "V2"), c("k", "N"))
dt[, k:= floor(qunif(k, 2, 100))][
  , N:= floor(qunif(N, 100, 1000))]















mat <- data.table(randtoolbox::sobol(N, length(params)))
mat <- setnames(mat, paste("V", 1:length(params), sep = ""), c("k", "N"))
mat <- mat[, k:= floor(qunif(k, 2, 100))][
  , N:= floor(qunif(N, 100, 1000))]


mat[, 1] <- floor(qunif(mat[, 1], 2, 100))
mat[, 2] <- floor(qunif(mat[, 2], 100, 1000))
colnames(mat) <- params






























# STI:
## jansen, sobol, homma: A, B, AB
## monod/janon: A, AB, BA
## azzini/rosati: A, B, AB, BA.


N <- 5000
k <- 3
params <- paste("X", 1:k, sep = "")
R <- 100
matrices <- c("A", "B", "AB", "BA")


A <- sobol_matrices(N = N, params = params, matrices = matrices)
Y <- sensobol::ishigami_Mapply(A)
ind <- sobol_indices(Y = Y, N = N, params = params, R = R, 
                     first = "azzini", total = "azzini")
plot_sobol(ind)





library(tidyverse)









sobol_matrices(N = 10, params = paste("X", 1:2, sep = ""), 
               matrices = matrices)








params <- c("X1", "X2")
N <- 249
order <- "first"
k <- length(params)
matrices <- c("A", "B", "AB")
df <- randtoolbox::sobol(n = N, dim = k * 2)
A <- df[, 1:k]
B <- df[, (k + 1) : (k * 2)]

first <- 1:ncol(A)
if(order == "first") {
  loop <- first
} else if(order == "second" | order == "third") {
  second <- c(first, utils::combn(1:ncol(A), 2, simplify = FALSE))
  loop <- second
} else if(order == "third") {
  third <- c(second, utils::combn(1:ncol(A), 3, simplify = FALSE))
  loop <- third
}
<- "AB" %in% matrices
BA.mat <- "BA" %in% matrices
X <- rbind(A, B)
for(i in loop) {
  AB <- A
  AB[, i] <- B[, i]
  X <- rbind(X, AB)
}
AB <- X[(2 * N + 1):nrow(X), ]



rm(list = ls())

params <- paste("X", 1:2, sep = "")

N <- 10
sobol_matrices(N = N, params = paste("X", 1:2, sep = ""))

df <- randtoolbox::sobol(n = 249, dim = 2 * 2)
A <- df[, 1:2]
B <- df[, (2 + 1) : (2 * 2)]
first <- 1:ncol(A)
loop <- first
X <- rbind(A, B)
for(i in loop) {
  AB <- A
  AB[, i] <- B[, i]
  X <- rbind(X, AB)
}
AB <- X[(2 * 249 + 1):nrow(X), ]










mat <- sobol_matrices(N = 249, params = paste("X", 1:2, sep = ""))



rm(list = ls())


