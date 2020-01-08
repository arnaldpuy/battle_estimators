

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
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA - (1 / (2 * N) * Rfast::colsums(Y_A + Y_BA)) ^ 2)
  } else {
    stop("Si should  be either jansen, saltelli, azzini or monod")
  } 
  if(first == "saltelli" | first == "jansen") {
    Si <- Vi[1:k] / VY
  } 
  if(first == "azzini") {
    Si <- abs(Vi / (1 / 2 * Rfast::colsums((Y_A - Y_B) ^ 2 + (Y_AB - Y_BA) ^ 2))) 
  } 
  if (first == "monod") { 
    Si <- Vi / (1 / (2 * N) * Rfast::colsums(Y_A ^ 2 + Y_BA ^ 2) - 
                  (1 / (2 * N) * Rfast::colsums(Y_A + Y_BA)) ^ 2)
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
    STi <- 1 - (1 / N * Rfast::colsums(Y_A * Y_AB) - 
                  (1/ N * Rfast::colsums((Y_A + Y_AB) / 2)) ^ 2) / 
      (1 / N * Rfast::colsums((Y_A ^ 2 + Y_AB ^ 2) / 2) - 
         (1/ N * Rfast::colsums((Y_A + Y_AB) / 2)) ^ 2)
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

# CHECK THAT ALL TI ESTIMATORS WORK -------------------------------------------

# Settings
estimators <- c("jansen", "sobol", "homma", "azzini", "monod")
test_functions <- c("Ishigami", "Sobol'G", "Morris")
N <- R <- 500

# Run model
ind <- Y <- mt <- list() 
for(i in estimators) {
  for(j in test_functions) {
    if(i == "jansen" | i == "sobol" | i == "homma") {
      matrices <- c("A", "B", "AB")
      first <- "saltelli"
    } else if(i == "azzini") {
      matrices <- c("A", "B", "AB", "BA")
      first = "azzini"
    } else if(i == "monod") {
      matrices <- c("A", "AB", "BA")
      first = "monod"
    }
    if(j == "Ishigami") {
      k <- 3
      modelRun <- sensobol::ishigami_Mapply
    } else if(j == "Sobol'G") {
      k <- 8
      modelRun <- sensobol::sobol_Fun
    } else if(j == "Morris") {
      k <- 20
      modelRun <- sensitivity::morris.fun
    }
    mt[[i]][[j]] <- sobol_matrices(N = N, params = paste("X", 1:k, sep = ""), matrices = matrices)
    Y[[i]][[j]] <- modelRun(mt[[i]][[j]])
    ind[[i]][[j]] <- sobol_indices(Y = Y[[i]][[j]], params = paste("X", 1:k, sep = ""), R = R, N = N, 
                                   first = first, total = i)
  }
}

# Plot sensitivity indices
lapply(ind, function(x) rbindlist(x, idcol = "Function")) %>%
  rbindlist(., idcol = "total") %>%
  .[sensitivity == "Ti"] %>%
  .[, parameters:= factor(parameters, levels = paste("X", 1:20, sep = ""))] %>%
  .[, Function:= factor(Function, levels = test_functions)] %>%
  ggplot(., aes(parameters, original, fill = total)) +
  geom_bar(stat = "identity", 
           position = position_dodge(0.5), 
           color = "black") +
  geom_errorbar(aes(ymin = low.ci,
                    ymax = high.ci),
                position = position_dodge(0.5)) +
  scale_fill_discrete(name = expression(paste("Sobol' ", T[italic(i)])),
                      labels = c("Azzini", "Homma", "Jansen", 
                                 "Monod", "Sobol")) +
  facet_grid(~Function, 
             scales = "free_x", 
             space = "free_x") +
  labs(x = "",
       y = expression(T[italic(i)])) +
  theme_AP() + 
  theme(axis.text.x = element_text(size = 6))


# CREATE METAFUNCTION ---------------------------------------------

function_list <- list(
  Linear = function(x) x,
  Quadratic = function(x) x ^ 2,
  Cubic = function(x) x ^ 3,
  Exponential = function(x) exp(1) ^ x / (exp(1) - 1),
  Periodic = function(x) sin(2 * pi * x) / 2,
  Discontinuous = function(x) ifelse(x > 0.5, 1, 0),
  Non.monotonic = function(x) 4 * (x - 0.5) ^ 2,
  Inverse = function(x) (10 - 1 / 1.1) ^ -1 * (x + 0.1) ^ - 1
)

metafunction <- function(k, n, N, coef1, coef2, coef3, f1, f2, f3, f4, f5, f6) {
  ## random sample size
  mt1 <- sobol_matrices(N = n, params = paste("X", 1:k, sep = ""), matrices = matrices)
  ## sample size of N=10000
  mt2 <- sobol_matrices(N = N, params = paste("X", 1:k, sep = ""), matrices = matrices)
  ## Compute metafunction
  output <- lapply(list(mt1, mt2), function(x) 
    # First order
    Rfast::rowsums(coef1 * function_list[[f1]](x)) + 
    # Second order
    Rfast::rowsums(coef2 * function_list[[f2]](x[, seq_len(2)]) * 
          function_list[[f3]](x[, seq_len(2)])) + 
    # Third order
    Rfast::rowsums(coef3 * function_list[[f4]](x[, seq_len(3)]) * 
          function_list[[f5]](x[, seq_len(3)]) * 
          function_list[[f6]](x[, seq_len(3)])))
  estimators <- c("jansen", "sobol", "homma", "azzini", "monod")
  ind <- Ind <- list()
  # Compute Sobol' indices
  for(i in estimators) {
    if(i == "jansen" | i == "sobol" | i == "homma") {
      first <- "saltelli"
      ind[[i]] <- lapply(list(output[[1]][1:(n * (k + 2))], 
                              output[[2]][1:(N * (k + 2))]), 
                         function(x) sobol_indices(Y = x, 
                                                   N = length(x) / (k + 2), 
                                                   params = paste("X", 1:k, sep = ""), 
                                                   first = first, 
                                                   total = i,
                                                   R = 5, 
                                                   parallel = "multicore", 
                                                   ncpus = n_cores) )
    } else if(i == "monod") {
      first <- "monod"
      ind[[i]] <- lapply(list(output[[1]][-c((n+1):(2*n))], 
                              output[[2]][-c((N+1):(2*N))]), 
                         function(x) sobol_indices(Y = x, N = length(x) / (1 + 2 * k), 
                                                   params = paste("X", 1:k, sep = ""), 
                                                   first = first, 
                                                   total = i,
                                                   R = 5, 
                                                   parallel = "multicore", 
                                                   ncpus = n_cores) )

    } else if(i == "azzini") {
      next
    }
  }
  return(ind)
}

sobol_matrices(N = 10, params = c("X", "Y", "Z"), matrices = c("A", "AB", "BA"))

# PLOT METAFUNCTION -----------------------------------------------------------

ggplot(data.frame(x = runif(100)), aes(x)) +
  map(1:length(function_list), function(nn) {
    stat_function(fun = function_list[[nn]], 
                  geom = "line", 
                  aes_(color = factor(names(function_list[nn])), 
                       linetype = factor(names(function_list[nn]))))
  }) + 
  labs(color= "Function", linetype = "Function", 
       x = expression(italic(x)), 
       y = expression(italic(y))) +
  theme_AP()

# DEFINE SETTINGS -------------------------------------------------------------

N <- 5
R <- 100 # Number of bootstrap replicas for Sobol' indices
coefficients <- paste("coef", 1:3, sep = "")
funs <- paste("f", 1:6, sep = "")
params <- c("k", "n", coefficients, funs) 
matrices <- c("A", "B", "AB", "BA")


# CREATE SAMPLE MATRIX --------------------------------------------------------

mat <- data.table(randtoolbox::sobol(N, length(params)))
setnames(mat, paste("V", 1:length(params), sep = ""), params)

mat[, k:= floor(qunif(k, 3, 100))][
  , n:= floor(qunif(n, 100, 1000))][
  , coef1:= qnorm(coef1, 0, 1)][
  , coef2:= qnorm(coef2, 10, 1)][
  , coef3:= qnorm(coef3, 3, 0.1)]

mat[, (funs):= lapply(.SD, function(x) floor(x * (8 - 1 + 1)) + 1), .SDcols = (funs)]
mat[, (funs):= lapply(.SD, function(x) ifelse(x == 1, names(function_list)[[1]],
                                              ifelse(x == 2, names(function_list)[[2]], 
                                                     ifelse(x == 3, names(function_list)[[3]], 
                                                            ifelse(x == 4, names(function_list)[[4]], 
                                                                   ifelse(x == 5, names(function_list)[[5]], 
                                                                          ifelse(x == 6, names(function_list)[[6]], 
                                                                                 ifelse(x == 7, names(function_list)[[7]], names(function_list)[[8]])))))))), 
    .SDcols = (funs)]

df <- data.frame(mat[, N:= 10000])

# RUN MODEL -------------------------------------------------------------------

# Set number of cores at 75%
n_cores <- floor(detectCores() * 0.75)

# Define parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Compute 
Y <- foreach(i=1:nrow(df)) %dopar%
  {
    metafunction(k = df[[i, "k"]], 
                 n = df[[i, "n"]], 
                 N = df[[i, "N"]],
                 coef1 = df[[i, "coef1"]], 
                 coef2 = df[[i, "coef2"]], 
                 coef3 = df[[i, "coef3"]], 
                 f1 = df[[i, "f1"]], 
                 f2 = df[[i, "f2"]], 
                 f3 = df[[i, "f3"]], 
                 f4 = df[[i, "f4"]], 
                 f5 = df[[i, "f5"]], 
                 f6 = df[[i, "f6"]])
  }

# Stop parallel cluster
stopCluster(cl)























str(Y.true)

# # EXTRACT OUTPUT ------------------------------------------------------------

estimators_AB <- c("jansen", "sobol", "homma")

# SI:
## jansen, saltelli: A, B, AB

# STI:
## jansen, sobol, homma: A, B, AB
## monod/janon: A, AB, BA
## azzini/rosati: A, B, AB, BA.

# Extract output
jansen <- monod <- list()
for(i in 1:nrow(df)) {
  # A, B and AB matrices
  jansen[[i]] <- Y[[i]][1:(df[i, "n"] * (df[i, "k"] + 2))]
}


which(do.call(rbind, lapply(A_B_AB, function(x) all(x == 0))))

# COMPUTE SOBOL' INDICES ------------------------------------------------------

ind.jansen.sobol.homma <- ind.monod <- ind.azzini <- list()
for(i in 1:nrow(df)) {
  ind.jansen.sobol.homma[[i]] <- lapply(estimators_AB, function(x) 
    sobol_indices(Y[[i]][1:(df[i, "n"] * (df[i, "k"] + 2))], # Extract A, B, AB matrices, 
                  N = df[i, "n"], 
                  params = paste("X", 1:df[i, "k"], sep = ""), 
                  total = x,
                  R = R, 
                  parallel = "multicore", 
                  ncpus = n_cores))
}













da <- sobol_indices(Y = Y[[1]], N = 550, params = paste("X", 1:51, sep = ""), R = R)

plot_sobol(da)

  ind.azzini[[i]] <- sobol_indices(Y[[i]], 
                                   N = df[i, "n"], 
                                   params = paste("X", 1:df[i, "k"], sep = ""), 
                                   first = "azzini", 
                                   total = "azzini", 
                                   R = R, 
                                   parallel = "multicore", 
                                   ncpus = n_cores)
}


# COMPUTE TRUE RANKINGS--------------------------------------------------------

df1 <- data.frame(do.call(rbind, lapply(Y, function(x) x$coefficients)))
colnames(df1) <- c("coef1", "coef2", "coef3")
df2 <- data.frame(do.call(rbind, lapply(Y, function(x) x$functions)))
colnames(df2) <- paste("f", 1:6, sep = "")

df <- cbind(mat, df1, df2)

# CREATE FUNCTION -------------------------------------------------------------

model_battle_trueranks <- function(k, n, coef1, coef2, coef3, 
                                   f1, f2, f3, f4, f5, f6) {
  mt <- sobol_matrices(N = n, params = paste("X", 1:k, sep = ""), matrices = matrices)
  out <- vector()
  for(i in 1:nrow(mt))
  out[[i]] <- sum(coef1 * function_list[[f1]](mt[i, ])) + 
    # Second order
    sum(coef2 * function_list[[f2]](mt[i, seq_len(2)]) * 
          function_list[[f3]](mt[i, seq_len(2)])) + 
    # Third order
    sum(coef3 * function_list[[f4]](mt[i, seq_len(3)]) * 
          function_list[[f5]](mt[i, seq_len(3)]) * 
          function_list[[f6]](mt[i, seq_len(3)]))
  return(out)
}

# RUN MODEL -------------------------------------------------------------------

# Set number of cores at 75%
n_cores <- floor(detectCores() * 0.75)

# Define parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Compute
Y.trueranks <- foreach(i=1:nrow(df)) %dopar%
  {
    model_battle_trueranks(k = df[[i, "k"]], 
                           n = df[[i, "n"]], 
                           coef1 = df[[i, "coef1"]], 
                           coef2 = df[[i, "coef2"]], 
                           coef3 = df[[i, "coef3"]], 
                           f1 = df[[i, "f1"]], 
                           f2 = df[[i, "f2"]], 
                           f3 = df[[i, "f3"]], 
                           f4 = df[[i, "f4"]], 
                           f5 = df[[i, "f5"]], 
                           f6 = df[[i, "f6"]])
  }
# Stop parallel cluster
stopCluster(cl)
































da <- list()
for(i in 1:nrow(prove)) {
  da[[i]] <- model_battle_trueranks(k = prove[i, "k"], 
                               N = prove[i, "N"], 
                               coef1 = prove[i, "coef1"], 
                               coef2 = prove[i, "coef2"], 
                               coef3 = prove[i, "coef3"], 
                               f1 = prove[i, "f1"], 
                               f2 = prove[i, "f2"], 
                               f3 = prove[i, "f3"], 
                               f4 = prove[i, "f4"], 
                               f5 = prove[i, "f5"], 
                               f6 = prove[i, "f6"])
}


function_list[["Exponential"]]




























plot_sobol(ind.jansen.sobol.homma[[499]][[3]])












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


