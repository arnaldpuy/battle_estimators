## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----preliminary steps, results="hide", message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------

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
loadPackages(c("Rcpp", "tidyverse", "parallel", "foreach", "doParallel", 
               "Rfast", "data.table", "scales", "cowplot", "benchmarkme", 
               "logitnorm"))

# Create custom theme
theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA))
}

# Set checkpoint

dir.create(".checkpoint")
library("checkpoint")

checkpoint("2020-01-23", 
           R.version ="3.6.1", 
           checkpointLocation = getwd())


## ----sample_matrices_functions, cache=TRUE----------------------------------------------------------------------------------------------------------------------

# FUNCTIONS TO CREATE SAMPLE MATRICES -----------------------------------------

scrambled_sobol <- function(matrices, A, B, order, cluster) {
  first <- 1:ncol(A)
  N <- nrow(A)
  if(order == "first") {
    loop <- first
  } else if(order == "second") {
    second <- c(first, utils::combn(1:ncol(A), 2, simplify = FALSE))
    loop <- second
  } else if(order == "third") {
    second <- c(first, utils::combn(1:ncol(A), 2, simplify = FALSE))
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
  final <- rbind(AB, BA)
  return(final)
}

sobol_matrices <- function(matrices = c("A", "B", "AB"), 
                           N, params, order = "first",
                           cluster = NULL) {
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

sobol_boot <- function(d, i, N, params, R, first, total, order, boot) {
  if(first == "monod" & !total == "monod" | !first == "monod" & total == "monod") {
    stop("monod should be used simultaneously in first and total indices 
         with an A, AB and BA matrices")
  }
  if(first == "azzini" & !total == "azzini" | !first == "azzini" & total == "azzini") {
    stop("azzini should be used simultaneously in first and total indices 
         with an A, B, AB and BA matrices")
  }
  k <- length(params)
  if(boot == TRUE) {
    m <- d[i, ]
  } else if(boot == FALSE) {
    m <- d
  }
  if(order == "second") {
    k <- length(params) + length(utils::combn(params, 2, simplify = FALSE))
  } else if(order == "third") {
    k <- length(params) + length(utils::combn(params, 2, simplify = FALSE)) +
      length(utils::combn(params, 3, simplify = FALSE))
  }
  # DEFINE FIRST-ORDER ESTIMATORS -----------------------
  if(first == "jansen" | first == "saltelli" &
     total == "jansen" | total == "sobol" | total == "homma") {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, -c(1, 2)]
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)
    if(first == "jansen") {
      Vi <- VY - 1 / (2 * N) * Rfast::colsums((Y_B - Y_AB) ^ 2)
      Si <- Vi[1:length(params)] / VY
    } else if(first == "saltelli") {
      Vi <- 1 / N * Rfast::colsums(Y_B * (Y_AB - Y_A))
      Si <- Vi[1:length(params)] / VY
    }
  }
  if(first == "monod" | total == "monod") {
    Y_A <- m[, 1]
    Y_AB <- m[, 2:(k + 1)]
    Y_BA <- m[,(k + 2):ncol(m)]
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA - (1 / (2 * N) * Rfast::colsums(Y_A + Y_BA)) ^ 2)
    VY <- (1 / (2 * N) * Rfast::colsums(Y_A ^ 2 + Y_BA ^ 2) - 
             (1 / (2 * N) * Rfast::colsums(Y_A + Y_BA)) ^ 2)
    Si <- Vi[1:length(params)] / VY[1:length(params)]
    STi <- 1 - (1 / N * Rfast::colsums(Y_A * Y_AB) - 
                  (1/ N * Rfast::colsums((Y_A + Y_AB) / 2)) ^ 2) / 
      (1 / N * Rfast::colsums((Y_A ^ 2 + Y_AB ^ 2) / 2) - 
         (1/ N * Rfast::colsums((Y_A + Y_AB) / 2)) ^ 2)
    STi <- STi[1:length(params)]
  }
  if(first == "azzini" | total == "azzini") {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, 3:(3 + k - 1)]
    Y_BA <- m[, (ncol(m) - k + 1):ncol(m)]
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA) - ((1 / N) * sum(Y_A * Y_B)) + 
      (1 / N) * Rfast::colsums(Y_B * Y_AB) - ((1 / N) * Rfast::colsums(Y_AB * Y_BA))
    VY <- (1 / (2 * N ) * (sum(Y_B * (Y_B - Y_A)) + sum(Y_A * (Y_A - Y_B)) + 
                             Rfast::colsums(Y_BA * (Y_BA - Y_AB)) + 
                             Rfast::colsums(Y_AB * (Y_AB - Y_BA))))
    Si <- Vi[1:length(params)] / VY[1:length(params)]
    STi <- 1 - abs(Rfast::colsums((Y_A - Y_BA) * (Y_B - Y_AB)) / 
                     (1 / 2 * Rfast::colsums((Y_A - Y_B) ^ 2 + (Y_AB - Y_BA) ^ 2)))
    STi <- STi[1:length(params)]
  }
  # DEFINE TOTAL-ORDER ESTIMATORS FOR THE REST ----------------
  if(total == "jansen") {
    STi <- (1 / (2 * N) * Rfast::colsums((Y_A - Y_AB) ^ 2)) / VY
  } else if(total == "homma") {
    STi <- (VY - (1 / N) * Rfast::colsums(Y_A * Y_AB) + f0 ^ 2) / VY
  } else if(total == "sobol") {
    STi <- ((1 / N) * Rfast::colsums(Y_A * (Y_A - Y_AB))) / VY
  }
  STi <- STi[1:length(params)]
  if(order == "second" | order == "third") {
    com2 <- utils::combn(1:length(params), 2, simplify = FALSE)
    mat2 <- t(mapply(c, Vi[(length(params) + 1):(length(params) + length(com2))], lapply(com2, function(x) Vi[x])))
    Vij <- apply(mat2, 1, function(x) Reduce("-", x))
    if(first == "monod" | first == "azzini") {
      Sij <- Vij / VY[(length(params) + 1):(length(VY) - length(utils::combn(params, 3, simplify = FALSE)))]
    } else {
      Sij <- Vij / VY
    }
  }
  if(order == "first") {
    Sij <- NULL
  }
  if(order == "third") {
    tmp <- do.call(rbind, com2)
    Vij.vec <- as.numeric(paste(tmp[, 1], tmp[, 2], sep = ""))
    Vij.named <- Vij
    names(Vij.named) <- Vij.vec
    com3 <- utils::combn(1:length(params), 3, simplify = FALSE)
    Vi.only <- do.call(rbind, lapply(com3, function(x) Vi[x])) # Extract Vi, Vj, Vk
    Vijk.only <- tail(Vi, length(com3)) # Extract Vijk
    tmp3 <- do.call(rbind, com3)
    first.pairwise <- lapply(paste(tmp3[, 1], tmp3[, 2], sep = ""), function(x) Vij.named[x])
    second.pairwise <- lapply(paste(tmp3[, 1], tmp3[, 3], sep = ""), function(x) Vij.named[x]) 
    third.pairwise <- lapply(paste(tmp3[, 2], tmp3[, 3], sep = ""), function(x) Vij.named[x])
    Vij.only <- t(mapply(cbind, first.pairwise, second.pairwise, third.pairwise))
    mat3 <- cbind(Vijk.only, Vij.only, Vi.only)
    Vijk <- apply(mat3, 1, function(x) Reduce("-", x))
    if(first =="monod" | first == "azzini") {
      Sijk <- Vijk / tail(VY, length(utils::combn(params, 3, simplify = FALSE)))
    } else {
      Sijk <- Vijk / VY
    }
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
sobol_indices <- function(Y, N, params, R = NULL, first = "saltelli", total = "jansen", 
                          order = "first", parallel = "no", boot = FALSE,
                          ncpus = 1, conf = 0.95, type = "norm") {
  k <- length(params)
  d <- matrix(Y, nrow = N)
  if(boot == FALSE) {
    tmp <- sobol_boot(d = d, N = N, params = params, first = first, total = total, 
                      order = order, boot = FALSE)
    out <- data.table::data.table(tmp)
    data.table::setnames(out, "tmp", "original")
  } else if(boot == TRUE) {
    tmp <- boot::boot(data = d, statistic = sobol_boot, R = R, N = N, params = params, 
                      first = first, total = total, order = order, 
                      parallel = parallel, ncpus = ncpus, boot = TRUE)
    out <- bootstats(tmp, conf = conf, type = type)
  } else {
    stop("boot has to be TRUE or FALSE")
  }
  if(order == "first") {
    parameters <- c(rep(params, times = 2))
    indices <- c(rep(c("Si", "Ti"), each = k))
  } else if(order == "second") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x) 
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    indices <- c(rep(c("Si", "Ti"), each = length(params)), 
                 rep("Sij", times = length(vector.second)))
  } else if(order == "third") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x) 
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    vector.third <- unlist(lapply(utils::combn(params, 3, simplify = FALSE), function(x) 
      paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.third)
    indices <- c(rep(c("Si", "Ti"), each = k), 
                 rep("Sij", times = length(vector.second)), 
                 rep("Sijk", times = length(vector.third)))
  }
  out[, sensitivity:= cbind(indices)][, parameters:= cbind(parameters)]
  return(out)
}

# FUNCTION TO PLOT SOBOL' INDICES ----------------------------------------------

plot_sobol <- function (x, dummy = NULL, type = 1) {
  sensitivity <- low.ci <- high.ci <- parameters <- original <- NULL
  if (type == 1) {
    if (is.null(dummy) == FALSE) {
      plot.dummy <- geom_rect(data = dummy, aes(ymin = 0, 
                                                ymax = high.ci, xmin = -Inf, xmax = Inf, fill = sensitivity), 
                              alpha = 0.2, inherit.aes = FALSE)
    }
    else {
      plot.dummy <- NULL
    }
    p <- x[sensitivity == "Si" | sensitivity == "Ti"]
    gg <- ggplot2::ggplot(p, aes(parameters, original, fill = sensitivity)) + 
      geom_bar(stat = "identity", position = position_dodge(0.6), 
               color = "black") + plot.dummy + geom_errorbar(aes(ymin = low.ci, 
                                                                 ymax = high.ci), position = position_dodge(0.6)) + 
      scale_fill_discrete(name = "Sobol' indices", labels = c(expression(S[italic(i)]), 
                                                              expression(S[italic(T[i])]))) + 
      labs(x = "", 
           y = "Variance") + 
      theme_bw() + 
      theme(legend.position = "top", 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.background = element_rect(fill = "transparent", 
                                             color = NA), 
            legend.key = element_rect(fill = "transparent", 
                                      color = NA))
  }
  else if (!type == 1) {
    if (type == 2) {
      plot.type <- "Sij"
    }
    else if (type == 3) {
      plot.type <- "Sijk"
    }
    else {
      stop("Type should be either 1, 2 or 3")
    }
    p <- x[sensitivity == plot.type]
    gg <- ggplot2::ggplot(p, aes(stats::reorder(parameters, 
                                                original), original)) + 
      geom_point() + 
      geom_errorbar(aes(ymin = low.ci, 
                        ymax = high.ci)) + 
      theme_bw() + labs(x = "", y = "Variance") + 
      geom_hline(yintercept = 0, lty = 2, color = "red") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.background = element_rect(fill = "transparent", 
                                             color = NA), 
            legend.key = element_rect(fill = "transparent", 
                                      color = NA), 
            axis.text.x = element_text(angle = 45, 
                                       hjust = 1))
  }
  return(gg)
}

## ----savage_scores, cache=TRUE----------------------------------------------------------------------------------------------------------------------------------

# SAVAGE SCORES ---------------------------------------------------------------

savage_scores <- function(x) {
  true.ranks <- rank(-x)
  p <- sort(1 / true.ranks)
  mat <- matrix(rep(p, length(p)), nrow = length(p), byrow = TRUE)
  mat[upper.tri(mat)] <- 0
  out <- sort(rowSums(mat), decreasing = TRUE)[true.ranks]
  return(out)
}


## ----ti_indices, cache=TRUE, dependson=c("savage_scores", "sample_matrices_functions")--------------------------------------------------------------------------

# COMPUTATION OF SOBOL' Ti INDICES --------------------------------------------

sobol_Ti <- function(d, N, params, total) {
  m <- matrix(d, nrow = N)
  k <- length(params)
  if(!total == "azzini") {
    Y_A <- m[, 1]
    Y_AB <- m[, -1]
    f0 <- (1 / length(Y_A)) * sum(Y_A)
    VY <- 1 / length(Y_A) * sum((Y_A - f0) ^ 2)
     # VY <- 1 / length(Y_A) * (sum(Y_A ^ 2) - 
    # (1 / N * sum(Y_A ^ 2))) ((Variance used by Becker))
  }
  if(total == "jansen") {
    Ti <- (1 / (2 * N) * Rfast::colsums((Y_A - Y_AB) ^ 2)) / VY
  } else if(total == "homma") {
    Ti <- (VY - (1 / N) * Rfast::colsums(Y_A * Y_AB) + f0 ^ 2) / VY
  } else if(total == "sobol") {
    Ti <- ((1 / N) * Rfast::colsums(Y_A * (Y_A - Y_AB))) / VY
  } else if(total == "monod") {
    Ti <- 1 - (1 / N * Rfast::colsums(Y_A * Y_AB) - 
                  (1/ N * Rfast::colsums((Y_A + Y_AB) / 2)) ^ 2) / 
      (1 / N * Rfast::colsums((Y_A ^ 2 + Y_AB ^ 2) / 2) - 
         (1/ N * Rfast::colsums((Y_A + Y_AB) / 2)) ^ 2)
  }
  if(total == "azzini") {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, 3:(3 + k - 1)]
    Y_BA <- m[, (ncol(m) - k + 1):ncol(m)]
    Ti <- 1 - abs(Rfast::colsums((Y_A - Y_BA) * (Y_B - Y_AB)) / 
                     (1 / 2 * Rfast::colsums((Y_A - Y_B) ^ 2 + (Y_AB - Y_BA) ^ 2)))
  } 
  output <- data.table(Ti)
  output[, `:=`(parameters = paste("X", 1:k, sep = ""), 
                ranks= rank(-Ti), 
                savage.scores = savage_scores(Ti))]
  return(output)
}


## ----check_ti, cache=TRUE, dependson="ti_indices"---------------------------------------------------------------------------------------------------------------

# CHECK THAT ALL TI ESTIMATORS WORK -------------------------------------------

# Settings
estimators <- c("jansen", "sobol", "homma", "azzini", "monod")
test_functions <- c("Ishigami", "Sobol'G", "Morris")
N <- 2^9

# Run model
ind <- Y <- mt <- list() 
for(i in estimators) {
  for(j in test_functions) {
    if(!i == "azzini") {
      matrices <- c("A", "AB")
    } else {
      matrices <- c("A", "B", "AB", "BA")
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
    ind[[i]][[j]] <- sobol_Ti(d = Y[[i]][[j]], params = paste("X", 1:k, sep = ""), 
                              N = N, total = i)
  }
}


## ----plot_prove, cache=TRUE, dependson="check_ti", dev="tikz", fig.height=3, fig.width=6.5----------------------------------------------------------------------

# PLOT SENSITIVITY INDICES ----------------------------------------------------

lapply(ind, function(x) rbindlist(x, idcol = "Function")) %>%
  rbindlist(., idcol = "estimator") %>%
  .[, parameters:= factor(parameters, levels = paste("X", 1:20, sep = ""))] %>%
  .[, Function:= factor(Function, levels = test_functions)] %>%
  ggplot(., aes(parameters, Ti, fill = estimator)) +
  geom_bar(stat = "identity", 
           position = position_dodge(0.7), 
           color = "black") +
  facet_grid(~Function, 
             scales = "free_x", 
             space = "free_x") +
  scale_fill_discrete(name = expression(paste("Sobol' ", T[italic(i)])),
                      labels = c("Azzini", "Homma", "Jansen", 
                                 "Monod", "Sobol")) +
  labs(x = "",
       y = expression(T[italic(i)])) +
  theme_AP() + 
  theme(axis.text.x = element_text(size = 6.5), 
        legend.position = "top")


## ----functions_metafunction, cache=TRUE-------------------------------------------------------------------------------------------------------------------------

# CREATE METAFUNCTION ---------------------------------------------------------

function_list <- list(
  Linear = function(x) x,
  Quadratic = function(x) x ^ 2,
  Cubic = function(x) x ^ 3,
  Exponential = function(x) exp(1) ^ x / (exp(1) - 1),
  Periodic = function(x) sin(2 * pi * x) / 2,
  Discontinuous = function(x) ifelse(x > 0.5, 1, 0),
  Non.monotonic = function(x) 4 * (x - 0.5) ^ 2,
  Inverse = function(x) (10 - 1 / 1.1) ^ -1 * (x + 0.1) ^ - 1, 
  No.effect = function(x) x * 0, 
  Trigonometric = function(x) cos(x)
)


## ----plot_functions_metafunction, cache=TRUE, dependson="functions_metafunction", dev = "tikz", fig.height=2.7, fig.width=4.6, fig.cap="Functions used in the metafunction of @Becker2019."----

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
  theme_AP() + 
  theme(legend.position = "right")


## ----source_cpp, warning = FALSE--------------------------------------------------------------------------------------------------------------------------------

Rcpp::sourceCpp("vector_multiplication.cpp")


## ----metafunction, cache=TRUE, dependson = "functions_metafunction", warning=FALSE------------------------------------------------------------------------------

# DEFINE METAFUNCTION ---------------------------------------------------------

metafunction <- function(X, k_2, k_3, epsilon) {
  k <- ncol(X)
  # Define functions according to k
  set.seed(epsilon) # Set seed before every call to a random number generator
  all_functions <- sample(names(function_list), k, replace = TRUE)
  # Define coefficients
  set.seed(epsilon) # Set seed before every call to a random number generator
  components <- sample(1:2, prob = c(0.5, 0.5), size = 200, replace = TRUE)
  mus <- c(0, 0)
  sds <- sqrt(c(0.5, 5))
  coefficients <- rnorm(100) * sds[components] + mus[components]
  # Coefficients for first
  set.seed(epsilon) # Set seed before every call to a random number generator
  coefD1 <- sample(coefficients, k)
  # Coefficients for pairs
  d2 <- t(utils::combn(1:k, 2))
  # Define functions according to k
  set.seed(epsilon) # Set seed before every call to a random number generator
  d2M <- d2[sample(nrow(d2), size = ceiling(k * k_2), replace = FALSE), ]
  sample.size.d2M <- ifelse(is.vector(d2M) == TRUE, 1, nrow(d2M))
  # Define functions according to k
  set.seed(epsilon) # Set seed before every call to a random number generator
  coefD2 <- sample(coefficients, sample.size.d2M, replace = TRUE)
  # Coefficients for triplets
  d3 <- t(utils::combn(1:k, 3))
  # Define functions according to k
  set.seed(epsilon) # Set seed before every call to a random number generator
  d3M <- d3[sample(nrow(d3), size = ceiling(k * k_3), replace = FALSE), ]
  sample.size <- ifelse(is.vector(d3M) == TRUE, 1, nrow(d3M))
  # Define functions according to k
  set.seed(epsilon) # Set seed before every call to a random number generator
  coefD3 <- sample(coefficients, sample.size, replace = TRUE)
  # Run sampled functions in each column
  output <- sapply(seq_along(all_functions), function(x) function_list[[all_functions[x]]](X[, x]))
  y1 <- Rfast::rowsums(mmult(output, coefD1))
  if(is.vector(d2M) == TRUE) {
    y2 <- sum(output[, d2M[1]] *  output[, d2M[2]] * coefD2)
  } else {
    y2 <- Rfast::rowsums(mmult(output[, d2M[, 1]] *  output[, d2M[, 2]], coefD2))
  }
  if(is.vector(d3M) == TRUE) {
    y3 <- sum(output[, d3M[1]] *  output[, d3M[2]] * output[, d3M[3]] * coefD3)
  } else {
    y3 <- Rfast::rowsums(mmult(output[, d3M[, 1]] *  output[, d3M[, 2]] * output[, d3M[, 3]], coefD3))
  }
  Y <- y1 + y2 + y3
}

# CREATE FUNCTION FOR RANDOM DISTRIBUTIONS -----------------------------------

sample_distributions <- list(
  "uniform" = function(x) x,
  "normal" = function(x) qnorm(x, 0.5, 0.2),
  "beta" = function(x) qbeta(x, 8, 2),
  "beta2" = function(x) qbeta(x, 2, 8),
  "beta3" = function(x) qbeta(x, 2, 0.5),
  "beta4" = function(x) qbeta(x, 0.5, 2),
  "logitnormal" = function(x) qlogitnorm(x, 0, 3.16)
  # Logit-normal, Bates too?
)

random_distributions <- function(X, phi) {
  names_ff <- names(sample_distributions)
  if(!phi == length(names_ff) + 1) {
    out <- sample_distributions[[names_ff[phi]]](X)
  } else {
    temp <- sample(names_ff, ncol(X), replace = TRUE)
    out <- sapply(seq_along(temp), function(x) sample_distributions[[temp[x]]](X[, x]))
  }
  return(out)
}

# PLOT DISTRIBUTIONS ------------------------------------------------------

names_ff <- names(sample_distributions)
prove <- randtoolbox::sobol(n = 1000, dim = length(names_ff))

out <- data.table(sapply(seq_along(names_ff), function(x) 
  sample_distributions[[names_ff[x]]](prove[, x])))

melt(out) %>%
  ggplot(., aes(value, group = variable, colour = variable)) + 
  geom_density() + 
  scale_color_discrete(labels = c("U(0, 1)", 
                                  "N(0.5, 0.2)", 
                                  "Beta(8, 2)", 
                                  "Beta(2, 8)", 
                                  "Beta(2, 0.5)", 
                                  "Beta(0.5, 2)", 
                                  "Logitnormal(0, 3.16)"), 
                       name = "Distribution") +
  labs(x = expression(italic(x)), 
       y = "Density") +
  theme_AP() 

## ----settings, cache=TRUE---------------------------------------------------------------------------------------------------------------------------------------

# DEFINE SETTINGS -------------------------------------------------------------

N <- 2 ^ 8 # Sample size of sample matrix
R <- 100 # Number of bootstrap replicas
n_cores <- detectCores() * 0.75
order <- "first"
params <- c("k", "C", "k_2", "k_3", "epsilon", "phi") 
N.high <- 2 ^ 12 # Maximum sample size of the large sample matrix


## ----sample_matrix, cache=TRUE, dependson="settings"------------------------------------------------------------------------------------------------------------

# CREATE SAMPLE MATRIX --------------------------------------------------------

mat <- sobol_matrices(N = N, params = params, order = order)
mat[, 1] <- floor(qunif(mat[, 1], 3, 200)) # k
mat[, 2] <- floor(qunif(mat[, 2], 10, 2000)) # C
mat[, 3] <- round(qunif(mat[, 3], 0.3, 0.5), 2) # k_2
mat[, 4] <- round(qunif(mat[, 4], 0.1, 0.3), 2) # k_3
mat[, 5] <- floor(qunif(mat[, 5], 1, 200)) # Epsilon
mat[, 6] <- floor(mat[, 6] * (6 - 1 + 1)) + 1 # Phi

colnames(mat) <- params

N.all <- apply(mat, 1, function(x) ceiling(x["C"] / (x["k"] + 1)))
N.azzini <- apply(mat, 1, function(x) ceiling(x["C"] / (2 * x["k"] + 2)))

tmp <- cbind(mat, N.all, N.azzini)
sel <- c("N.all", "N.azzini")

mat <- as.matrix(data.table(tmp)[, (sel):= lapply(.SD, function(x) 
  ifelse(x == 1, 2, x)), .SDcols = (sel)])


## ----define_model, cache=TRUE, dependson=c("sample_matrices_functions", "metafunction", "ti_indices", "savage_scores")------------------------------------------

# DEFINE MODEL ----------------------------------------------------------------

model_Ti <- function(k, N.all, N.azzini, N.high, k_2, k_3, epsilon, phi) {
  ind <- list()
  estimators <- c("jansen", "sobol", "homma", "monod", "azzini")
  all.but.azzini <- sobol_matrices(N = N.all, params = paste("X", 1:k, sep = ""), 
                                   matrices = c("A", "AB"))
  azzini <- sobol_matrices(N = N.azzini, params = paste("X", 1:k, sep = ""), 
                           matrices = c("A", "B", "AB", "BA"))
  large.matrix <- sobol_matrices(N = N.high, params = paste("X", 1:k, sep = ""), 
                                 matrices = c("A", "AB"))
  set.seed(epsilon)
  all.matrices <- random_distributions(X = rbind(all.but.azzini, azzini, large.matrix), 
                                       phi = phi)
  output <- metafunction(X = all.matrices, 
                         k_2 = k_2, 
                         k_3 = k_3, 
                         epsilon = epsilon)
  full.ind <- sobol_Ti(d = tail(output, nrow(large.matrix)), 
                       N = N.high, 
                       params = paste("X", 1:k, sep = ""), 
                       total = "jansen")
  full.ind[, sample.size:= "N"]
  for(i in estimators) {
    if(!i == "azzini") {
      y <- output[1:nrow(all.but.azzini)]
      n <- N.all
    } else {
      y <- output[-c(1:nrow(all.but.azzini), 
                     (nrow(all.but.azzini) + nrow(azzini) + 1):length(output))]
      n <- N.azzini
    }
    ind[[i]] <- sobol_Ti(d = y, N = n, params = paste("X", 1:k, sep = ""), total = i)
    ind[[i]][, sample.size:= "n"]
    ind[[i]] <- rbind(ind[[i]], full.ind)
  }
  return(ind)
}


## ----model_run, cache=TRUE, dependson=c("define_model", "settings", "sample_matrix", "source_cpp")--------------------------------------------------------------

# RUN MODEL -------------------------------------------------------------------

Y.ti <- list()
for(i in 1:nrow(mat)) {
  Y.ti[[i]] <- model_Ti(k = mat[[i, "k"]], 
                        k_2 = mat[[i, "k_2"]], 
                        k_3 = mat[[i, "k_3"]],
                        epsilon = mat[[i, "epsilon"]],
                        phi = mat[[i, "phi"]],
                        N.all = mat[[i, "N.all"]], 
                        N.azzini = mat[[i, "N.azzini"]], 
                        N.high = N.high)
}

## ----arrange_output, cache=TRUE, dependson="model_run"----------------------------------------------------------------------------------------------------------

# ARRANGE OUTPUT --------------------------------------------------------------

out <- lapply(Y.ti, function(x) rbindlist(x, idcol = "estimator")) %>%
  rbindlist(., idcol = "row")

out_wide <- spread(out[, .(sample.size, savage.scores, parameters, estimator, row)], 
                   sample.size, savage.scores)

out_cor <- out_wide[, .(correlation = cor(N, n)), .(estimator, row)]
mt.dt <- data.table(mat) %>%
  .[, row:= 1:.N]

full_output <- merge(mt.dt, out_cor) %>%
  .[, Nt:= ifelse(estimator == "azzini", N.azzini * (2 * k + 2), N.all * (k + 1))] %>%
  .[, estimator:= ifelse(estimator %in% "azzini", "Azzini and Rosati", 
                        ifelse(estimator %in% "homma", "Homma and Saltelli",
                               ifelse(estimator %in% "monod", "Janon/Monod", 
                                      ifelse(estimator %in% "jansen", "Jansen", "Sobol'"))))] %>%
  .[, ratio:= Nt / k]

# Define A matrix
A <- full_output[,.SD[1:N], estimator]

## ----export_output, cache=TRUE, dependson="arrange_output"------------------------------------------------------------------------------------------------------

# EXPORT OUTPUT ---------------------------------------------------------------

fwrite(out, "out.csv")


## ----plot_full, cache=TRUE, dependson="arrange_output", dev = "tikz", fig.height=8, fig.width=4-----------------------------------------------------------------

# PLOT OUTPUT -----------------------------------------------------------------

# Scatterplot
a <- ggplot(A[correlation > 0], aes(Nt, k, color = correlation)) +
  geom_point(size = 0.6) + 
  scale_colour_gradientn(colours = c("purple", "red", "orange", "lightgreen"), 
                         name = expression(italic(r))) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = expression(italic(N[t])), 
       y = expression(italic(k))) + 
  facet_wrap(~estimator, 
             ncol = 1) + 
  theme_AP() + 
  theme(legend.position = "none")

# Get legend
legend <- get_legend(a + theme(legend.position = "top"))

# Ratio
b <- ggplot(A[correlation > 0], aes(ratio, correlation)) +
  geom_point(alpha = 0.15, size = 0.6) +
  facet_wrap(~estimator, 
             ncol = 1) +
  labs(x = expression(italic(N[t]/k)), 
       y = expression(italic(r))) +
  scale_x_log10() +
  theme_AP()

# Merge plot
bottom <- plot_grid(a, b, ncol = 2, labels = "auto")
plot_grid(legend, bottom, ncol = 1, rel_heights = c(0.15, 1))


## ----plot_boxplot, cache=TRUE, dependson="arrange_output", dev = "tikz", fig.width=4, fig.height=3--------------------------------------------------------------

# PLOT BOXPLOT ----------------------------------------------------------------

ggplot(A[correlation > 0], aes(estimator, correlation)) +
  geom_boxplot() + 
  labs(x = "", 
       y = expression(italic(r))) + 
  theme_AP() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# SENSITIVITY ANALYSIS --------------------------------------------------------

indices <- full_output[, sobol_indices(Y = correlation, 
                                       N = N,
                                       params = params,
                                       R = R, 
                                       boot = TRUE, 
                                       order = order), 
                       estimator]

# PLOT SOBOL' INDICES --------------------------------------------------------

plot_sobol(indices) +
  facet_wrap(~estimator, 
             ncol = 1)

## ----session_information----------------------------------------------------------------------------------------------------------------------------------------

# SESSION INFORMATION ---------------------------------------------------------

sessionInfo()

## Return the machine CPU
cat("Machine:     "); print(get_cpu()$model_name)

## Return number of true cores
cat("Num cores:   "); print(detectCores(logical = FALSE))

## Return number of threads
cat("Num threads: "); print(detectCores(logical = TRUE))

## Return the machine RAM
cat("RAM:         "); print (get_ram()); cat("\n")

