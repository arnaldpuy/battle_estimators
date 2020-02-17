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
               "Rfast", "data.table", "scales", "cowplot", "benchmarkme"))

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

metafunction <- function(X) {
  k <- ncol(X)
  # Define functions according to k
  all_functions <- sample(names(function_list), k, replace = TRUE)
  # Define coefficients
  components <- sample(1:2, prob = c(0.5, 0.5), size = 200, replace = TRUE)
  mus <- c(0, 0)
  sds <- sqrt(c(0.5, 5))
  coefficients <- rnorm(100) * sds[components] + mus[components]
  # Coefficients for first
  coefD1 <- sample(coefficients, k)
  # Coefficients for pairs
  d2 <- t(utils::combn(1:k, 2))
  d2M <- d2[sample(nrow(d2), size = ceiling(k * 0.5), replace = FALSE), ]
  coefD2 <- sample(coefficients, nrow(d2M), replace = TRUE)
  # Coefficients for triplets
  d3 <- t(utils::combn(1:k, 3))
  d3M <- d3[sample(nrow(d3), size = ceiling(k * 0.2), replace = FALSE), ]
  sample.size <- ifelse(is.vector(d3M) == TRUE, 1, nrow(d3M))
  coefD3 <- sample(coefficients, sample.size, replace = TRUE)
  # Run sampled functions in each column
  output <- sapply(seq_along(all_functions), function(x) function_list[[all_functions[x]]](X[, x]))
  y1 <- Rfast::rowsums(mmult(output, coefD1))
  y2 <- Rfast::rowsums(mmult(output[, d2M[, 1]] *  output[, d2M[, 2]], coefD1))
  if(is.vector(d3M) == TRUE) {
    y3 <- sum(output[, d3M[1]] *  output[, d3M[2]] * output[, d3M[3]] * coefD3)
  } else {
    y3 <- Rfast::rowsums(mmult(output[, d3M[, 1]] *  output[, d3M[, 2]] * output[, d3M[, 3]], coefD3))
  }
  Y <- y1 + y2 + y3
  return(Y)
}


## ----settings, cache=TRUE---------------------------------------------------------------------------------------------------------------------------------------

# DEFINE SETTINGS -------------------------------------------------------------

N <- 2 ^ 10 # Sample size of sample matrix
params <- c("k", "C") 
N.high <- 2 ^ 13 # Maximum sample size of the large sample matrix


## ----sample_matrix, cache=TRUE, dependson="settings"------------------------------------------------------------------------------------------------------------

# CREATE SAMPLE MATRIX --------------------------------------------------------

mat <- randtoolbox::sobol(N, length(params))
mat[, 1] <- floor(qunif(mat[, 1], 3, 200))
mat[, 2] <- floor(qunif(mat[, 2], 10, 2000))
colnames(mat) <- params

N.all <- apply(mat, 1, function(x) ceiling(x["C"] / (x["k"] + 1)))
N.azzini <- apply(mat, 1, function(x) ceiling(x["C"] / (2 * x["k"] + 2)))

tmp <- cbind(mat, N.all, N.azzini)
sel <- c("N.all", "N.azzini")

mat <- as.matrix(data.table(tmp)[, (sel):= lapply(.SD, function(x) 
  ifelse(x == 1, 2, x)), .SDcols = (sel)])


## ----define_model, cache=TRUE, dependson=c("sample_matrices_functions", "metafunction", "ti_indices", "savage_scores")------------------------------------------

# DEFINE MODEL ----------------------------------------------------------------

model_Ti <- function(k, N.all, N.azzini, N.high) {
  ind <- list()
  estimators <- c("jansen", "sobol", "homma", "monod", "azzini")
  all.but.azzini <- sobol_matrices(N = N.all, params = paste("X", 1:k, sep = ""), 
                                   matrices = c("A", "AB"))
  azzini <- sobol_matrices(N = N.azzini, params = paste("X", 1:k, sep = ""), 
                           matrices = c("A", "B", "AB", "BA"))
  large.matrix <- sobol_matrices(N = N.high, params = paste("X", 1:k, sep = ""), 
                                 matrices = c("A", "AB"))
  output <- metafunction(rbind(all.but.azzini, azzini, large.matrix))
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
                                      ifelse(estimator %in% "jansen", "Jansen", "Sobol'"))))]


## ----export_output, cache=TRUE, dependson="arrange_output"------------------------------------------------------------------------------------------------------

# EXPORT OUTPUT ---------------------------------------------------------------

fwrite(out, "out.csv")


## ----plot_full, cache=TRUE, dependson="arrange_output", dev = "tikz", fig.height=8, fig.width=4-----------------------------------------------------------------

# PLOT OUTPUT -----------------------------------------------------------------

full.output <- full_output[, ratio:= Nt / k]

# Scatterplot
a <- ggplot(full_output[correlation > 0], aes(Nt, k, color = correlation)) +
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
b <- ggplot(full_output[correlation > 0], aes(ratio, correlation)) +
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

ggplot(full_output[correlation > 0], aes(estimator, correlation)) +
  geom_boxplot() + 
  labs(x = "", 
       y = expression(italic(r))) + 
  theme_AP() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


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

