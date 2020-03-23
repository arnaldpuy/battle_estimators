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

install.packages("devtools") # if you have not installed devtools package already
devtools::install_github("arnaldpuy/sensobol")

# Load the packages
loadPackages(c("Rcpp", "RcppArmadillo", "tidyverse", "parallel", "foreach", 
               "doParallel", "Rfast", "data.table", "scales", "cowplot", 
               "benchmarkme", "logitnorm", "sensobol"))

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
  output[, `:=`(parameters = paste("X", 1:k, sep = ""))]
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
      modelRun <- sensobol::ishigami_Fun
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

a <- ggplot(data.frame(x = runif(100)), aes(x)) +
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

a

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

# PLOT DISTRIBUTIONS ----------------------------------------------------------

names_ff <- names(sample_distributions)
prove <- randtoolbox::sobol(n = 1000, dim = length(names_ff))

out <- data.table(sapply(seq_along(names_ff), function(x) 
  sample_distributions[[names_ff[x]]](prove[, x])))

b <- data.table::melt(out) %>%
  ggplot(., aes(value, group = variable, colour = variable)) + 
  geom_density() + 
  scale_color_discrete(labels = c("$U(0, 1)$", 
                                  "$\\mathcal{N}(0.5, 0.2)$", 
                                  "$Beta(8, 2)$", 
                                  "$Beta(2, 8)$", 
                                  "$Beta(2, 0.5)$", 
                                  "$Beta(0.5, 2)$", 
                                  "$Logitnormal(0, 3.16)$"), 
                       name = "") +
  labs(x = expression(italic(x)), 
       y = "Density") +
  theme_AP() 

b

# MERGE METAFUNCTION PLOT AND DISTRIBUTIONS PLOT -------------------------------

plot_grid(a, b, ncol = 1, labels = "auto", align = "hv")

# EXEMPLIFY THE METAFUNCTION WITH AN EXAMPLE ----------------------------------

N <- 5000 # Sample size
R <- 100 # Number of bootstrap replications
k <- 55 # Number of model inputs
k_2 <- 0.5 # Fraction of active pairwise interactions
k_3 <- 0.3 # Fraction of active three-wise interactions
epsilon <- 5 # to reproduce the results
params <- paste("X", 1:k, sep = "")
A <- sobol_matrices(N = N, params = params)

# Compute
Y <- metafunction(data = A, k_2 = k_2, k_3 = k_3, epsilon = epsilon)
ind <- sobol_indices(Y = Y, N = N, params = params, R = R, boot = TRUE)
ind <- ind[, parameters:= factor(parameters, levels = paste("X", 1:k, sep = ""))]

plot_sobol(ind) + 
  theme(axis.text.x = element_text(size = 4))

## ----settings, cache=TRUE---------------------------------------------------------------------------------------------------------------------------------------

# DEFINE SETTINGS -------------------------------------------------------------

N <- 2 ^ 4 # Sample size of sample matrix
R <- 100 # Number of bootstrap replicas
n_cores <- ceiling(detectCores() * 0.5)
order <- "first"
params <- c("k", "N_t", "k_2", "k_3", "epsilon", "phi", "delta") 
N.high <- 2 ^ 10 # Maximum sample size of the large sample matrix


## ----sample_matrix, cache=TRUE, dependson="settings"------------------------------------------------------------------------------------------------------------

# CREATE SAMPLE MATRIX --------------------------------------------------------

mat <- sobol_matrices(N = N, params = params, order = order)
mat[, 1] <- floor(qunif(mat[, 1], 3, 100)) # k
mat[, 2] <- floor(qunif(mat[, 2], 10, 1000)) # N_t
mat[, 3] <- round(qunif(mat[, 3], 0.3, 0.5), 2) # k_2
mat[, 4] <- round(qunif(mat[, 4], 0.1, 0.3), 2) # k_3
mat[, 5] <- floor(qunif(mat[, 5], 1, 200)) # Epsilon
mat[, 6] <- floor(mat[, 6] * (8 - 1 + 1)) + 1 # Phi
mat[, 7] <- floor(mat[, 7] * (3 - 1 + 1)) + 1 # Phi

colnames(mat) <- params

N.all <- apply(mat, 1, function(x) ceiling(x["N_t"] / (x["k"] + 1)))
N.azzini <- apply(mat, 1, function(x) ceiling(x["N_t"] / (2 * x["k"] + 2)))

tmp <- cbind(mat, N.all, N.azzini)
sel <- c("N.all", "N.azzini")

mat <- as.matrix(data.table(tmp)[, (sel):= lapply(.SD, function(x) 
  ifelse(x == 1, 2, x)), .SDcols = (sel)])

da <- data.table(mat)[, `:=` (Nt.all = N.all * (k + 1), 
                        Nt.azzini = N.azzini * (2 * k + 2))]


melt(da, measure.vars = c("Nt.all", "Nt.azzini")) %>%
  ggplot(., aes(value, k, color = variable)) +
  geom_line()


## ----define_model, cache=TRUE, dependson=c("sample_matrices_functions", "metafunction", "ti_indices", "savage_scores")------------------------------------------

# DEFINE MODEL ----------------------------------------------------------------

model_Ti <- function(k, N.all, N.azzini, N.high, k_2, k_3, epsilon, phi, delta) {
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
  output <- sensobol::metafunction(data = all.matrices, 
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
  # Arrange data
  out <- rbindlist(ind, idcol = "estimator") 
  out.wide <- dcast(out, estimator + parameters ~ sample.size, value.var = "Ti")
  if(delta == 1) { # Regular Pear
    final <- out.wide[, .(correlation = cor(N, n)), estimator]
  } else if(delta == 2) { # kendall tau
    final <- out.wide[, .(correlation = pcaPP::cor.fk(N, n)), estimator]
  } else { # Savage ranks
    final <- out.wide[, lapply(.SD, savage_scores), .SDcols = c("N", "n"), estimator][
      , .(correlation = cor(N, n)), estimator]
  }
  return(final)
}


## ----model_run, cache=TRUE, dependson=c("define_model", "settings", "sample_matrix", "source_cpp")--------------------------------------------------------------

# RUN MODEL -------------------------------------------------------------------

# Define parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Compute
Y.ti <- foreach(i=1:nrow(mat), 
                .packages = c("sensobol", "data.table", "pcaPP", 
                                 "logitnorm")) %dopar%
  {
    model_Ti(k = mat[[i, "k"]], 
             k_2 = mat[[i, "k_2"]], 
             k_3 = mat[[i, "k_3"]],
             epsilon = mat[[i, "epsilon"]],
             phi = mat[[i, "phi"]],
             delta = mat[[i, "delta"]],
             N.all = mat[[i, "N.all"]], 
             N.azzini = mat[[i, "N.azzini"]], 
             N.high = N.high)
  }

# Stop parallel cluster
stopCluster(cl)

## ----arrange_output, cache=TRUE, dependson="model_run"----------------------------------------------------------------------------------------------------------

# ARRANGE OUTPUT --------------------------------------------------------------

out_cor <- rbindlist(Y.ti, idcol = "row")

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

ggplot(indices, aes(parameters, original, fill = sensitivity)) +
  geom_bar(stat = "identity", 
           position = position_dodge(0.6), 
           color = "black") +
  geom_errorbar(aes(ymin = low.ci, 
                    ymax = high.ci), 
                position = position_dodge(0.6)) +
  facet_wrap(~estimator, 
             ncol = 1) +
  labs(x = "", 
       y = "Sobol' index") +
  scale_fill_discrete(name = "Sobol' indices", 
                      labels = c(expression(S[i]), 
                                 expression(T[i]))) +
  theme_AP() +
  theme(legend.position = "top")


indices[sensitivity == "Si", sum(original), estimator]


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












