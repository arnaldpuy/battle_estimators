
# The Battle of Total-Order Sensitivity Estimators

<!-- badges: start -->
<!-- badges: end -->

This is the preliminary code of an investigation on total-order sensitivity indices. We assess which estimator (Sobol', Jansen, Homma and Saltelli, Azzini and Rosati, Janon/Monod, Owen, Glen and Isaac or Razavi and Gupta) performs better in an 8-dimension hypercube where the following parameters are treated as uncertain:

* The sampling method to create the sampliing matrix.
* The total number of model runs available.
* The function dimensionality.
* The proportion of active two-order effects.
* The proportion of active three-order effects.
* The randomness in the test function.
* The distribution of the model inputs.
* The performance measure used.

