# spmodel 0.7.0

## Minor Updates

* Added `AUROC()`  functions to compute the area under the receiver operating characteristic (AUROC) curve for `spglm()` and `spgautor()` models when `family` is `"binomial"` and the response is binary (i.e., represents a single success or failure).
* Added a `BIC()` function to compute the Bayesian Information Criterion (BIC) for `splm()`, `spautor()`, `spglm()`, and `spgautor()` models when `estmethod` is `"reml"` (restricted maximum likelihood; the default) or `"ml"` (maximum likelihood).
* Added a `type` argument to `loocv()` when `cv_predict = TRUE` and using `spglm()` or `spgautor()` models so that predictions may be obtained on the link or response scale.
* Added a warning message when `data` is an `sf` object and a geographic (i.e., degrees) coordinate system is used instead of a projected coordinate system.
* Changed the default behavior of `local` in `predict.spmodel` so that it depends only on the observed data sample size. Now, when the observed data sample size exceeds 10,000 `local` is set to `TRUE` by default. This change was made because prediction for big data depends almost exclusively on the observed data sample size, not the number of predictions desired.
* Minor external data updates (for package testing).
* Minor vignette updates.
* Minor documentation updates.
* Minor error message updates.

## Bug Fixes

* Fixed a bug that prohibited proper indexing when calling `predict()` with the `local` method `"distance"` on a model object fit with a random effect or partition factor.

# spmodel 0.6.0

## Minor Updates

* Improved efficiency of handling random effects in big data models fit using `splm(..., local)` and `spglm(..., local)`.
* Changed `Matrix::rankMatrix(X, method = "tolNorm2")` to `Matrix::rankMatrix(X, method = "qr")` when determining linear independence in `X`, the design matrix of explanatory variables.
* Replaced an error message with a warning message when `X` has perfect collinearities (i.e., is not full rank). If this warning message occurs, it is possible that a subsequent error occurs while model fitting resulting from a covariance matrix that is not positive definite (i.e., a covariance matrix that is singular or computationally singular).
* Improved efficiency of `splm()` when `spcov_type` is `"none"` and there are no random effects [(#15)](https://github.com/USEPA/spmodel/issues/15).
* Added a `range_positive` argument to `spautor()` and `spgautor()` that when `TRUE` (the new default), restricts the range parameter to be positive. When `FALSE` (the prior default), the range parameter may be negative or positive.
* Updated the initial parameter grid search for `spautor()` and `spgautor()` to include range parameter values near the lower and upper boundaries.
* Minor documentation updates

## Bug Fixes

* Fixed a bug that yielded improper predictions when performing local prediction (specifying `local` in a call to `predict(object, newdata, ...)`) when the model object (`object`) was fit using `splm(formula, ...)` or `spglm(formula, ...)` and `formula` contained at least one call to `poly(..., raw = FALSE)`.
* Fixed a bug that caused big data models fit using `splm(..., local)` and `spglm(..., local)` to fail when a user-specified local index was passed to `local` that was a factor variable and at least one factor level not was observed in the local index.
* Fixed a bug that caused models fit using `splm(..., partition_factor)` and `spglm(..., partition_factor)` to fail when the partition factor variable was a factor variable and at least one factor level was not observed in the data.
* Fixed a bug in `spgautor()` that inflated the covariance matrix of the fixed effects (accessible via `vcov()`).
* Fixed a bug in `sp*(spcov_params, ...)` simulation functions that caused an error when `spcov_params` had class `"car"` or `"sar"` and `W` was provided by the user.

# spmodel 0.5.1

## Minor Updates

* Set a default value of `newdata_size = 1` when `newdata_size` was omitted while predicting `type = "response"` for binomial families.
* Improved computational efficiency of `loocv(object)` when `object` was created using `splm()` or `spglm()`, `spcov_type` was `"none"`, and there were no random effects specified via `random`.
* Changed the number of k-means iterations from 10 to 30 (when fitting models using the `local` argument to `splm()` or `spglm()`).
* Added bias and root-mean-squared-prediction error to `loocv(object)`. When `object` was created using `splm()` or `spautor()`, `loocv(object)` added the squared correlation between the observed data and leave-one-out predictions, regarded as a prediction r-squared.
* Improved prediction efficiency (using `predict()` or `augment()`) for `splm()` objects when `spcov_type` was `"none"` and there were no random effects.
* Minor error message updates.

## Bug Fixes

* Fixed a bug that caused local prediction to fail when the fitted model used a partition factor ([#13](https://github.com/USEPA/spmodel/issues/13)).
* Fixed a bug that caused significant increases in computational and memory demands when calling `loocv(object, local, ...)` if `object` was created using `splm(..., random)` or `spglm(..., random)` (i.e., when random effects were specified via the `random` argument to `splm()` or `spglm()`).
* Fixed a bug that caused significant increases in computational and memory demands when calling `loocv(object, local, ...)` if `object` was created using `splm(..., partition_factor)` or `spglm(..., partition_factor)` (i.e., when a partition factor was specified via the `partition_factor` argument to `splm()` or `spglm()`).


# spmodel 0.5.0

## Minor updates

* Predictions can now be made for prediction locations whose random effect levels are not present in the observed data
    * When this occurs, the random-effect covariance between the observed data and these prediction locations is assumed to be zero.
* The default for `local = TRUE` in `splm()` and `spglm()` now uses the `kmeans` assignment method with group sizes approximately equal to 100.
    * Previously, the `random` assignment method was used with group sizes approximately equal to 50.
* The default for `local = TRUE` in `predict()` and `augment()` now uses 100 local neighbors.
    * Previously, 50 local neighbors were used.
* Moved the "A Detailed Guide to `spmodel`" and "Technical Details" vignettes to the package website.
* Added a "Spatial Generalized Linear Models in `spmodel`" vignette to the package website.
* Changed name of "An Overview of Basic Features in `spmodel`" vignette to "An Introduction to `spmodel`" and changed output type from PDF to HTML.
* Other minor vignette updates.
* Minor documentation updates.

## Bug fixes

* Fixed a bug that occurred with prediction for success/failure binomial data (e.g., Bernoulli data) when `local` in `predict()` was `TRUE`.
* Fixed a bug that could affect simulating data using `sprbinom()` when the `size` argument was different from `1`.
* Fixed a bug that could cause local prediction to fail when only one level of a random effect was present in the prediction site's local neighborhood.
* Fixed a bug that could cause an error when local estimation was used for the `"sv-wls"` estimation method.
* Fixed a bug that caused undesirable behavior from `tidy()` when `conf.level` was less than zero or greater than one.

# spmodel 0.4.0

## Major updates

* Added an `spglm()` function to fit spatial generalized linear models for point-referenced data (i.e., generalized geostatistical models).
    * `spglm()` syntax is very similar to `splm()` syntax.
    * Poisson, negative binomial, binomial, beta, gamma, and inverse Gaussian families are accommodated.
    * `spglm()` fitted model objects use the same generics as `splm()` fitted model objects.
* Added an `spgautor()` function to fit spatial generalized linear models for areal data (i.e., spatial generalized autoregressive models).
    * `spgautor()` syntax is very similar to `spautor()` syntax.
    * Poisson, negative binomial, binomial, beta, gamma, and inverse Gaussian families are accommodated.
    * `spgautor()` fitted model objects use the same generics as `spautor()` fitted model objects.

## Minor updates

* In `augment()`, made the `level` and `local` arguments explicit (rather than being passed to `predict()` via `...`).
* Added `offset` support for relevant modeling functions.
* Minor documentation updates.
* Minor vignette updates.

## Bug fixes

* Fixed a bug in `spcov_params()` that yielded output with improper names when a named vector was used as an argument.
* Fixed a bug in `spautor()` that did not properly coerce `M` if given as a matrix (instead of a vector).
* Fixed a bug in `esv()` that prevented coercion of `POLYGON`geometries to `POINT` geometries if `data` was an `sf` object.
* Fixed a bug in `esv()` that did not remove `NA` values from the response.
* Fixed a bug in `splm()` and `spautor()` that caused an error when random effects or partition factors were ordered factors.
* Fixed a bug in `spautor()` that prevented an error from occurring when a partition factor was not categorical or not a factor
* Fixed a bug in `covmatrix(object, newdata)` that returned a matrix with improper dimensions when `spcov_type` was `"none"`.
* Fixed a bug in `predict()` that caused an error when at least one level of a fixed effect factor was not observed within a local neighborhood (when the `local` method was `"covariance"` or `"distance")`.
* Fixed a bug in `cooks.distance()` that used the Pearson residuals instead of the standarized residuals.

# spmodel 0.3.0

## Minor updates

* Added the `varcomp` function to compare variance components.
* Added an error message when there are `NA` values in predictors.
* Added an error message when the design (model) matrix is not invertible (i.e., perfect collinearities are detected).
* Added support for plotting anisotropic level curves of equal correlation when the `which` argument to `plot()` contains `8`.
* Renamed `residuals()` type `raw` to `response` to match `stats::lm()`.
* Changed class of `splm()` output to `"splm"` from `"spmod"` or `"splm_list"` from `"spmod_list"`.
* Changed class of `spautor()` output to `"spautor"` from `"spmod"` or `"spautor_list"` from `"spautor_list"`.
* Changed class of `splmRF()` output to `"splmRF"` from `"spmodRF"` or `"splmRF_list"` from `"spmodRF_list"`.
* Changed class of `spautorRF()` output to `"spautorRF"` from `"spmodRF"` or `"spautorRF_list"` from `"spmodRF_list"`.
* Methods corresponding to a generic function defined outside of `spmodel` are now all documented using
    an `.spmodel` suffix, making it easier to find documentation of a particular
    `spmodel` method for the generic function of interest.
* Added an error when random effect grouping variables or partition factors are numeric.
* Added an error when random effect or partition factor levels in `newdata` are not also in `data`.
* Updated citation information.

## Bug fixes

* Fixed a bug that produced irregular spacing in an error message for `spcov_initial()`.
* Fixed a bug that prevented proper display of row names when calling `predict()`
    with `interval = "confidence"`.
* Fixed a bug that sometimes caused miscalculations in model-fitting and prediction
    when random effect or partition factor variables were improperly coerced to a different type.
* Fixed bugs that sometimes caused miscalculations in certain model diagnostics.
* Fixed inconsistencies in several non-exported generic functions.
* Fixed a bug that prevented names from appearing with output from certain model diagnostics.

# spmodel 0.2.0

* `spmodel` v0.3.0 changed the names of `spmod`, `spmodRF`, `spmod_list`, and `spmodRF_list` objects.

## Minor updates

* `splm()` and `spautor()` allow multiple models to be fit when the `spcov_type` argument is a vector of length greater than one or the `spcov_initial` argument is a list (with length greater than one) of `spcov_initial` objects.
    * The resulting object is a list with class `spmod_list`. Each element of the list holds a different model fit.
    * `glances()` is used on an `spmod_list` object to glance at each model fit.
    * `predict()` is used on an `spmod_list` object to predict at the locations in `newdata` for each model fit.
* Added the `splmRF()` and `spautorRF()` functions to fit random forest spatial residual models.
    * The resulting object has class `spmodRF` (one spatial covariance) or `spmodRF_list` (multiple spatial covariances)
    * These objects are built for use with `predict()` to perform prediction.
* Added the `covmatrix()` function to extract covariance matrices from an `spmod` object fit using `splm()` or `spautor()`.
* Minor vignette updates.
* Minor documentation updates.

## Bug fixes

* Fixed a bug that prevents display of spatial covariance type in summary of `spmod` objects.
* Fixed a bug that prevented prediction of factor variables when all levels of all factor variables did not appear in `newdata`.

# spmodel 0.1.1

## Minor updates

* Updated unit tests so that they are compatible with an upcoming version of `Matrix`.

# spmodel 0.1.0

This is the initial release of spmodel.
