# spmodel 0.2.0

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

## Bug Fixes

* Fixed a bug that prevents display of spatial covariance type in summary of `spmod` objects.
* Fixed a bug that prevented prediction of factor variables when all levels of all factor variables did not appear in `newdata`.

# spmodel 0.1.1

## Minor updates

* Updated unit tests so that they are compatible with an upcoming version of `Matrix`.

# spmodel 0.1.0

This is the initial release of spmodel.
