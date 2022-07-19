<!-- badges -->
[![R-CMD-check](https://github.com/USEPA/spmodel/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/USEPA/spmodel/actions/workflows/check-standard.yaml)

# spmodel: Spatial statistical models in R

spmodel is an R package used to fit, summarize, and predict for a variety spatial
statistical models. Parameters are estimated using various methods. Additional
modeling features include anisotropy, random effects, partition factors,
big data approaches, and more. Model-fit statistics are used to summarize, visualize,
and compare models. Predictions at unobserved locations are readily obtainable.

# Installation

Install and load the most recent approved version from CRAN by running
```r
# install the most recent approved version from CRAN
install.packages("spmodel")
# load the most recent approved version from CRAN
library(spmodel)
```

Install and load the most recent development version of`spmodel` from GitHub by running
```r
# Installing from GitHub requires you first install the remotes package
install.packages("remotes")

# install the most recent development version from GitHub
remotes::install_github("USEPA/spmodel", ref = "main")
# load the most recent development version from GitHub
library(spmodel)
```

Install the most recent development version of `spmodel` from GitHub with package vignettes by running
```r
install the most recent development version from GitHub with package vignettes
devtools::install_github("USEPA/spmodel", ref = "main", build_vignettes=TRUE)
```

View the vignettes in RStudio by running
```r
vignette("basics", "spmodel") # an overview of basic features
vignette("guide", "spmodel") # a detailed guide with to spmodel
vignette("technical", "spmodel") # technical details regarding many functions
```

To view the vignettes in a web format, visit [here](https://cran.r-project.org/package=spmodel).

Further detail regarding spmodel is contained in the package's documentation manual available for download [here](https://cran.r-project.org/package=spmodel).

## Citation

If you use spmodel in a formal publication or report, please cite it. Citing spmodel lets us devote more resources to it in the future. View the spmodel citation by running
```{r}
citation(package = "spmodel")
```

```
#> 
#> To cite spmodel in publications use:
#> 
#>   Dumelle, Michael., Higham, Matt., and Ver Hoef, Jay M. (2022).
#>   spmodel: Spatial Statistical Modeling and Prediction. R package version 0.1.0.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {spmodel: Spatial Statistical Modeling and Prediction},
#>     author = {Michael Dumelle and Matt Higham and Jay M. {Ver Hoef}},
#>     year = {2022},
#>     note = {R package version 0.1.0},
#>   }
```

## Code Coverage

The [covr](https://cran.r-project.org/package=covr) package measures the percentage of code being exercised by a set of tests as an indirect measure of test quality and completeness. spmodel's code coverage is 95.66%.

## Package Contributions

We encourage users submit GitHub issues and enhancement requests so we may
continue to improve spmodel.

## EPA Disclaimer

The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.

### License

This project is licensed under the GNU General Public License, [GPL-3](https://cran.r-project.org/web/licenses/GPL-3).  
