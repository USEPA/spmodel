This is a minor update that adds a few features and bug fixes. Thank you.

-------

## Resubmission

This is a resubmission.

## R CMD check results

Here is the output from `devtools::check(manual = TRUE)` on
the Windows 10 x64 operating system

0 errors | 1 warnings | 0 notes

The warning: "WARNING qpdf is needed for checks on size reduction of PDFs". 
I have received this warning while running `devtools::check()`
with my other CRAN packages (spsurvey and sptotal), but compression of some kind seems to happen on
CRAN's end, as these packages have made it through. I will note that locally 
installing and inspecting the PDFs, their actual file sizes are significantly
smaller than what `devtools::check()` says, so compression of some kind
seems to happen on local installation. 
In short, it is my understanding that this warning in `devtools::check()`
does not accurately reflect the size of the final PDFs installed upon package build.
I will also note that the GitHub actions checks return no warnings when passing
`build_args = "--compact-vignettes=both" to R CMD check.

## Downstream dependencies

There are no downstream dependencies.
