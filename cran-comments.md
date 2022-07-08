This is a new release

-------

## New Release

This is a new release.

## Test environments

* GitHub actions [here](https://github.com/USEPA/spmodel/actions/workflows/check-standard.yaml)
    * Status: Passing on macOS-latest (release), windows-latest (release),
      ubuntu-latest (devel), ubuntu-latest (release), ubuntu-latest (oldrel-1)

* `rhub::check(platform = "macos-highsierra-release-cran")`
* Platform: macOS 10.13.6 High Sierra, R-release, CRAN's setup
        * Status: SUCCESS (NOTE related to manual PDF size see next section)
        * Build ID: https:
        
* rhub tests were available on Solaris but not tested, as CRAN does not appear to
  perform Solaris checks anymore
  
* Because rhub has temporariliy stopped Windows and Linux builders to sort out
  billing issues, no rhub builds were tested on Windows or Linux. See 
  [here](https://twitter.com/rhub_/status/1542039387369885698) for more
  information.

## R CMD check results

Here is the output from `devtools::check(manual = TRUE)` on
the Windows 10 x64 operating system

0 errors | 1 warnings | 0 notes

The warning: "WARNING qpdf is needed for checks on size reduction of PDFs". 
I have received this warning while running `devtools::check()`
with my other CRAN packages, but compression of some kind seems to happen on
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
