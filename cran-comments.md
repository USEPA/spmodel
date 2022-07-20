This is a new release

-------

## New Release

This is a new release.

## Test environments

Several tests were run using rhub. Below are the results for windows, linux, and mac os builds.

* `rhub::check_on_windows()`
    * Platform: Windows Server 2022, R-release, 32/64 bit
        * Status: SUCCESS (no ERRORs, WARNINGs, or NOTEs)
        * Build ID: https://builder.r-hub.io/status/spmodel_0.1.0.tar.gz-7c88621231cf41ffa4f661dd79353b34
        
* `rhub::check_on_linux()`
    * Platform: Ubuntu Linux 20.04 1 LTS, R-release, GCC
        * Status: SUCCESS (no ERRORs, WARNINGs, or NOTEs)
        * Build ID: https://builder.r-hub.io/status/spmodel_0.1.0.tar.gz-ccfe2a1de93b44b79dc478ec606d315b
        
* `rhub::check(platform = "debian-gcc-release")`
    * Platform: Debian Linux, R-release, GCC
        * Status: SUCCESS (no ERRORs, WARNINGs, or NOTEs)
        * Build ID: https://builder.r-hub.io/status/spmodel_0.1.0.tar.gz-fd988dd5be0842f3ae5f131297204c4f

* `rhub::check(platform = "macos-highsierra-release-cran")`
    * Platform: macOS 10.13.6 High Sierra, R-release, CRAN's setup
        * Status: SUCCESS (no ERRORs, WARNINGs, or NOTEs)
        * Build ID: https://builder.r-hub.io/status/spmodel_0.1.0.tar.gz-b7b851d01fb0465abb64fff84255c3b4
        
* rhub tests were available on Solaris but not tested, as CRAN does not appear to
  perform Solaris checks anymore.
  
* The CRAN check may find the following `NOTE`s
    * Possibly misspelled words in `DESCRIPTION`:
        * anisotropy
        * This word is spelled correctly.
    * Problems when formatting `CITATION` entries:
        * x:1: unexpected '}'
        * An author has a two-word last name, so the `{}` are used in the `CITATION` file to help ensure proper BibTeX formatting.
        
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
