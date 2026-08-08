# These files were renamed from "test-*.R" to "old-*.R" so that neither
# devtools::test() nor testthat::test_dir() picks them up automatically --
# both only ever look for files matching "^test.*\\.[Rr]$" (this is hardcoded
# inside testthat::find_test_scripts() and can't be overridden via test_dir()'s
# arguments), so a "test_dir()" call pointed at this folder will report
# "No test files found" even though the files are right here.
#
# To run every archived test in this folder, from the package root:
#
#   devtools::load_all()
#   source("tests/testthat/deprecated/run_deprecated_tests.R")
#
# (testthat::test_file() runs a single file regardless of its name, so looping
# it over the directory listing sidesteps the "test" prefix requirement above.)

devtools::load_all()
grDevices::pdf(NULL) # avoid leaving an Rplots.pdf behind from plot() calls
invisible(lapply(
  list.files("tests/testthat/deprecated", pattern = "^old-.*\\.R$", full.names = TRUE),
  testthat::test_file
))
