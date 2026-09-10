# Several tests call plot() on fitted models. Without an active graphics
# device, R opens a default pdf() device and leaves an Rplots.pdf file behind
# in the working directory after every test run. pdf(NULL) opens a real but
# fileless pdf device (see ?pdf, "file" argument), so plotting still works but
# nothing is written to disk.
grDevices::pdf(NULL)
