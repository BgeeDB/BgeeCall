context('testing kallisto')

test_that("test download of kallisto", {
  kallisto <- new("KallistoMetadata")
  kallisto@kallisto_dir <- file.path(getwd(), "kallisto") 
  testthat::expect_message(download_kallisto(kallisto), "Kallisto has been succesfully installed.\n")
  unlink(kallisto@kallisto_dir, recursive = TRUE)
})
