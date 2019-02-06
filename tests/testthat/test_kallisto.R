context('testing kallisto')

test_that("test download of kallisto", {
  kallisto <- new("KallistoMetadata")
  user <- new("UserMetadata")
  kallisto@kallisto_dir <- file.path(getwd(), "kallisto") 
  testthat::expect_message(download_kallisto(kallisto, user), "Kallisto has been succesfully installed.")
  unlink(kallisto@kallisto_dir, recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})
