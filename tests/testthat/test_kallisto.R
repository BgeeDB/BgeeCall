context('testing kallisto')

test_that("test download of kallisto", {
  kallisto <- new("KallistoMetadata")
  user <- new("UserMetadata")
  testthat::expect_message(download_kallisto(kallisto, user), "Kallisto has been succesfully installed.")
  unlink(BgeeCall:::get_kallisto_dir_path(kallisto, user), recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})
