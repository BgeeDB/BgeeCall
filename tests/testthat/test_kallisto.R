context('testing kallisto')

test_that("test download of kallisto", {
  kallisto <- new("KallistoMetadata")
  user <- new("UserMetadata")
  expect_failure(expect_error(download_kallisto(kallisto, user), NULL))
  unlink(BgeeCall:::get_kallisto_dir_path(kallisto, user), recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})
