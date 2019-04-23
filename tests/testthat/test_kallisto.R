context('testing kallisto')

test_that("test kallisto URL", {
    library(httr)
    kallisto <- new("KallistoMetadata")
    user <- new("UserMetadata")
    os_version <- get_os()
    if (os_version == "linux") {
        test_url <- http_status(GET(kallisto@kallisto_linux_url))$message
        expect_equal(test_url, "Success: (200) OK")
    } else if(os_version == "osx") {
        test_url <- http_status(GET(kallisto@kallisto_osx_url))$message
        expect_equal(test_url, "Success: (200) OK")
    } else if(os_version == "windows") {
        test_url <- http_status(GET(kallisto@kallisto_windows_url))$message
        expect_equal(test_url, "Success: (200) OK")
    }
})

test_that("test download of kallisto", {
  kallisto <- new("KallistoMetadata")
  user <- new("UserMetadata")
  expect_failure(expect_error(download_kallisto(kallisto, user), NULL))
  unlink(BgeeCall:::get_kallisto_dir_path(kallisto, user), recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})
