context('testing generation of transcript to biotype file')
library(AnnotationHub)

test_that("test generation of transcript to biotype file", {
  bgee <- new("BgeeMetadata")
  kallisto <- new ("KallistoMetadata")
  user <- new("UserMetadata")
  user@species_id <- "6239" # C. elegans
  ah = AnnotationHub()
  # query the annotation hub
  potential_datasets <- query(ah, c("GTF","Ensembl", "Caenorhabditis elegans", "Caenorhabditis_elegans.WBcel235.84"))
  # retrieve dataset locally and keep path to local file
  user <- setAnnotationFromObject(user, potential_datasets[["AH50789"]], "testAnnotation")
  expect_failure(expect_error(test <- load_transcript_to_biotype(kallisto, bgee, user), NULL))
  expect_equal(nrow(test), 62111)
  expect_equal(ncol(test), 3)
  unlink(BgeeCall:::get_bgee_release_path(myBgeeMetadata = bgee, myUserMetadata = user), recursive = TRUE)
})