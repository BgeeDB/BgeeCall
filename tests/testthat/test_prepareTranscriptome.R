context('testing preparation of GTF file')

test_that ("test merging of transcriptome and intergenic regions", {
  bgee <- new("BgeeMetadata", intergenic_release = "0.1")
  kallisto <- new ("KallistoMetadata")
  user <- new("UserMetadata")
  user@species_id <- "6239" # C. elegans
  transcriptome_file <- system.file("extdata", "transcriptome.fa", package = "BgeeCall")
  # retrieve dataset locally and keep path to local file
  user <- setTranscriptomeFromFile(user, transcriptome_file, "WBcel235")
  user@working_path <- getwd()
  # test that no error is produced
  expect_failure(expect_error(merge_transcriptome_and_intergenic(kallisto, bgee, user), NULL))
  unlink(BgeeCall:::get_intergenic_release_path(myBgeeMetadata = bgee, myUserMetadata = user), recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})
