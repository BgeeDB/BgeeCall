context('testing preparation of GTF file')

test_that ("test merging of transcriptome and intergenic regions", {
  bgee <- new("BgeeMetadata", intergenic_release = "0.1")
  kallisto <- new ("KallistoMetadata")
  user <- new("UserMetadata")
  user@species_id <- "6239" # C. elegans
  ah = AnnotationHub()
  # query the annotation hub
  ah_resources <- AnnotationHub::query(ah, c("Ensembl", "Caenorhabditis elegans", "84"))
  transcriptome_object <- rtracklayer::import.2bit(ah_resources[["AH50453"]])
  # retrieve dataset locally and keep path to local file
  user <- setTranscriptomeFromObject(user, transcriptome_object, "WBcel235")
  user@working_path <- getwd()
  merge_transcriptome_and_intergenic(kallisto, bgee, user)
  # test that no error is produced
  expect_failure(expect_error(merge_transcriptome_and_intergenic(kallisto, bgee, user), NULL))
  unlink(BgeeCall:::get_intergenic_release_path(myBgeeMetadata = bgee, myUserMetadata = user), recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})
