context('testing preparation of GTF file')

test_that ("test merging of transcriptome and intergenic regions", {
  bgee <- new("BgeeMetadata", intergenic_release = "0.1")
  kallisto <- new ("KallistoMetadata")
  user <- new("UserMetadata")
  user@species_id <- "6239" # C. elegans
  ah = AnnotationHub()
  # query the annotation hub
  potential_datasets <- query(ah, c("FaFile","Ensembl", "Caenorhabditis elegans", "Caenorhabditis_elegans.WBcel235.cdna.all.fa"))
  # retrieve dataset locally and keep path to local file
  user <- setTranscriptomeFromFile(user, potential_datasets[["AH49057"]]$path)
  user@working_path <- getwd()
  merge_transcriptome_and_intergenic(kallisto, bgee, user)
  # test that no error is produced
  expect_failure(expect_error(merge_transcriptome_and_intergenic(kallisto, bgee, user), NULL))
  unlink(BgeeCall:::get_intergenic_release_path(myBgeeMetadata = bgee, myUserMetadata = user), recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})