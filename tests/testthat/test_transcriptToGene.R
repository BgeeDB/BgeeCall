context('testing generation of gene to biotype file')
library(AnnotationHub)

test_that("test generation of intergenic transcript to gene", {
  bgee <- new("BgeeMetadata", intergenic_release = "0.1")
  user <- new("UserMetadata")
  user@species_id <- "6239" # C. elegans
  expect_failure(expect_error(intergenic_txToGene <- intergenic_tx2gene(myBgeeMetadata = bgee, myUserMetadata = user), NULL))
  expect_equal(nrow(intergenic_txToGene), 4277)
  expect_equal(ncol(intergenic_txToGene), 2)
  unlink(BgeeCall:::get_intergenic_release_path(myBgeeMetadata = bgee, myUserMetadata = user), recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})

test_that("test creation of transcript to gene file", {
  bgee <- new("BgeeMetadata", intergenic_release = "0.1")
  kallisto <- new("KallistoMetadata")
  user <- new("UserMetadata")
  user@species_id <- "6239" # C. elegans
  user@working_path <- getwd()
  ah <- AnnotationHub()
  # query the annotation hub
  potential_datasets <- query(ah, c("GTF","Ensembl", "Caenorhabditis elegans", "Caenorhabditis_elegans.WBcel235.84"))
  # retrieve dataset locally and keep path to local file
  user <- setAnnotationFromObject(user, potential_datasets[["AH50789"]], "testAnnotation")
  expect_failure(expect_error(txToGene <- read.csv(create_tx2gene(myAbundanceMetadata = kallisto, myBgeeMetadata = bgee, myUserMetadata = user), header = TRUE, sep = "\t"), NULL))
  expect_equal(nrow(txToGene), 62111)
  expect_equal(ncol(txToGene), 2)
  unlink(BgeeCall:::get_intergenic_release_path(myBgeeMetadata = bgee, myUserMetadata = user), recursive = TRUE)
  unlink(file.path(user@working_path, "release.tsv"))
})

