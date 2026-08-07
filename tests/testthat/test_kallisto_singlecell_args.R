context('testing single-cell command line construction')

# bustools available flags
bustools_count_flags <- c(
    "-o", "--output", "-g", "--genemap", "-e", "--ecmap", "-t", "--txnames",
    "-s", "--split", "-m", "--multimapping", "--genecounts", "--umi-gene",
    "--cm")

test_that("bus_output_paths names bustools files in the right format", {
    paths <- BgeeCall:::bus_output_paths("/out/dge", "ERR2854358")

    expect_equal(paths$mtx, "/out/dge/ERR2854358.mtx")
    expect_equal(paths$genes, "/out/dge/ERR2854358.genes.txt")
    expect_equal(paths$barcodes, "/out/dge/ERR2854358.barcodes.txt")
    expect_equal(paths$prefix, "/out/dge/ERR2854358")
    expect_false(grepl("/[.]mtx$", paths$mtx))
    expect_false(any(grepl("//", unlist(paths))))
})

test_that("check_run_id rejects identifiers that would make bustools crash", {
    # An empty run_id makes file.path() return character(0), BgeeCall will try
    # to run bustools with an empty -o value 
    expect_error(BgeeCall:::check_run_id(character(0)))
    expect_error(BgeeCall:::check_run_id(""))
    expect_error(BgeeCall:::check_run_id(NA_character_))
    expect_error(BgeeCall:::check_run_id(c("a", "b")))
    expect_error(BgeeCall:::check_run_id(1))
    # a path separator would redirect the output into a subdirectory
    expect_error(BgeeCall:::check_run_id("a/b"))
    expect_error(BgeeCall:::check_run_id("a\\b"))
    expect_error(BgeeCall:::check_run_id("ERR2854358"), NA)
})

test_that("Order_fastq_r1_r2 interleaves R1 and R2 files", {
    # kallisto bus needs R1_1 R2_1 R1_2 R2_2, not R1_1 R1_2 R2_1 R2_2. Passing
    # them grouped pairs one lane's barcode reads with another lane's
    # biological reads, and kallisto reports no error.
    expect_equal(
        BgeeCall:::Order_fastq_r1_r2(c("a1", "b1"), c("a2", "b2"),
            check_exists = FALSE),
        c("a1", "a2", "b1", "b2"))
    expect_equal(
        BgeeCall:::Order_fastq_r1_r2("x1", "x2", check_exists = FALSE),
        c("x1", "x2"))
})

test_that("Order_fastq_r1_r2 validates its inputs", {
    expect_error(BgeeCall:::Order_fastq_r1_r2(character(0), character(0),
        check_exists = FALSE))
    expect_error(BgeeCall:::Order_fastq_r1_r2(c("a1", "b1"), "a2",
        check_exists = FALSE))
    expect_error(BgeeCall:::Order_fastq_r1_r2("a1", character(0),
        check_exists = FALSE))

    # test on missing files
    tmp_dir <- file.path(tempdir(), "bgeecall_sc_args")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    present <- file.path(tmp_dir, c("lane1_R1.fq", "lane1_R2.fq"))
    file.create(present)
    absent <- file.path(tmp_dir, c("lane2_R1.fq", "lane2_R2.fq"))
    expect_error(
        BgeeCall:::Order_fastq_r1_r2(c(present[1], absent[1]),
            c(present[2], absent[2])),
        "lane2_R1.fq")
    expect_error(BgeeCall:::Order_fastq_r1_r2(present[1], present[2]), NA)
    unlink(tmp_dir, recursive = TRUE)
})

test_that("bustools_count_args emits only flags that bustools accepts", {
    args <- BgeeCall:::bustools_count_args(t2g_path = "tx2gene_sc.tsv",
        out_prefix = "/out/dge/RUN1", ecmap_path = "/out/dge/matrix.ec",
        txnames_path = "/out/dge/transcripts.txt",
        input_bus = "/out/dge/output_sorted.bus")

    expect_false("-em" %in% args)
    expect_true(all(grep("^-", args, value = TRUE) %in% bustools_count_flags))
    
    expect_equal(args[which(args == "-e") + 1], "/out/dge/matrix.ec")
    expect_equal(args[which(args == "-g") + 1], "tx2gene_sc.tsv")
    expect_equal(args[which(args == "-t") + 1], "/out/dge/transcripts.txt")
    expect_equal(args[which(args == "-o") + 1], "/out/dge/RUN1")
    expect_equal(args[1], "count")
    # the sorted bus file is positional and must come last
    expect_equal(args[length(args)], "/out/dge/output_sorted.bus")
})

test_that("bustools_count_args toggles the optional flags", {
    default_args <- BgeeCall:::bustools_count_args("t", "p", "e", "x", "b")
    
    expect_true("--genecounts" %in% default_args)
    expect_false("-m" %in% default_args)
    expect_false("--umi-gene" %in% default_args)
    expect_false("--cm" %in% default_args)

    expect_false("--genecounts" %in%
        BgeeCall:::bustools_count_args("t", "p", "e", "x", "b",
            genecounts = FALSE))
    expect_true("-m" %in%
        BgeeCall:::bustools_count_args("t", "p", "e", "x", "b",
            multimapping = TRUE))
    expect_true("--umi-gene" %in%
        BgeeCall:::bustools_count_args("t", "p", "e", "x", "b",
            umi_gene = TRUE))
    expect_true("--cm" %in%
        BgeeCall:::bustools_count_args("t", "p", "e", "x", "b",
            count_multiplicities = TRUE))
})

test_that("kallisto_bus_args builds a valid kallisto bus command", {
    args <- BgeeCall:::kallisto_bus_args(index_path = "transcriptome.idx",
        out_dir = "/out/dge", technology = "10xv3",
        fastq_files = c("r1.fq", "r2.fq"), threads = 2)

    expect_equal(args[1], "bus")
    expect_equal(args[which(args == "-i") + 1], "transcriptome.idx")
    expect_equal(args[which(args == "-o") + 1], "/out/dge")
    expect_equal(args[which(args == "-x") + 1], "10xv3")
    # kallisto expects its FASTQ files last, after every option
    expect_equal(args[c(length(args) - 1, length(args))], c("r1.fq", "r2.fq"))
    # numeric thread counts must reach system2() as characters
    expect_equal(args[which(args == "-t") + 1], "2")
    expect_true(is.character(args))
    expect_false(any(grepl("//", args)))

    expect_error(BgeeCall:::kallisto_bus_args("i", "o", "", "a"))
    expect_error(BgeeCall:::kallisto_bus_args("i", "o", character(0), "a"))
    expect_error(BgeeCall:::kallisto_bus_args("i", "o", NA_character_, "a"))
})

test_that("bustools_sort_args and bustools_correct_args build valid commands", {
    expect_equal(
        BgeeCall:::bustools_sort_args("sorted.bus", "output.bus", threads = 4),
        c("sort", "-t", "4", "-o", "sorted.bus", "output.bus"))
    expect_equal(
        BgeeCall:::bustools_correct_args("whitelist.txt", "corrected.bus",
            "output.bus"),
        c("correct", "-w", "whitelist.txt", "-o", "corrected.bus",
            "output.bus"))
})

test_that("run_binary reports failures", {
    skip_on_os("windows")
    expect_error(BgeeCall:::run_binary("false", character(0), "bustools count"),
        "bustools count")
    expect_error(BgeeCall:::run_binary("false", character(0), "bustools count"),
        "exit status")
    expect_equal(
        as.character(BgeeCall:::run_binary("echo", "hello", "test step")),
        "hello")
})
