test_that("get_protect_stats_table_from_bams works for QuasR bam file", {
  #library(GenomicRanges)
	quasr_protect_stats <- get_protect_stats_table_from_bams(bamfiles = "QuasR_test.bam",
																													 samplenames = "quasr_pe",
																													 regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
																													 																 strand = "+",
																													 																 IRanges::IRanges(start = 200,end = 500)),
																													 genome = "random_genome_700bp.fa",
																													 remove_nonunique = F,
																													 clip_until_nbg = 0L,
																													 max_bisC_meth = 1)
	## order as currently the order is random due to unordered_map
	quasr_protect_stats$ProtectStats[[1]] <- quasr_protect_stats$ProtectStats[[1]][order(quasr_protect_stats$ProtectStats[[1]][,1],
																																											 quasr_protect_stats$ProtectStats[[1]][,2]),]
	
	quasr_expect_protstats <- readRDS("QuasR_expected_protectstats.rds")
	
	expect_true(all(quasr_protect_stats$ProtectStats[[1]] == quasr_expect_protstats))
})


test_that("get_protect_stats_table_from_bams works for Bismark bam file", {
  #library(GenomicRanges)
	bismark_protect_stats <- get_protect_stats_table_from_bams(bamfiles = "Bismark_test.bam",
																														 samplenames = "bismark_pe",
																														 regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
																														 																 strand = "+",
																														 																 IRanges::IRanges(start = 200,end = 500)),
																														 genome = "random_genome_700bp.fa",
																														 remove_nonunique = F,
																														 clip_until_nbg = 0L,
																														 max_bisC_meth = 1)
	## order as currently the order is random due to unordered_map
	bismark_protect_stats$ProtectStats[[1]] <- bismark_protect_stats$ProtectStats[[1]][order(bismark_protect_stats$ProtectStats[[1]][,1],
																																													 bismark_protect_stats$ProtectStats[[1]][,2]),]
	
	expect_protstats <- readRDS("QuasR_expected_protectstats.rds")
	expect_protstats[,3] <- 2*expect_protstats[,3]
	expect_true(all(bismark_protect_stats$ProtectStats[[1]] == expect_protstats))
})