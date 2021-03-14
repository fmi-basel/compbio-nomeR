test_that("get_cooccurrence_ctable_from_bams works for QuasR bam file", {
	#library(GenomicRanges)
	
	
	quasR_coocc_matrix <- get_cooccurrence_ctable_from_bams(bamfiles = "QuasR_test.bam",
																													samplenames = "quasr_pe",
																													regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
																																													 strand = "+",
																																													 IRanges::IRanges(start = 200,end = 500)),
																													genome = "random_genome_700bp.fa",
																													max_spacing = 300,
																													remove_nonunique = F,
																													clip_until_nbg = 0L,
																													max_bisC_meth = 1)
	
	quasR_expect_coocc_mat <- readRDS("QuasR_expected_coocc_mat.rds") 
	
	expect_true(all(quasR_coocc_matrix$CoocCountTable[[1]]==quasR_expect_coocc_mat))
	
})

test_that("get_cooccurrence_ctable_from_bams works for Bismark bam file", {
	
  #library(GenomicRanges)
	
	bismark_coocc_matrix <- get_cooccurrence_ctable_from_bams(bamfiles = "Bismark_test.bam",
																														samplenames = "bismark_pe",
																														regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
																																														 strand = "+",
																																														 IRanges::IRanges(start = 200,end = 500)),
																														genome = "random_genome_700bp.fa",
																														max_spacing = 300,
																														remove_nonunique = F,
																														clip_until_nbg = 0L,
																														max_bisC_meth = 1)
	
	expect_coocc_mat <- readRDS("QuasR_expected_coocc_mat.rds") 
	
	expect_true(all(bismark_coocc_matrix$CoocCountTable[[1]] == expect_coocc_mat*2))
	
})