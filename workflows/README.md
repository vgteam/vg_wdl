## `vg_multi_map_call.wdl`

1. Raw reads (pair of FASTQ files) are split into chunks by the *splitReads* functions.
1. Each chunk is mapped using either *runVGMAP* or *runVGMPMAP*.
1. Variants are called using either:
	- Surject mode:
		1. The chunked GAMs are surjected to chunked BAMs with *runSurject*/*sortMDTagBAMFile*/*runPICARD*.
		1. The chunked BAMs are merged with *mergeAlignmentBAMChunksVGMAP* or *mergeAlignmentBAMChunksVGMPMAP*.
		1. The BAM is chunked by chromosome with *splitBAMbyPath*.
		1. Indel realignment on each chrom with *runGATKIndelRealigner*.
		1. The chunked BAMs are merged back with *mergeIndelRealignedBAMs*.
		1. Variants are called with one of:
			- *runGATKHaplotypeCaller* (run per chrom after *runGATKIndelRealigner*)
			- *runGATKHaplotypeCallerGVCF* in gVCF mode.
			- *runDragenCaller* to use the Dragen module.
			- *runDragenCallerGVCF* to use the Dragen module in gVCF mode.
	- vg call:
		1. The chunked GAMs are merged and sorted with *mergeAlignmentGAMChunks*.
		1. The GAM file is chunked by region with *chunkAlignmentsByPathNames*.
		1. Variants are called in each chunk with *runVGCaller*.
		1. Chunked VCFs are clipped using *runVCFClipper*.
	- vg call using packed coverage:
		1. The chunked GAMs are merged without sorting with *TODO*.
		1. For each chromosome, the GAMs and graphs are chunked and SV called with *TODO*.
1. Chunked VCFs are merged and compressed with *concatClippedVCFChunks* and **bgzipMergedVCF*.
1. Optional: VCFs are annotated with SnpEff using *snpEffAnnotateVCF*.

Generic functions:

- *concatClippedVCFChunks* concatenates and sort VCF files with bcftools.
- *bgzipMergedVCF* compresses and indexes a VCF file.
- *normalizeVCF* normalizes a VCF file with bcftools.

Cleanup functions removes temporary files for different steps:

- *cleanUpGoogleFilestore*
- *cleanUpUnixFilesystem*
- *cleanUpVGMapperInputsGoogle*
- *cleanUpVGMapperInputsUnix*
- *cleanUpMergeAlignmentBAMChunksGoogle*
- *cleanUpMergeAlignmentBAMChunksUnix*
- *cleanUpGoogleFilestore*
- *cleanUpUnixFilesystem*
