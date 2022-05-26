# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.2.9 - 2022-05-26

### Added
- Add a motif calling module. Supports both `meme` and `streme` motif callers
- Added fimo calling with discovered motifs
- Add ability get coverage around the set of peaks for each sample
- Combine bed files from each chromosome into one master bed file
- Create an annotation file that maps unique name to locus tags and gene names

### Changed
- Split macs2 broad and macs2 default into their own peak callers. Always call
  summits when running macs2 without broad specified

## 0.2.8 - 2022-01-20

### Added
- Ability to collapse multiple bigwig files into one using summary stats
- Ability to perform operations between two groups of bigwig files

## 0.2.7 - 2022-01-10

### Added
- Added support for single end libraries

## 0.2.6 - 2021-12-08

### Added
- Added relative polymerase progression calculation to `bwtools_query` (http://dx.doi.org/10.1016/j.molcel.2014.02.005)
- Add the ability to smooth raw signals using convolution with a Gaussian or flat kernel. 
- Add the ability to smooth raw signals using a Savitzky-Golay filter
- Ability to calculate the region-level spearman correlations for a set of samples in `bwtools_query`
- Add calculations of the traveling ratio (10.1016/j.molcel.2008.12.021)
- Add discovery of a local max within a region

### Changed
- Output of `bwtools_query` is now compressed (gzip). File extension changed to `.tsv.gz`

## 0.2.5 - 2021-10-25

### Added
- Allow mix and match between `querysub`, `queryscale` `fixedsub`, and `fixedscale`.

## 0.2.4 - 2021-10-13

### Added
- Create a compiled html table of all mutations called across all samples in `variant_calling`

### Bug fixes
- Fixed a bug in renaming `sample/output/output.gd` files to `renamed_output/sample.gd` in `variant_calling`

## 0.2.3 - 2021-10-12

### Added
- A `variant_calling` module to call changes from the reference genome with ChIP input samples ([#2](https://github.com/mikewolfe/ChIPseq_pipeline/issues/2))
- Add option to control what fraction of a region can have NaNs in it and still report a summary value in `bwtools query`

## 0.2.2 - 2021-09-26

### Added
- Support for spike-in based normalization
- Ability to specify summary function for all region-based normalizations
- Ability to specify multiple quality fields for genbank parsing to name
  parsed regions.

### Bug fixes
- Fixed bug in error reporting for missing parameters in config file

## 0.2.1 - 2021-09-22

### Added
- Support for division of the input with no log2 transformation
- Support for subtract and scaling by a known and fixed set of regions

### Bug fixes
- Cleaned up multiqc report config file to handle automatic sample name generation. 
  Fixes issue ([#10](https://github.com/mikewolfe/ChIPseq_pipeline/issues/10))

## 0.2.0 - 2021-09-08

### Added
- Support for subtracting rather than dividing the input from the extracted
- A summarize method for `bwtools` that gives total signal per contig

### Changed
- Refactored the config file syntax to be more general for `peak_calling`
- Moved output of ext to inp comparisions from `results/coverage_and_norm/deeptools_log2ratio/` to `bwtools_compare/`

## 0.1.1 - 2021-06-03

### Bug Fixes
- Added biopython dependency to the alignment environment

## 0.1.0 - 2021-05-26

### Added
- Full support for genomes with multiple contigs/chromosomes ([#5](https://github.com/mikewolfe/ChIPseq_pipeline/issues/5))
- Full support of masking regions of the genome from analysis ([#11](https://github.com/mikewolfe/ChIPseq_pipeline/issues/11))
- Full support for normalization based on a known background region ([#6](https://github.com/mikewolfe/ChIPseq_pipeline/issues/6))
- Calculate summaries over regions or locations ([#14](https://github.com/mikewolfe/ChIPseq_pipeline/issues/14))
- Ability to specify parameter strings for certain rules either as a single
  value for every sample or through specifying a column in the metadata sheet
  for by sample control
- Ability to add pseudocount before within sample normalization. Only implemented for median normalization
- Support for normalization based on averages over a given set of regions.
  Enables Occupancy Apparent calculations.
- Support for normalization for occupancy apparent calculations by dynamically
  determining the highest and lowest N genes by average signal
- Support for reciprocal ratio calculations to keep ratios in a linear scale.
  Inverts all ratios less than one and takes negative value. Then subtracts or
  adds 1 to recenter to zero.
- Added a specification to exclude Infs and NaNs from bigwig files

### Changed
- Config syntax change for specifying genome inputs
- Removed general pseudocount specification. Never adds a pseudocount for ratios
- Substantially updated the `bwtools.py` module to enable summary calculations
- Config syntax change for coverage and normalization specification

### Bug fixes
- Fixed an issue when running `coverage_and_norm` module only ([#3](https://github.com/mikewolfe/ChIPseq_pipeline/issues/3))
- Fixed an issue with the `tbb` version being too high on some systems causing bowtie2 to fail to run (`workflow/envs/alignment.yaml` file).
- Fixed an issue where genome size was not properly read as input to Macs2 in `workflow/rules/peak_calling.smk`
- Fixed an issue where filename paths with spaces in them would not input correctly to cutadapt `workflow/rules/preprocessing.smk`
- Upgraded `numpy` to version 1.20 in peak calling environment to fix issue from incompatibility in `workflow/rules/peak_calling.smk` ([#16](https://github.com/mikewolfe/ChIPseq_pipeline/issues/16))
- Fixed an issue where raw fastqs with spaces in the path were not read properly
