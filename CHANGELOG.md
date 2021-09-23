# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Next

### Added
- Support for spike-in based normalization
- Ability to specify summary function for all region-based normalizations

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
