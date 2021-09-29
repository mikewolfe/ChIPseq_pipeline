[![Snakemake](https://img.shields.io/badge/snakemake-≥5.24.2-brightgreen.svg)](https://snakemake.bitbucket.io) [![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)
# A fully portable and reproducible ChIP-seq pipeline

The goal of this project is to create a reproducible analysis pipeline for ChIP-seq data 
that can be easily deployed in any computational environment. Currently this pipeline is
aimed at processing ChIP-seq samples with the following 
characteristics:

- Paired-end short-read illumina sequencing data
- Smaller genome sizes such as bacterial genomes

# What the pipeline does
Starting from raw fastq files this pipeline does the following:

- preprocessing to remove adapters and low quality sequences using [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/)
  and [`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic)
- alignment to (a) reference genome(s) using [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- read coverage calculation with a variety of normalization options within or between samples 
  with [`deeptools`](https://deeptools.readthedocs.io/en/latest/) and custom python scripts 
- peak calling using either [`macs2`](https://github.com/macs3-project/MACS) or a custom python 
  implementation of [`CMARRT`](https://github.com/mikewolfe/CMARRT_python)
- variant calling against a reference genome on ChIP input samples using
  [`breseq`](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing)


In addition this pipeline uses [`multiqc`](https://multiqc.info/) to compile the 
following quality control metrics into an interactive html report:
- Read QC using [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) both before and after preprocessing
- A large number of ChIP quality control metrics calculated by [`deeptools`](https://deeptools.readthedocs.io/en/latest/)

# Installing the pipeline

This pipeline has three dependencies:
- The package manager [`miniconda`](https://docs.conda.io/en/latest/miniconda.html)
- The workflow management system [`snakemake`](https://snakemake.readthedocs.io/en/stable/index.html)
- The API for working with Portable Encaspulated Projects (PEP) [`peppy`](http://peppy.databio.org/en/latest/)

`miniconda` can be installed following the installation instructions for your system [here](https://docs.conda.io/en/latest/miniconda.html).

Once `miniconda` is installed, both `snakemake` and `peppy` can be installed in their own environment easily using:
```
conda create -n ChIPseq_pipeline snakemake=5.24.2 peppy
conda activate ChIPseq_pipeline
```

**Note** If you are using a computational cluster that requires job management
software, you may want to install that with your environment as well.
For example, if you are using an `htcondor`-managed server you would
instead create your environment like so:
```
conda create -n ChIPseq_pipeline snakemake=5.24.2 peppy htcondor=8.9.5
conda activate ChIPseq_pipeline
```

Now you can pull the pipeline from github using:
```
git clone --recurse-submodules https://github.com/mikewolfe/ChIPseq_pipeline/
```

And you can change into the newly cloned `ChIPseq_pipeline` directory and test your installation with:
```
snakemake --use-conda --cores 10
```

Or if using a cluster with job management software you can run this
with an environment-specific
[profile](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles).
For example:
```
snakemake --use-conda --cores 10 --profile htcondor
```


This will run the entire pipeline using the provided test data consisting of small example fastqs heavily downsampled from real ChIP data. 

The first time you run the pipeline it will need to create dedicated `conda` environments for each module which will take some time. Afterwards, it will run quickly. For more information on using `conda` with `snakemake` including how to set things up to run offline check out the documentation [here](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management).

If everything runs smoothly you can then clean up and remove the results from the test data using:
```
snakemake clean_all --cores 1
```

# Running the pipeline

This pipeline uses `snakemake` to manage the workflow and familiarity with `snakemake` will help with getting the most out of the pipeline. Fortunately, `snakemake` has an excellent tutorial that can be found [here](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) if you are unfamiliar with this workflow management system.

## Input

This pipeline takes as input a [Portable Encaspulated Project (PEP)](http://pep.databio.org/en/latest/) which is essentially a `.csv` of samples and metadata together with a `.yaml` file allowing for extensions to the existing metadata.

The following required fields are needed for this pipeline
- `sample_name` - a unique identifier for each sample
- `filenameR1` - the base file name for read1 of a set of paired fastq files
- `filenameR2` - the base file name for read2 of a set of paired fastq files
- `file_path` - the path to where the files for a given sample live
- `input_sample` - which unique `sample_name` acts as the input for each extracted sample.
Input samples should leave this field blank.
- `genome` - what reference genome should this sample be aligned to

An example of a sample sheet is included at [pep/test_samples.csv](pep/test_samples.csv).

Additionally, the sample sheet can be augmented with a required `config.yaml` file. In the included test example this is used to replace the `file_path` field with a specific location. This example can be found at [pep/config.yaml](pep/config.yaml).

The `pep/config.yaml` should be edited to point towards your `pep/samples.csv` file. If you create your `samples.csv` file with excel be sure to save it as `comma seperated values` not any of the other encodings for `.csv`.

## Configuration

The pipeline itself, including parameters controlling specific tools, is controlled by a `.yaml` file in [config/config.yaml](config/config.yaml). The included `config.yaml` has all possible options specified with comments describing what those options control.

## Rules

The pipeline is organized into modules each of which runs a specific task needed for ChIP-seq analysis.

- [workflow/rules/preprocessing.smk](workflow/rules/preprocessing.smk) includes rules for trimming raw reads for adapters and quality
- [workflow/rules/alignment.smk](workflow/rules/alignment.smk) includes rules for aligning samples to their reference genome
- [workflow/rules/coverage_and_norm.smk](workflow/rules/coverage_and_norm.smk) includes rules for calculating read coverage over the genome and performing within and between sample normalization
- [workflow/rules/peak_calling.smk](workflow/rules/peak_calling.smk) includes rules for calling ChIP-seq peaks
- [workflow/rules/quality_control.smk](workflow/rules/quality_control.smk) includes rules for performing summarizing quality control on the reads themselves and ChIP-seq specific quality control
- [workflow/rules/postprocessing.smk](workflow/rules/postprocessing.smk) includes rules for getting summaries of ChIP signals over specified regions or locations
- [workflow/rules/variant_calling.smk](workflow/rules/variant_calling.smk) includes rules for checking for mutations against the reference genome. Typically run on ChIP input samples. Will take awhile to run.

Each of these rules can be run individually using:
```
snakemake run_module_name --use-conda --cores 10
```

For example:
```
snakemake run_preprocessing --use-conda --cores 10
```

Additionally to remove the output of a given module run:
```
snakemake clean_module_name --use-conda --cores 1
```

For example:
```
snakemake clean_preprocessing --use-conda --cores 1
```

Many of later modules are dependent on the earlier modules and running a later module will run the required rules in an earlier module automatically.

# Issues with the pipeline

If you run into any issues with the pipeline and would like help please submit it to the Issues page. Please include your `config/config.yaml` file, your `pep/config.yaml` file, your `pep/samples.yaml` file, and the output from `snakemake` that includes your error.

# Version history

Currently at version 0.2.2

See the [Changelog](CHANGELOG.md) for version history and upcoming
features.
