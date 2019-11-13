# ampatchlab/nf-rnasnv

[![Build Status](https://codebuild.ap-southeast-2.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoicDJIMmhoK3JibzFVVUxwOXBDcXlISk5qV0tGRmYrZVZzUlVRS1VNZlh2VW9tWlpBeis5Q2VJb2piUy85WngrbHZYK0JzaVZySVExbkVuQ25WOFoyNGJ3PSIsIml2UGFyYW1ldGVyU3BlYyI6InFWWWxQVTBwQ2RYeUpBazMiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=master)](https://ap-southeast-2.console.aws.amazon.com/codesuite/codebuild/projects/nf-rnasnv/history)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.07.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

RNA-Seq variant calling nextflow pipeline

## Usage

```
Usage:
    nextflow run -profile <profile> -revision <revision> ampatchlab/nf-rnasnv [options]


Nextflow execution options:

    -profile STR
        Nextflow configuration profile to use. Available profiles include:
        'awsbatch', 'conda', 'docker' and 'singularity'

    -revision STR
        Git branch/tag (version) of this workflow to use

    -work-dir DIR
        Directory where intermediate result files are stored

    -help
        Show additional execution options and exit


Input options:

    --csv FILE
        Comma-separated list of sample and readgroup inputs

    --paired_end
        Expect entries for 'fastq1' and 'fastq2' in the input CSV


Reference genome options:

    --genome STR
        Reference genome name [Either: GRCh38, GRCm38; Default: GRCh38]

    --fasta FILE
        Override the reference genome FASTA with FILE [Default: null]

    --gtf FILE
        Override the reference genome GTF with FILE [Default: null]


Sequencing adapter options:

    --adapters STR
        The adapters to trim [Either: TruSeq, BGISeq, NONE; Default: TruSeq]

    --r1_adapter STR
        Override the sequence of the R1 adapter with STR [Default: null]

    --r2_adapter STR
        Override the sequence of the R2 adapter with STR [Default: null]


Output options:

    --refdir DIR
        Path where the reference index files will be saved [Default: ./reference]

    --outdir DIR
        Path where the results will be saved [Default: ./results]


Standard options:

    --advanced
        Show advanced usage and exit

    --help
        Show this message and exit

    --version
        Show the pipeline version and exit
```

## Advanced options

```
STAR genome generate options:

    --star_genome_chr_bin_n_bits INT
        Size of the bins for genome storage [Default: 18]

    --star_genome_sa_index_n_bases INT
        Length (bases) of the SA pre-indexing string [Default: 14]

    --star_sjdb_overhang INT
        Length of the donor/acceptor sequence on each side of the junctions [Default: 100]


Cutadapt options:

    --cutadapt_base_qual_cutoff [INT,]INT
        Trim low-quality bases from each read [Default: 20]

    --cutadapt_min_read_length INT[:INT]
        Discard reads shorter than INT [Default: 20]


BCFtools mpileup options:

    --mpileup_num_regions INT
        Number of regions per mpileup process [Default: 10000]

    --mpileup_max_depth INT
        Maximum number of reads to pileup at a given position [Default: 20000]

    --mpileup_min_bq INT
        Minimum base quality for a base to be considered [Default: 13]


Ensembl VEP options:

    --vep_cache_type STR
        Alternate cache to use [Either: refseq, merged; Default: null]

    --vep_indexed_cache_file FILE
        Override the Ensembl VEP indexed cache file with FILE [Default: null]

    --vep_species_name STR
        Override the Ensembl VEP species name with STR [Default: null]

    --vep_assembly_name STR
        Override the Ensembl VEP assembly name with STR [Default: null]


MultiQC options:

    --multiqc_config FILE
        MultiQC YAML config file [Default: [:]/assets/multiqc_config.yaml]


Report options

    --execution_report STR
        Name of the Nextflow execution report to generate [Default: ./reports/execution_report.html]

    --trace_report STR
        Name of the Nextflow trace report to generate [Default: ./reports/trace_report.txt]

    --timeline_report STR
        Name of the Nextflow timeline report to generate [Default: ./reports/timeline_report.html]

    --flowchart STR
        Name of the Nextflow flowchart to generate [Default: ./reports/flowchart.png]


AWS Batch options

    --aws_queue STR
        AWS Batch JobQueue definition [Default: false]

    --aws_region STR
        AWS Region definition [Default: false]
```

## Inputs

For paired-end data, the input CSV must have the following required columns:

 * sample: Unique sample name or ID (required)
 * readgroup: Unique readgroup name or ID (optional)
 * fastq1: Absolute path of the 'R1' FASTQ file (required)
 * fastq2: Absolute path of the 'R2' FASTQ file (required)

For single-end data, the CSV must have the following columns:

 * sample: Unique sample name or ID (required)
 * readgroup: Unique readgroup name or ID (optional)
 * fastq: Absolute path of the FASTQ file (required)

If a particular sample has multiple FASTQ files (or pairs of FASTQ files), then these may
be specified on additional lines with a unique readgroup identifier. All readgroups belonging
to a particular sample will be aligned and merged using STAR.

## Notes

### GRCh37

Input FASTA files must be sorted by karyotypic order or RNA-SeQC v1.1.8 (GATK) will complain.
Unfortunately, the current GRCh37 assembly available from GENCODE is sorted lexicographically.
This appears to be an artifact of the liftover process. To use the GRCh37 reference, please
re-sort the FASTA file using the instructions below. It is hoped that this issue is resolved
in later releases.

1. Obtain the curent GRCh37 FASTA and GTF files:
```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz
```

2. Split the GRCh37 FASTA into separate chromosomes/contigs:
```
$ zcat GRCh37.primary_assembly.genome.fa.gz | awk '/^>/ { f=substr($1,2) ".fa" } { print > f }'
```

3. Replace the existing assembly with an ordered set of chromosomes/contigs:
```
$ cat chr{{1..22},X,Y,M}.fa GL*.fa | gzip > GRCh37.primary_assembly.genome.fa.gz
```

4. Clean up:
```
$ rm *.fa
```
