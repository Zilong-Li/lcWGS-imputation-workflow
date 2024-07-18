# General settings

To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample and reference panel sheet

* Add samples to `config/samples.tsv`. For each sample, the columns `sampleid`, `bam` and `depth` have to be defined. `sampleid` has to be in the header SM tag of corresponding `bam`. `depth` represents the sequencing coverage of the `bam`. If `downsample` is defined in `config/cofing.yaml`, the `depth` is used to calculate the fraction of downsampling. 
* Add reference panel information to `config/refpanel.tsv`. For each chromosome, the columns `chr`, `vcf`, `start` and `end` have to be defined. The `start` and `end` refer to the position of variable SNPs in the reference VCF not for the reference genome, which is used to split the imputation tasks into multiple different chunks given `chunksize` in the `config/config.yaml`. One can also use optional column `quilt_chunk`, `glimpse_chunk` to specify the chunks from QUILT and GLIMPSE. When making accuracy plots, column `truth` needs to be defined which is a VCF file with truth genotypes of target samples. Optionally, column `af` can be used to specify the allele frequency of each variable site instead of calculated from the reference panel. The `af` represents a plain tab separated file of 5 columns `chr,pos,ref,alt,af` **without header**. Also, optional column `quilt_map`,`glimpse_map` refers to the genetic map file for *QUILT* and *GLIMPSE* respectively. Lastly, optional column `exclude_samples` can be used to remove extra samples in addition to target samples from the reference panel.

* Imputation test on single region. Add `region` column in `config/refpanel.tsv` which is of `chrom:start-end` form.

# Benchmarking on sequencing coverage

Modify the list of `downsample` in `config/config.yaml`.

# Benchmarking on reference panel size

Modify the list of `refsize` in `config/config.yaml`. In default `refsize=[0]`, all samples in the reference panel excluding the target samples in `config/samples.tsv` are used.

# Different scenario run

Modify the variable `scenario` in `config/config.yaml` to decide what analyses to run.

- `all` : run all comparisions across all programs and settings. requires `truth` to be configured.
- `speed`: run only speed comparisons across all programs and settings. 
- `quilt1`: run only `QUILT`.
- `quilt2`: run only `QUILT2`, which is the default.
- `glimpse1`: run only `GLIMPSE`.
- `glimpse2`: run only `GLIMPSE2`.
