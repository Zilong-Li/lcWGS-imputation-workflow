# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample and reference panel sheet

* Add samples to `config/samples.tsv`. For each sample, the columns `sampleid`, `bam` and `depth` have to be defined. `sampleid` has to be in the header SM tag of corresponding `bam`. `depth` represents the sequencing coverage of the `bam`. If `downsample` is defined in `config/cofing.yaml`, the `depth` is used to calculate the fraction of downsampling. 
* Add reference panel information to `config/refpanel.tsv`.For each chromosome, the columns `chr`, `vcf`, `start` and `end` have to be defined. The `start` and `end` is for the variable SNPs in the reference VCF not for the reference genome, which is used to split the imputation tasks into multiple different chunks given `chunksize` in the `config/config.yaml`. Optional column `truth` represents the VCF file with truth genotypes of target samples for evaluating the accuracy.

# Benchmarking on sequencing coverage

Modify the list of `downsample` in `config/config.yaml`.

# Benchmarking on reference panel size

Modify the list of `refsize` in `config/config.yaml`. In default `refsize=[0]`, all samples in the reference panel excluding the target samples in `config/samples.tsv` are used.
