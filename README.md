# PA-h2
Proportional Amplification Heritability: Partition variance of interaction-QTL effects between proportional-amplification-mediated effects and non-proportional-amplification mediated effects

## Installing dependencies for PA-h2

To install dependencies: 
```
git clone https://github.com/BennyStrobes/PA-h2
cd PA-h2
conda env create -f environment.yml
conda activate PA-h2
```

## Running PA-h2

PA-h2 takes as input the same file-type as used by [TensorQTL](https://github.com/broadinstitute/tensorqtl).

To run PA-h2:

```
python ${pa_h2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--output-stem $output_stem
```
where:
- ${pa_h2_code_dir} is the absolute path (including a "/" at the end) to the downloaded PA-h2 package.
- ${expression_bed_file} is the path to an expression bed file (same format as TensorQTL [example](https://github.com/broadinstitute/tensorqtl/blob/master/example/data/GEUVADIS.445_samples.expression.bed.gz))
- ${E_var_file} is the absolute path to a covariate file with a single binary covariate (encoded {0, 1}; in same format as TensorQTL [example](https://github.com/broadinstitute/tensorqtl/blob/master/example/data/GEUVADIS.445_samples.covariates.txt), but with only a single covariate).
- ${plink2-per-chrom-stem} is the absolute path to the stem of a plink2 file. For example: ${plink2-per-chrom-stem}"22.pgen" is the corresponding pgen file for chromosome 22.
- ${output_stem} is the absolute path to output files. PA-h2 makes two output files: ${output_stem}"_PA_H2_summary.txt" and ${output_stem}"_standard_interaction_h2_summary.txt"


Some notes:
- Sample names in ${expression_bed_file} and ${E_var_file} need to be identical (in both name and order)
- Currently only works with plink2 files already seperated by chromosome. 
- plink2 samples (iid column of psam) do not need to be in same order as samples in ${expression_bed_file}. There can even more samples in the plink2 file than there are in the ${expression_bed_file}. Though all samples in ${expression_bed_file} need to be found in plink2.


## To DO

- allow IID column in psam file not to be exclusively column 1 (base 0)
- allow gzipped bed file
- allow covariates to csv or tsv

