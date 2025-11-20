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
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--binary-E-interaction-covariate-file $E_var_file \
	--output-stem $output_stem
```

where ${pa_h2_code_dir}
