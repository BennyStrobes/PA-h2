import numpy as np
import os
import sys
import pdb
import argparse
import pgenlib as pg
from tqdm import tqdm
import statsmodels.api as sm
import pickle
import random

def load_in_genotype(plink2_genotype_stem, chrom_num):
	pgen_path = plink2_genotype_stem + str(chrom_num) + '.pgen'
	pvar_path = plink2_genotype_stem + str(chrom_num) + '.pvar'
	psam_path = plink2_genotype_stem + str(chrom_num) + '.psam'

	# Load in thing
	pgen = pg.PgenReader(pgen_path.encode())

	n_samples  = pgen.get_raw_sample_ct()
	n_variants = pgen.get_variant_ct()

	pvar = np.loadtxt(pvar_path, dtype=str, delimiter='\t')
	psam = np.loadtxt(psam_path, dtype=str, delimiter='\t')

	if pvar.shape[0] != n_variants:
		print('assumption erroror')
		pdb.set_trace()
	if psam.shape[0] != n_samples:
		print('assumpton error')
		pdb.set_trace()

	return pgen, pvar, psam

def load_in_E_variable(covariate_file):
	raw_cov = np.loadtxt(covariate_file, dtype=str,delimiter='\t')
	sample_names = raw_cov[0, 1:]
	interaction_var = raw_cov[1,1:].astype(float)

	if np.array_equal(np.sort(np.unique(interaction_var)),np.asarray([0.0, 1.0])) == False:
		print('Fatal error: Interaction variable must be binary and encoded as 0 / 1')
		sys.exit(1)

	return interaction_var, sample_names

def reorder_genotype_samples(psam, E_var_sample_names):
	geno_samples = psam[:, 1]
	mapping = {}
	for ii, geno_sample in enumerate(geno_samples):
		mapping[geno_sample] = ii

	reordering_vec = []
	for sample_name in E_var_sample_names:
		if sample_name not in mapping:
			print('Fatal error: Sample in interaction (E) variable that is not found in plink2 IID')
			sys.exit(1)
		reordering_vec.append(mapping[sample_name])
	reordering_vec = np.asarray(reordering_vec)

	# Quick error check
	if np.array_equal(geno_samples[reordering_vec], E_var_sample_names) == False:
		print('Fatal error: something went wrong with reordering genotype samples')
		sys.exit(1)

	return reordering_vec


def extract_dosages_from_reader(ra, cis_snp_indices, dtype=np.float32):
	"""
	Extract dosage genotypes from an existing pgenlib.PgenReader.

	Parameters
	----------
	ra : pgenlib.PgenReader
		Already opened reader.
	cis_snp_indices : np.ndarray[bool]
		Boolean mask of length K (#variants).
	dtype : numpy dtype
		Output dtype (default float32 for dosage values).

	Returns
	-------
	D : np.ndarray of shape (M, N)
		Dosage matrix for selected variants.
	idx : np.ndarray[int] of length M
		The variant indices selected (to align with .pvar).
	"""
	n_samples  = ra.get_raw_sample_ct()
	n_variants = ra.get_variant_ct()

	# Convert mask → indices
	idx = np.flatnonzero(cis_snp_indices)
	M = idx.size

	# Allocate output matrix
	D = np.empty((M, n_samples), dtype=dtype)

	# Temporary buffer for a single variant
	buf = np.empty(n_samples, dtype=dtype)

	# Read dosages one variant at a time
	for j, v_idx in enumerate(idx):
		ra.read_dosages(int(v_idx), buf)
		D[j] = buf

	return D

def standardize_genotype(G, axis=1, eps=1e-12):
	"""
	Standardize a genotype matrix SNP-wise (rows = SNPs, columns = samples).

	Parameters
	----------
	G : np.ndarray (num_snps, num_samples)
		Genotype or dosage matrix.
	axis : int
		Axis along which to compute mean/variance (default: 1 = per SNP).
	eps : float
		Small constant to avoid division by zero.

	Returns
	-------
	G_std : np.ndarray
		Standardized genotype matrix with zero-variance SNPs removed.
	keep_mask : np.ndarray[bool]
		Boolean mask indicating which SNPs were retained.
	"""

	# Compute mean and std per SNP
	mean = G.mean(axis=axis, keepdims=True)
	std  = G.std(axis=axis, keepdims=True)

	# Identify SNPs with nonzero variance
	keep_mask = (std.squeeze() > eps)

	# Subset to variance>0 SNPs
	G_keep = G[keep_mask]

	# Standardize
	mean_keep = mean[keep_mask]
	std_keep  = std[keep_mask]

	G_std = (G_keep - mean_keep) / std_keep

	return G_std, keep_mask

def get_HE_regression_summary_stats_for_single_gene(YY, GG, EE):
	# Standardize EE
	stand_EE = (EE - np.mean(EE))/np.std(EE)
	# Create interaction effects
	EE_GG = GG * stand_EE[:, np.newaxis]


	pred_SS = np.ones(len(EE))
	pred_SS[EE==1.0] = np.var(YY[EE==1.0],ddof=1)
	pred_SS[EE==0.0] = np.var(YY[EE==0.0],ddof=1)		
	pred_SS = pred_SS/np.mean(pred_SS)
	PA_GG = GG * np.sqrt(pred_SS[:, np.newaxis])


	# Covariances
	Y_t_Y = np.dot(YY.reshape(-1,1), YY.reshape(1,-1))
	E_t_E = np.dot(stand_EE.reshape(-1,1), stand_EE.reshape(1,-1))
	G_T_G = np.dot(GG, np.transpose(GG))
	EG_T_EG = np.dot(EE_GG, np.transpose(EE_GG))
	SG_T_SG =np.dot(PA_GG, np.transpose(PA_GG))

	# Get upper-traingle elements in vector form of each covariance matrix
	i_idx, j_idx = np.triu_indices(len(YY), k=1)
	Y_t_Y_vec = Y_t_Y[i_idx, j_idx]
	E_t_E_vec = E_t_E[i_idx, j_idx]
	G_T_G_vec = G_T_G[i_idx, j_idx]
	EG_T_EG_vec = EG_T_EG[i_idx, j_idx]
	SG_T_SG_vec = SG_T_SG[i_idx, j_idx]

	# Put all into one matrix
	X_mat = np.transpose(np.vstack((E_t_E_vec, G_T_G_vec, EG_T_EG_vec, SG_T_SG_vec)))

	# Compute summary stats
	X_t_X = np.dot(np.transpose(X_mat), X_mat)
	X_t_Y = np.dot(np.transpose(X_mat), Y_t_Y_vec)
	ratio = np.square(np.mean(np.sqrt(pred_SS)))
	
	return X_t_X, X_t_Y, ratio

def extract_total_number_of_genes(filer):
	f = open(filer)
	autosomal_chroms = {}
	for chrom_num in range(1,23):
		autosomal_chroms['chr' + str(chrom_num)] =1
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] in autosomal_chroms:
			counter = counter + 1
	f.close()
	return counter

def extract_per_gene_HE_regression_summary_stats(args):

	# First get total number of genes
	total_n_genes = extract_total_number_of_genes(args.expression_bed)

	# Initialize list to keep track of summary stats
	gene_ss_array = []

	# Extract interaction QTL E variable
	EE, E_var_sample_names = load_in_E_variable(args.binary_E_interaction_covariate_file)

	# Make progress bar
	pbar = tqdm(total=total_n_genes, desc="Processing genes")

	# Loop through chromosomes
	for chrom_num in range(1,23):

		# Convert from integer chrom to chrom_string
		chrom_string = 'chr' + str(chrom_num)

		# Load in genotype data for this chromosome
		pgen, pvar, psam = load_in_genotype(args.plink2_per_chrom_stem, chrom_num)
		# Extract positions of variants (all on same chrom)
		variant_positions = pvar[:, 1].astype(float)
		# Check to make sure pgen file has variants on only one chromosome
		if len(np.unique(pvar[:,0])) != 1:
			print('Fatal error: plink2 file contains variants on multiple chromosomes')
			sys.exit(1)
		if np.unique(pvar[:,0])[0] != str(chrom_num):
			print('Fatal error: variants in plink2 file do not patch chromosome of plink2 label')
			sys.exit(1)

		# E_var_sample_names must be a subset of psam.
		# But psam can have more samples (and in different order) than E_var_sample_names
		# Extract vector containing ordering of genotype samples
		genotype_sample_reordering = reorder_genotype_samples(psam, E_var_sample_names)

		# Now loop through genes
		f = open(args.expression_bed)
		head_count = 0
		gene_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')

			# Treat header line seperately
			if head_count == 0:
				head_count = head_count + 1
				expression_sample_names = np.asarray(data[4:])
				# We assume expression samples have same order as E_var_sample_names
				if np.array_equal(expression_sample_names, E_var_sample_names) == False:
					print('Fatal error: Expression sample names do not match interaction variable sample names (need 1 to 1 match)')
				continue

			#########
			# Standard line
			# Parse relevent fields
			gene_chrom_string = data[0]
			gene_tss = float(data[1])
			gene_id = data[3]

			# Skip genes not on current chromosome
			if gene_chrom_string != chrom_string:
				continue

			########
			# Get expression levels for this gene
			YY = np.asarray(data[4:]).astype(float)
			# filter out genes with 0 variance
			if np.std(YY) == 0.0:
				continue
			# Standardize YY 
			if args.standardize_expression:
				YY = (YY - np.mean(YY))/np.std(YY)

			#########
			# Extract cis-genotype matrix for this gene
			# Filter genotype data for this simulation
			cis_snp_indices = (variant_positions >= gene_tss - args.cis_radius) & (variant_positions < gene_tss + args.cis_radius)
			# Skip genes with too few snps
			n_cis_snps = np.sum(cis_snp_indices)
			if n_cis_snps < args.min_snps_per_gene:
				continue
			# Extract cis snp genotype
			cis_genotype_mat = extract_dosages_from_reader(pgen, cis_snp_indices)
			# Filter to sample subset
			cis_genotype_mat = cis_genotype_mat[:, genotype_sample_reordering]
			# Standardize genotype and filter out snps with 0 variance
			standardized_cis_genotype_mat, no_var_filter_mask = standardize_genotype(cis_genotype_mat)
			# updated number of cis snps
			n_cis_snps = standardized_cis_genotype_mat.shape[0]
			GG = np.transpose(standardized_cis_genotype_mat)
			del standardized_cis_genotype_mat

			# Get regression summary stats for a single gene
			gene_XtX, gene_XtY, gene_ratio = get_HE_regression_summary_stats_for_single_gene(YY, GG, EE)


			# Add to array
			gene_ss_array.append((gene_id, gene_chrom_string, gene_tss, n_cis_snps, gene_XtX, gene_XtY, gene_ratio))

			# Update progress bar
			pbar.update(1)
		f.close()
	return gene_ss_array


def compute_variance_parameters(per_gene_HE_ss):
	global_XtX = None
	global_XtY = None
	gene_nsnps = []
	gene_ratios = []

	# Loop through genes
	for gene_tuple in per_gene_HE_ss:
		gene_nsnps.append(gene_tuple[3])
		gene_ratios.append(gene_tuple[6])
		if global_XtX is None:
			global_XtX = np.copy(gene_tuple[4])
			global_XtY = np.copy(gene_tuple[5])
		else:
			global_XtX = global_XtX + gene_tuple[4]
			global_XtY = global_XtY + gene_tuple[5]

	# Compute regression coefficints
	coef = np.dot(np.linalg.inv(global_XtX), global_XtY)

	# Organize into arrays
	gene_ratios = np.asarray(gene_ratios)
	gene_nsnps = np.asarray(gene_nsnps)
	n_genes = len(per_gene_HE_ss)

	# Get gene level estimates
	gene_constant_genetic_variance = coef[1]*np.sum(gene_nsnps)/n_genes
	gene_traditional_interaction_genetic_variance = coef[2]*np.sum(gene_nsnps)/n_genes
	gene_PA_constant_genetic_variance = coef[3]*np.sum(gene_ratios*gene_nsnps)/n_genes
	gene_PA_interaction_genetic_variance = coef[3]*np.sum((1.0-gene_ratios)*gene_nsnps)/n_genes

	gene_total_constant_genetic_variance = gene_constant_genetic_variance + gene_PA_constant_genetic_variance
	
	return gene_constant_genetic_variance, gene_PA_constant_genetic_variance, gene_PA_interaction_genetic_variance, gene_traditional_interaction_genetic_variance, gene_total_constant_genetic_variance


def compute_standard_interaction_variance_parameters(per_gene_HE_ss):
	global_XtX = None
	global_XtY = None
	gene_nsnps = []

	# Loop through genes
	for gene_tuple in per_gene_HE_ss:
		gene_nsnps.append(gene_tuple[3])

		gene_XtX = gene_tuple[4]
		gene_XtY = gene_tuple[5]

		# Remove PA term from gene_XtX
		gene_XtX_no_PA = gene_XtX[:3, :][:, :3]
		gene_XtY_no_PA = gene_XtY[:3]

		if global_XtX is None:
			global_XtX = np.copy(gene_XtX_no_PA)
			global_XtY = np.copy(gene_XtY_no_PA)
		else:
			global_XtX = global_XtX + gene_XtX_no_PA
			global_XtY = global_XtY + gene_XtY_no_PA

	# Compute regression coefficints
	coef = np.dot(np.linalg.inv(global_XtX), global_XtY)

	# Organize into arrays
	gene_nsnps = np.asarray(gene_nsnps)
	n_genes = len(per_gene_HE_ss)

	# Get gene level estimates
	gene_constant_genetic_variance = coef[1]*np.sum(gene_nsnps)/n_genes
	gene_traditional_interaction_genetic_variance = coef[2]*np.sum(gene_nsnps)/n_genes

	return gene_constant_genetic_variance, gene_traditional_interaction_genetic_variance


def bootstrap_variance_parameters(per_gene_HE_ss, n_boots=1000):
	# Keep track of bootsrapped values
	bs_mat = []
	for bs_iter in range(n_boots):
		# Randomly sample genes with replacement
		bs_per_gene_HE_ss = random.choices(per_gene_HE_ss, k=len(per_gene_HE_ss)) 

		# Run genome wide regression on bootstrapped samples
		bs_constant_genetic_variance, bs_pa_constant_genetic_variance, bs_pa_interaction_genetic_variance, bs_traditional_interaction_genetic_variance, bs_total_constant_genetic_variance = compute_variance_parameters(bs_per_gene_HE_ss)
		
		# Add to global bootstrap mat
		bs_mat.append([bs_constant_genetic_variance, bs_pa_constant_genetic_variance, bs_pa_interaction_genetic_variance, bs_traditional_interaction_genetic_variance, bs_total_constant_genetic_variance])

	bs_mat = np.asarray(bs_mat)

	bs_means = np.mean(bs_mat,axis=0)
	bs_ses = np.std(bs_mat,axis=0, ddof=1)

	return bs_means[0], bs_ses[0], bs_means[1], bs_ses[1], bs_means[2], bs_ses[2], bs_means[3], bs_ses[3], bs_means[4], bs_ses[4]


def bootstrap_standard_interaction_variance_parameters(per_gene_HE_ss, n_boots=1000):
	# Keep track of bootsrapped values
	bs_mat = []
	for bs_iter in range(n_boots):
		# Randomly sample genes with replacement
		bs_per_gene_HE_ss = random.choices(per_gene_HE_ss, k=len(per_gene_HE_ss)) 

		# Run genome wide regression on bootstrapped samples
		bs_constant_genetic_variance, bs_traditional_interaction_genetic_variance = compute_standard_interaction_variance_parameters(bs_per_gene_HE_ss)
		
		# Add to global bootstrap mat
		bs_mat.append([bs_constant_genetic_variance, bs_traditional_interaction_genetic_variance])

	bs_mat = np.asarray(bs_mat)

	bs_means = np.mean(bs_mat,axis=0)
	bs_ses = np.std(bs_mat,axis=0, ddof=1)

	return bs_means[0], bs_ses[0], bs_means[1], bs_ses[1]

def print_pah2_cat():
	"""Print a fun ASCII cat logo for PA-h2."""
	print(r"""
 /\_/\
( o.o )   PA-h2
 > ^ <
	""")
	return

def print_pah2_bear():
	# Credit: https://github.com/ivolo/animals/blob/master/data/animals.txt
	print(r"""
   _,-""`""-~`)
(`~           \
 |     a   a   \
 ;        o     ; ___  _,,,,_     _.-~'.
  \      `^`    /`_.-"~      `~-;`      \
   \_      _  .'                 `,     |
     |`-                           \'__/
    /                      ,_       \  `'-.
   /    .-""~~--.            `"-,   ;_    /
  |              \               \  | `""`
   \__.--'`"-.   /_               |'
              `"`  `~~~---..,     |
 PA-h2                       \ _.-'`-.
                              \       \
                               '.     /
                                 `"~"`
	""")
	return
################################################################################
# main
################################################################################
def main():
	######################
	# Command line args
	######################
	# Necessary
	parser = argparse.ArgumentParser()
	parser.add_argument('--expression-bed', default='', type=str,
						help='Bed file containing gene expression')
	parser.add_argument('--plink2-per-chrom-stem', default='', type=str,
						help='Path to genotype stem')
	parser.add_argument('--binary-E-interaction-covariate-file', default='', type=str,
						help='Path to binary E covariate file')
	parser.add_argument('--output-stem', default='', type=str,
						help='Path to output file stem')

	# Defaults
	parser.add_argument('--cis-radius', default=500000, type=int,
						help='cis window around TSS to consider snps')
	parser.add_argument('--min-snps-per-gene', default=50, type=int,
						help='Minimum number of snps per gene. else we throw out the gene.')
	parser.add_argument('--standardize-expression', default=True, type=bool,
						help='Minimum number of snps per gene. else we throw out the gene.')
	args = parser.parse_args()

	print_pah2_bear()
	######################
	# Load in data
	# Create a list of length number of genes
	# Where each element of list contains (gene_name, gene_chrom, gene_position, n_cis_snps, XTX, XTY, gene_ratio)
	# Where XTX and XTY correspond to HE regression summary stats
	######################
	per_gene_HE_ss = extract_per_gene_HE_regression_summary_stats(args)
	'''
	# Temp saving
	f = open("mydata.pkl", "wb")
	pickle.dump(per_gene_HE_ss, f)
	f.close()
	'''
	'''
	# Temp Loading
	f = open("mydata.pkl", "rb")
	per_gene_HE_ss = pickle.load(f)
	f.close()
	'''

	######################
	# Run genome wide PA-h2 regression from per gene summary stats
	######################
	print('\n')
	print('Running genome-wide regression + bootstrapping\n')
	constant_genetic_variance, pa_constant_genetic_variance, pa_interaction_genetic_variance, traditional_interaction_genetic_variance, total_constant_genetic_variance = compute_variance_parameters(per_gene_HE_ss)
	bs_cgv_mean, bs_cgv_se, bs_pa_cgv_mean, bs_pa_cgv_se, bs_pa_igv_mean, bs_pa_igv_se, bs_t_igv_mean, bs_t_igv_se, bs_tot_cgv_mean, bs_tot_cgv_se = bootstrap_variance_parameters(per_gene_HE_ss, n_boots=5000)

	# Print to output
	t = open(args.output_stem + '_PA_h2_summary.txt','w')
	t.write('variance_class\tvariance_estimate\tstandard_error\tZ_score\n')
	t.write('constant_genetic_variance\t' + str(constant_genetic_variance) + '\t' + str(bs_cgv_se) + '\t' + str(constant_genetic_variance/bs_cgv_se) + '\n')
	t.write('PA_constant_genetic_variance\t' + str(pa_constant_genetic_variance) + '\t' + str(bs_pa_cgv_se) + '\t' + str(pa_constant_genetic_variance/bs_pa_cgv_se) + '\n')
	t.write('PA_interaction_genetic_variance\t' + str(pa_interaction_genetic_variance) + '\t' + str(bs_pa_igv_se) + '\t' + str(pa_interaction_genetic_variance/bs_pa_igv_se) + '\n')
	t.write('traditional_interaction_genetic_variance\t' + str(traditional_interaction_genetic_variance) + '\t' + str(bs_t_igv_se) + '\t' + str(traditional_interaction_genetic_variance/bs_t_igv_se) + '\n')
	t.write('total_constant_genetic_variance\t' + str(total_constant_genetic_variance) + '\t' + str(bs_tot_cgv_se) + '\t' + str(total_constant_genetic_variance/bs_tot_cgv_se) + '\n')
	t.close()

	######################
	# Run genome wide standard interaction h2 regression from per gene summary stats
	######################
	std_constant_genetic_variance, std_interaction_genetic_variance = compute_standard_interaction_variance_parameters(per_gene_HE_ss)
	bs_std_cgv_mean, bs_std_cgv_se, bs_std_igv_mean, bs_std_igv_se = bootstrap_standard_interaction_variance_parameters(per_gene_HE_ss, n_boots=5000)

	# Print to output
	t = open(args.output_stem + '_standard_interaction_h2_summary.txt','w')
	t.write('variance_class\tvariance_estimate\tstandard_error\tZ_score\n')
	t.write('constant_genetic_variance\t' + str(std_constant_genetic_variance) + '\t' + str(bs_std_cgv_se) + '\t' + str(std_constant_genetic_variance/bs_std_cgv_se) + '\n')
	t.write('traditional_interaction_genetic_variance\t' + str(std_interaction_genetic_variance) + '\t' + str(bs_std_igv_se) + '\t' + str(std_interaction_genetic_variance/bs_std_igv_se) + '\n')
	t.close()

	print('PA-h2 summary file: ' + args.output_stem + '_PA_h2_summary.txt\n')
	print('standard interaction h2 summary file: ' + args.output_stem + '_standard_interaction_h2_summary.txt')

	return



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
	main()