#     file: gene_reduction.py
#   author: Jesse Eaton
#  created: April 29, 2017
# modified: April 29, 2017
#  purpose: Remove genes we are 1-alpha percent confident are not differentially expressed
#             between kataegis positive and kataegis negative samples


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import numpy as np       # for manipulating matricies
import scipy.stats as st # for calculating Z value of normal(0, 1) distribution

# local modules
sys.path.insert(0, '../helper/')
import kataegis_splitter as ks # for splitting RNA seq data to kataegis pos and neg samples


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

GENE_WHITELIST = ['APOBEC1|339', 'APOBEC2|10930', 'APOBEC3A|200315', 'APOBEC3B|9582',
                  'APOBEC3C|27350', 'APOBEC3D|140564', 'APOBEC3F|200316', 'APOBEC3G|60489',
                  'APOBEC3H|164668', 'APOBEC4|403314']


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)

	kats = np.genfromtxt(args['kataegis_file'], dtype = str, delimiter = '\t')
	rnas = np.genfromtxt(args['rna_seq_file'], dtype = str, delimiter = '\t')

	# keep all samples with TCGA barcodes in both
	kats, rnas = ks.keep_same_barcode(kats, rnas)

	rnas_pos, rnas_neg = ks.kat_split(rnas, kats, args['q_value_cutoff'])
	
	rnas_out = remove_genes(rnas_pos, rnas_neg, rnas[0], args['confidence_level'], GENE_WHITELIST)

	np.savetxt(sys.stdout, rnas_out, fmt = '%s', delimiter = '\t', newline = '\n')

#  input: rnas_pos (np.array) [num_samples_kataegis_positive, num_genes]
#         rnas_neg (np.array) [num_samples_kataegis_negative, num_genes]
#         rnas_header (np.array) [num_genes+1] header information for output RNA sequence data
#         cl (float) (1-alpha) percent confidence level. used to remove non-differentially expressed genes
#         whitelist (list of string) genes that should never be removed
# output: rnas (np.array) [num_samples+1, num_confident_genes+1] with header information
def remove_genes(rnas_pos, rnas_neg, rnas_header, cl, whitelist):
	Z = st.norm.ppf(cl + (1.0-cl)/2.0) # Z-score for standard normal
	col_head_pos = rnas_pos[:, 0] # save column header of TCGA barcodes
	col_head_neg = rnas_neg[:, 0]
	rnas_pos = np.delete(rnas_pos, 0, axis = 1) # remove column header
	rnas_neg = np.delete(rnas_neg, 0, axis = 1)
	gene_avgs_pos = np.mean(rnas_pos.astype(float), axis = 0) # average expression for kataegis pos genes
	gene_avgs_neg = np.mean(rnas_neg.astype(float), axis = 0) #                                 neg
	gene_vars_pos = np.var(rnas_pos.astype(float), axis = 0)  # variance expression for kataegis pos genes
	gene_vars_neg = np.var(rnas_neg.astype(float), axis = 0)  #                                  neg

	gene_avgs_diff = gene_avgs_pos - gene_avgs_neg
	gene_stds_diff = np.sqrt(gene_vars_pos / float(len(rnas_pos)) + gene_vars_neg / float(len(rnas_neg)))
	gene_los = gene_avgs_diff - Z * gene_stds_diff
	gene_his = gene_avgs_diff + Z * gene_stds_diff

	should_remove = np.zeros_like(gene_los)
	for i in xrange(0, len(gene_los)): # for each gene
		if is_between(0, gene_los[i], gene_his[i]) and not rnas_header[i+1] in whitelist:
			should_remove[i] = 1
		else:
			should_remove[i] = 0

	# remove columns that are not differentially expressed
	rnas_pos = np.delete(rnas_pos, np.where(should_remove), axis = 1)
	rnas_neg = np.delete(rnas_neg, np.where(should_remove), axis = 1)
	should_remove = np.insert(should_remove, 0, False) # do not remove TCGA barcode header
	rnas_header = np.delete(rnas_header, np.where(should_remove))

	# insert headers back
	rnas_pos = np.insert(rnas_pos, 0, col_head_pos, axis = 1)
	rnas_neg = np.insert(rnas_neg, 0, col_head_neg, axis = 1)

	# merge kataegis positive, negative, and header for RNA sequence data
	return np.insert(np.insert(rnas_neg, 0, rnas_pos, 0), 0, rnas_header, 0)

def is_between(x, low, high):
	return low <= x and x <= high


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'barcode_and_q_values.py', description = "creates a tab separated value sheet with TCGA barcodes and q-value enrichment")
	parser.add_argument('-r', '--rna_seq_file', type = lambda x: is_valid_file(parser, x), required = True, help = 'file containing all RNA sequence data. should contain all values for each RNA for each sample in TCGA BRCA')
	parser.add_argument('-k', '--kataegis_file', type = lambda x: is_valid_file(parser, x), required = True, help = 'file containing q-value enrichment for each sample')
	parser.add_argument('-c', '--confidence_level', type = lambda x: bounded_float(parser, x, 0.0, 1.0), required = True, help = 'float between 0.0 and 1.0. probability a gene is differentially expressed between kataegis positve and negative samples. Increasing this will remove more genes.')
	parser.add_argument('-q', '--q_value_cutoff', type = lambda x: bounded_float(parser, x, 0.0, 1.0), default = 0.05, help = 'float between 0.0 and 1.0. any samples with q value <= cutoff are kataegis positive. q value > cutoff are kataegis negative')
	return vars(parser.parse_args(argv))

def bounded_float(parser, arg, low, high):
	err_msg = 'confidence level must be float between 0.0 and 1.0'
	try:
		arg = float(arg)
	except:
		parser.error(err_msg)
	if not (arg > low and arg < high):
		parser.error(err_msg)
	return arg

def is_valid_file(parser, arg):
	if not os.path.exists(arg):
		parser.error('The file \"' + str(arg) + '\" could not be found.')
	else:
		return open(arg, 'r')


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
