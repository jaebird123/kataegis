#     file: plot_expression_diff.py
#   author: Jesse Eaton
#  created: April 28, 2017
# modified: April 28, 2017
#  purpose: plots histogram of number of genes VS difference in mean expression levels between
#             kataegis positive and kataegis negative samples


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import numpy as np              # for manipulating matricies
import matplotlib.pyplot as plt # for plotting histogram
import matplotlib.mlab as mlab


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

Q_VAL_CUTOFF = 0.05
STR_LENGTH   = 20   # set maximum length of a string in a numpy array


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	
	kats = np.genfromtxt(args['kataegis_file'], dtype = str, delimiter = '\t')
	rnas = np.genfromtxt(args['rna_seq_file'], dtype = str, delimiter = '\t')

	rnas_pos, rnas_neg = split_by_kat(rnas, kats)

	expr_diffs = get_expr_diffs(rnas_pos, rnas_neg)

	buckets = np.histogram(expr_diffs, bins = 100, range = [-200, 200])[0]
	plt.plot(buckets)
	plt.show()

#  input: rnas (np.array) [num_samples+1, num_genes]
#         kats (np.array) [num_samples, 2] has TCGA barcodes and q-value
# output: rnas_pos (np.array) [num_samples_kataegis_positive]
#         rnas_neg (np.array) [num_samples_kataegis_negative]
def split_by_kat(rnas, kats):
	rnas = np.delete(rnas, 0, axis = 0) # remove first (header) row
	num_genes = rnas.shape[1]

	# dictionary. key: barcode, val: True/Fals for Kat pos/neg
	has_kat = get_has_kat_dic(kats)
	num_pos, num_neg = count_pos_neg(has_kat) # count number of kataegis pos vs neg

	# create two separate arrays. one for kat pos. other for kat neg
	rnas_pos = np.zeros([num_pos, num_genes], dtype = 'S' + str(STR_LENGTH))
	rnas_neg = np.zeros([num_neg, num_genes], dtype = 'S' + str(STR_LENGTH))
	i, j = 0, 0
	for sample in rnas:
		barcode = sample[0]
		if barcode in has_kat:
			if has_kat[barcode]:
				rnas_pos[i, :] = sample
				i += 1
			else:
				rnas_neg[j, :] = sample
				j += 1
	
	# remove all zero rows (rows that were not used due to not finding a barcode in RNA seq data)
	return rnas_pos[np.all(rnas_pos != '', axis=1)], rnas_neg[np.all(rnas_neg != '', axis=1)]

#  input: kats (np.array) [num_samples, 2] has TCGA barcodes and q-value
# output: dic (dict) key is TCGA barcode. val is boolean True if q-value <= 0.05
#                                                        False if q-value > 00.05
def get_has_kat_dic(kats):
	dic = {}
	for barcode, q_val in kats:
		if float(q_val) <= Q_VAL_CUTOFF:
			dic[barcode] = True
		else:
			dic[barcode] = False
	return dic

# returns the number of values that are True and number of values that are False
def count_pos_neg(dic):
	num_pos, num_neg = 0, 0
	for key, is_pos in dic.iteritems():
		if is_pos:
			num_pos += 1
		else:
			num_neg += 1
	return num_pos, num_neg

#  input: rnas_pos (np.array) [num_samples_with_kataegis, num_genes]
#         rnas_neg (np.array) [num_samples_w/o_kataegis, num_genes]
# output: expr_diffs (np.array) [num_genes] difference in mean expression betwen kat pos and kat neg
def get_expr_diffs(rnas_pos, rnas_neg):
	pos_genes = np.delete(np.transpose(rnas_pos), 0, axis = 0) # transpose: rows are now genes
	neg_genes = np.delete(np.transpose(rnas_neg), 0, axis = 0) # also remove TCGA header
	num_genes, num_pos = pos_genes.shape
	expr_diffs = np.empty(num_genes, dtype = float)

	# convert to float
	pos_genes = pos_genes.astype(np.float)
	neg_genes = neg_genes.astype(np.float)

	mean_pos = np.mean(pos_genes, axis = 1)
	mean_neg = np.mean(neg_genes, axis = 1)
	return mean_pos - mean_neg


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'plot_expression_diff.py', description = "plots histogram of number of genes VS difference in mean expression levels between kataegis positive and kataegis negative samples")
	parser.add_argument('-r', '--rna_seq_file', type = lambda x: is_valid_file(parser, x), required = True, help = 'file containing all RNA sequence data. should contain all values for each RNA for each sample in TCGA BRCA')
	parser.add_argument('-k', '--kataegis_file', type = lambda x: is_valid_file(parser, x), required = True, help = 'file containing q-value enrichment for each sample')
	return vars(parser.parse_args(argv))

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
