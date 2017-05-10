#     file: rnaseq_filter.py
#   author: Jesse Eaton and Jacob West-Roberts
#  created: April 27, 2017
# modified: April 27, 2017
#  purpose: transposes RNA sequence data so output rows are samples (with a TCGA barcode) and cols
#             are RNAs. removes any sample not seen in kataegis data


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import numpy as np # for manipulating matricies


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

KAT_FILE_BARCODE_COL = 0 # column with TCGA barcode
KAT_FILE_Q_VAL_COL   = 1 # column that contains q values of CG enrichment
NUM_IDS_IN_BARCODE   = 3 # number of '-' separated strings to keep in TCGA barcode. keep first 3
RNA_FILE_BARCODE_COL = 0 # column with TCGA barcode in RNA sequence file


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	
	kats = get_kats(args['kataegis_file'])
	rnaSeq = get_rnaSeq(args['rna_seq_file'])

	rnaSeq = rnaSeq_with_kat_data(rnaSeq, kats)
	np.savetxt(sys.stdout, rnaSeq, fmt = '%s', delimiter = '\t', newline = '\n')

#  input: f (file) kataegis file containing TCGA barcodes and q value enrichment for each sample
#                    along with other values which will be ignored
# output: kats (list of KatSample) each with TCGA barcode and q value
def get_kats(f):
	return np.genfromtxt(f, dtype = str, delimiter = '\t')

#  input: f (file) RNA sequence data. rows are RNAs. cols are TCGA samples
# output: rnaSeq (np.array) rows are TCGA samples. cols are RNAs. vals are normalized expression counts
#           first row and first col are meta information: RNA IDs and TCGA barcodes respectively
def get_rnaSeq(f):
	rnaSeq = np.genfromtxt(f, dtype = str, delimiter = '\t')
	rnaSeq = np.transpose(rnaSeq)    # flip so rows are TCGA samples and cols are different RNAs
	rnaSeq = np.delete(rnaSeq, 1, 1) # remove first column (filled with 'normalized_count')
	rnaSeq[:, RNA_FILE_BARCODE_COL] = clean_barcodes(rnaSeq[:, RNA_FILE_BARCODE_COL])
	return rnaSeq

#  input: rnas (np.array) [num_rna_samples, num_rnas]
#         kats (np.array) [num_kat_samples, 2] cols are [TCGA_barcode, q_val]
# output: rnas_out (np.array) [num_mutual_samples, num_rnas]
def rnaSeq_with_kat_data(rnas, kats):
	rna_header = rnas[0]                 # remove RNA sequence header row. will be replaced later
	rnas = np.delete(rnas, 0, axis = 0)
	kat_barcodes = get_barcode_set(kats) # find all TCGA barcodes found in both data files
	rna_barcodes = get_barcode_set(rnas)
	mut_barcodes = np.intersect1d(kat_barcodes, rna_barcodes)      # mutual barcodes in common
	to_keep = np.in1d(rnas[:, RNA_FILE_BARCODE_COL], mut_barcodes) # mask. find RNA indixes with mut_barcodes
	return np.insert(rnas[to_keep], 0, rna_header, axis = 0)       # keep those RNAs. replace header

#  input: arr (np.array) first column must be TCGA barcodes
# output: barcodes (set)
def get_barcode_set(arr):
	return map(np.unique, arr[:, 0].flatten())


#
#   H E L P E R   F U N C T I O N S
#

# input: barcodes (np.array) 1D (string)
def clean_barcodes(barcodes):
	bar_outs = np.empty_like(barcodes)
	# print bar_outs
	# exit()
	for i, barcode in enumerate(barcodes):
		if 'TCGA' in barcode:
			bar_outs[i] = clean_barcode(barcode)
		else:
			bar_outs[i] = barcode
	return bar_outs

# shortens to first 3 IDs in barcode and removes quotes
def clean_barcode(barcode):
	return split_nth_occurance(barcode.replace('"', ''), NUM_IDS_IN_BARCODE, '-')

# input: s (string) string to split
#        n (int) number of occurances before we split string
#        d (string) deliminator that we look for
def split_nth_occurance(s, n, d):
	groups = s.split(d)
	return d.join(groups[:n])


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'rnaseq_filter.py', description = "transposes RNA sequence data so output rows are samples (with a TCGA barcode) and cols are RNAs. removes any sample not seen in kataegis data")
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
