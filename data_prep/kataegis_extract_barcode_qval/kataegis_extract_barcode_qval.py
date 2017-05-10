#     file: kataegis_extract_barcode_qval.py
#   author: Jesse Eaton and Jacob West-Roberts
#  created: April 27, 2017
# modified: April 27, 2017
#  purpose: Outputs only cleaned TCGA barcode and q value enrichment


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

KAT_FILE_BARCODE_COL = 0  # column with TCGA barcode
KAT_FILE_Q_VAL_COL   = 58 # column that contains q values of CG enrichment
NUM_IDS_IN_BARCODE   = 3  # number of '-' separated strings to keep in TCGA barcode. keep first 3


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	
	s = ''
	for line in args['input_file']:
		cols = line.split('\t')
		barcode = clean_barcode(cols[KAT_FILE_BARCODE_COL])
		q_val = cols[KAT_FILE_Q_VAL_COL]
		if 'TCGA' in barcode: # exclude header and footer rows
			s += barcode + '\t' + q_val + '\n'
	print s

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
	parser = argparse.ArgumentParser(prog = 'barcode_and_q_values.py', description = "creates a tab separated value sheet with TCGA barcodes and q-value enrichment")
	parser.add_argument('input_file', help = 'input .txt file', type = lambda x: is_valid_file(parser, x))
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
