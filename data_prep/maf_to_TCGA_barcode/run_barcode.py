#     file: run_TCGA_barcode.py
#   author: Jesse Eaton
#  created: April 9, 2017
# modified: April 9, 2017
#  purpose: uses the Broad Institute's MAF sorted sum of all Fisher correlation values
#             for each sample and prints the sample barcodes to standard output


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

TCGA_PREFIX = 'TCGA'   # prefix of string to identify a TCGA barcode in input file
NUM_IDS_IN_BARCODE = 3 # number of first '-' separated sections we say is the barcode.
                       #   there is more but we ignore the rest after the NUM_IDS_IN_BARCODE occurance of a '-'
DELIMINATOR = ','      # character used to separate each barcode in output


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	barcodes = get_barcodes(args['input_file'])
	print DELIMINATOR.join(barcodes)

#  input: fname (string) name of input file (should have .txt extension)
# output: barcodes (list of strings) TCGA barcodes for each sample listed in input file
def get_barcodes(file):
	barcodes = []
	for line in file:
		# split by tabs \t, get first column [0], and remove quotes using translate
		first_col_val = line.split('\t')[0].translate(None, '\"')
		if TCGA_PREFIX in first_col_val:
			barcodes.append(split_nth_occurance(first_col_val, NUM_IDS_IN_BARCODE, '-'))
	return barcodes

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
	parser = argparse.ArgumentParser(prog = 'run_TCGA_barcode.py', description = "uses the Broad Institute's MAF sorted sum of all Fisher correlation values for each sample and prints the sample barcodes to standard output")
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
