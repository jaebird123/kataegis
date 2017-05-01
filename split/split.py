#     file: split.py
#   author: Jesse Eaton
#  created: April 30, 2017
# modified: April 30, 2017
#  purpose: Splits RNA sequence data to kataegis positive and kataegis negative files


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import numpy as np # for manipulating matricies

# local modules
sys.path.insert(0, '../helper/')
import kataegis_splitter as ks # for splitting RNA seq data to kataegis pos and neg samples


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

# here


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

	# insert file header (with gene names)
	rnas_pos = np.insert(rnas_pos, 0, rnas[0], 0)
	rnas_neg = np.insert(rnas_neg, 0, rnas[0], 0)

	# get file name with no directory to file and no .txt extension
	fname_prefix = os.path.splitext(os.path.basename(args['rna_seq_file'].name))[0]
	fname_pos = args['output_directory'] + fname_prefix + '.kat_pos.txt'
	fname_neg = args['output_directory'] + fname_prefix + '.kat_neg.txt'

	touch(fname_pos)
	touch(fname_neg)

	np.savetxt(open(fname_pos, 'w'), rnas_pos, fmt = '%s', delimiter = '\t', newline = '\n')
	np.savetxt(open(fname_neg, 'w'), rnas_neg, fmt = '%s', delimiter = '\t', newline = '\n')

# creates file if file does not already exist
def touch(fname, times = None):
	with open(fname, 'a'):
		os.utime(fname, times)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'barcode_and_q_values.py', description = "creates a tab separated value sheet with TCGA barcodes and q-value enrichment")
	parser.add_argument('-r', '--rna_seq_file', type = lambda x: is_valid_file(parser, x), required = True, help = 'file containing all RNA sequence data. should contain all values for each RNA for each sample in TCGA BRCA')
	parser.add_argument('-k', '--kataegis_file', type = lambda x: is_valid_file(parser, x), required = True, help = 'file containing q-value enrichment for each sample')
	parser.add_argument('-o', '--output_directory', type = lambda x: valid_directory(parser, x), required = True, help = 'directory for two output RNA sequence files to go')
	parser.add_argument('-q', '--q_value_cutoff', type = lambda x: bounded_float(parser, x, 0.0, 1.0), default = 0.05, help = 'float between 0.0 and 1.0. any samples with q value <= cutoff are kataegis positive. q value > cutoff are kataegis negative')
	return vars(parser.parse_args(argv))

def is_valid_file(parser, arg):
	if not os.path.exists(arg):
		parser.error('The file \"' + str(arg) + '\" could not be found.')
	else:
		return open(arg, 'r')

# returns string as directory. adds error to parser if no valid directory
def valid_directory(parser, arg):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	return directorize(arg)

# add "/" to end of directory name if necessary
def directorize(dir_name):
	if dir_name.endswith('/'):
		return dir_name
	return dir_name + '/'


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
