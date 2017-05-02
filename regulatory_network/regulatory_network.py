#     file: regulatory_network.py
#   author: Jesse Eaton
#  created: May 1, 2017
# modified: May 1, 2017
#  purpose: Create gene regulatory network for each sample from RNA sequence input


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys         # for command line arguments
import os          # for manipulating files and folders
import argparse    # for command line arguments
import numpy as np # for manipulating matricies
from pypanda import Panda          # for inferring a single gene regulatory network for all samples
from pypanda import Lioness        # for inferring gene regulatory networks for each sample
from pypanda import AnalyzePanda   # for plotting gene regulatory network from PANDA
from pypanda import AnalyzeLioness # for plotting gene regulatory network from LIONESS
import pandas as pd # haha I kind of don't know. maybe remove


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

TEMP_DIR = 'tmp/'
TEMP_PANDA_INPUT = 'expression_data.txt'


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	rnas = np.genfromtxt(args['input_file'], dtype = str, delimiter = '\t')

	# write input files for Panda
	fname = write_panda_input(rnas, TEMP_DIR, TEMP_PANDA_INPUT)

	# run Panda (create gene regulatory network)
	p = Panda(fname, None, None, remove_missing = False)

	# save Panda results
	p.save_panda_results(file = args['panda_output_file'])

	if args['top_genes_plot']:
		num_genes = args['top_genes_plot']
		plot = AnalyzePanda(p)
		plot_fname = os.path.splitext(args['panda_output_file'])[0] + '.top_' + str(num_genes) + '_genes.png'
		plot.top_network_plot(top = num_genes, file = plot_fname)

	if args['lion_output_file']:

		# run Lioness (infer many gene regulatory networks. one for each sample)
		l = Lioness(p)

		# save Lioness
		l.save_lioness_results(file = args['lion_output_file'])

	# plot = AnalyzeLioness(l)
	# plot.top_network_plot(column= 0, top = 100, file = 'top_100_genes.png')

def write_panda_input(rnas, temp_dir, fname):
	rnas = np.transpose(np.delete(rnas, 0, 1))
	if not os.path.exists(temp_dir):
		os.makedirs(temp_dir)
	temp_fname = temp_dir + fname
	touch(temp_fname)
	np.savetxt(open(temp_fname, 'w'), rnas, fmt = '%s', delimiter = '\t', newline = '\n')
	return temp_fname

# creates file if file does not already exist
def touch(fname, times = None):
	with open(fname, 'a'):
		os.utime(fname, times)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'barcode_and_q_values.py', description = "creates a tab separated value sheet with TCGA barcodes and q-value enrichment")
	parser.add_argument('input_file', help = '.txt file containing all RNA sequence data. should contain all values for each RNA for each sample in TCGA BRCA', type = lambda x: is_valid_file(parser, x))
	parser.add_argument('-p', '--panda_output_file', type = lambda x: has_extension(parser, x, '.panda.pairs'), required = True, help = 'name of file for PANDA output to go. should have .panda.pairs extension')
	parser.add_argument('-l', '--lion_output_file', type = lambda x: has_extension(parser, x, '.lion.pairs'), help = 'name of file for LIONESS output to go. should have .lion.pairs extension')
	parser.add_argument('-t', '--top_genes_plot', type = lambda x: int_between(parser, x, 1, 500), help = 'number of top genes to plot a gene regulatory network for. will only plot PANDA gene regulatory network')
	return vars(parser.parse_args(argv))

def int_between(parser, arg, low, high):
	err_msg = 'should be integer between ' + str(low) + ' and ' + str(high)
	try:
		arg = int(arg)
	except:
		parser.error(err_msg)
	if not (low <= arg and arg <= high):
		parser.error(err_msg)
	return arg

def has_extension(parser, arg, ext):
	if not arg.endswith(ext):
		parser.error('file should end with ' + ext + ' extension')
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
