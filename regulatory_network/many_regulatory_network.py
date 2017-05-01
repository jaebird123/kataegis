#     file: many_regulatory_network.py
#   author: Jesse Eaton
#  created: May 1, 2017
# modified: May 1, 2017
#  purpose: Iteratively calls regulatory_network.py on many directories


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys        # for command line arguments
import os         # for manipulating files and folders
import argparse   # for command line arguments
import subprocess # for calling other python scripts


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

IN_KAT_EXTS = ['.kat_pos.txt', '.kat_neg.txt']
OUT_EXTS = ['.panda.pairs', '.lion.pairs']
REG_NET_RUN_FILE = 'regulatory_network.py'
GENE_TOP_NUM = 100



# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)

	in_file_dict = get_file_dict(args['input_directory'], IN_KAT_EXTS)
	out_file_dict = touch_output_directories_and_files(args['output_directory'], args['input_directory'], in_file_dict, OUT_EXTS)

	num_times_run = len(in_file_dict) * 2
	print_now('Running ' + REG_NET_RUN_FILE + ' ' + str(num_times_run) + ' times:')
	counter = 0
	for subdir, in_files in in_file_dict.iteritems():
		for kat_ext in ['.kat_pos', '.kat_neg']:
			counter += 1
			in_file = elements_with(in_files, kat_ext)[0]
			out_files = elements_with(out_file_dict[subdir], kat_ext)
			panda_file = elements_with(out_files, 'panda')[0]
			lion_file = elements_with(out_files, 'lion')[0]
			print_now('\n\t' + str(counter) + ' of ' + str(num_times_run) + ' - ' + in_file)
			run_reg_net(args['input_directory'] + subdir + in_file, args['output_directory'] + subdir + panda_file, args['output_directory'] + subdir + lion_file, GENE_TOP_NUM)

# def run_reg_net(in_file, out_file):
def run_reg_net(in_file, panda_file, lion_file, top_num):
	cmd = ' '.join(['python', REG_NET_RUN_FILE, in_file, '-p', panda_file, '-l', lion_file, '-t', str(top_num)])
	try:
		subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT) # run build_phylogeny command
	except subprocess.CalledProcessError as e:
		eprint('\nAn error in ' + REG_NET_RUN_FILE + ' occured:\n') # error if error occured while running script
		eprint(e.output)
		exit(e.returncode)

#  input: arr (list of string)
#         substr (string) substring
# output: out_arr (list of string) any elements in arr containing substring
def elements_with(arr, substr):
	out_arr = []
	for elem in arr:
		if substr in elem:
			out_arr.append(elem)
	return out_arr


#  input: in_dir (string) name of input directory
#         exts (list of string) file extensions. should including '.' in string. ex. '.vcf'
# output: dic (dictionary) keys are subdirectories (folder names) of input directory
#                          vals (list of strings) files with ext extension inside subdirectory
def get_file_dict(in_dir, exts):
	dic = {}
	for subdir, dirs, files in walklevel(in_dir):
		for d in dirs:
			dic[directorize(d)] = []
			for ext in exts:
				dic[directorize(d)].append(files_with_extension(in_dir + d, ext)[0])
	return dic

#  input: out_dir (string) directory to make subdirectories of same name to subdirectories in in_dir
#         in_dir (string) master directory for bulk input. ex: 'data/bulk/bulk_output/'
#         in_file_dict (dict) keys are (string) subdirectory names
#                             vals are (string) input file name in subdirectories
#         out_exts (list of string) extensions of output files to touch in each output subdirectory. ex: ['.xml']
#   does: creates all output subdirectories. touches all output files in subdirectories.
# output: out_file_dict (dict) keys are (string) output subdirectory names
#                              vals are (string) output file names in subdirectories
def touch_output_directories_and_files(out_dir, in_dir, in_file_dict, out_exts):
	out_file_dict = {}
	for in_subdir, in_files in in_file_dict.iteritems():
		out_file_dict[in_subdir] = []
		for in_file in in_files:
			for out_ext in out_exts:
				out_file = os.path.splitext(in_file)[0] + out_ext
				if not os.path.exists(out_dir + in_subdir):
					os.makedirs(out_dir + in_subdir)
				touch(out_dir + in_subdir + out_file)
				out_file_dict[in_subdir].append(out_file)
	return out_file_dict

# creates file if file does not already exist
def touch(fname, times = None):
	with open(fname, 'a'):
		os.utime(fname, times)

# prints when called (not after script is finished running)
def print_now(s):
	sys.stdout.write(s)
	sys.stdout.flush()

# print to standard error
def eprint(s):
	temp = sys.stdout
	sys.stdout = sys.stderr
	print s
	sys.stdout = temp


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'barcode_and_q_values.py', description = "creates a tab separated value sheet with TCGA barcodes and q-value enrichment")
	parser.add_argument('-i', '--input_directory', type = lambda x: valid_master_directory(parser, x, IN_KAT_EXTS), help = 'directory containing subdirectories each with .kat_pos.txt and .kat_neg.txt files inside', required = True)
	parser.add_argument('-o', '--output_directory', type = lambda x: valid_directory(parser, x), help = 'directory where output data will go', required = True)
	parser.add_argument('-t', '--top_genes_plot', type = lambda x: int_between(parser, x, 1, 500), default = 100, help = 'number of top genes to plot a gene regulatory network for. will only plot PANDA gene regulatory network')
	return vars(parser.parse_args(argv))

# returns directory name with "/" suffix if directory has subdirectories with a single file from with extension from exts
# errors otherwise
def valid_master_directory(parser, arg, exts):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	arg = directorize(arg)
	for subdir, dirs, files in walklevel(arg):
		if not dirs:
			parser.error('The directory \"' + str(subdir) + '\" has no subdirectories.')
		for d in dirs:
			for ext in exts:
				d_files = files_with_extension(arg + d, ext)
				if len(d_files) != 1:
					parser.error('The subdirectory \"' + arg + d + '\" should have exactly one ' + ext + ' file.')
	return arg

# returns all files in the directory with extension ext (ex. ext = '.vcf')
def files_with_extension(directory, ext):
	files = []
	for file in os.listdir(directory):
		if file.endswith(ext):
			files.append(file)
	return files

# os.walk but only goes recursively down level levels
def walklevel(some_dir, level = 0):
	some_dir = some_dir.rstrip(os.path.sep)
	assert os.path.isdir(some_dir)
	num_sep = some_dir.count(os.path.sep)
	for root, dirs, files in os.walk(some_dir):
		yield root, dirs, files
		num_sep_this = root.count(os.path.sep)
		if num_sep + level <= num_sep_this:
			del dirs[:]

# returns string as directory. adds error to parser if no valid directory
def valid_directory(parser, arg):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	return directorize(arg)

def is_valid_file(parser, arg):
	if not os.path.exists(arg):
		parser.error('The file \"' + str(arg) + '\" could not be found.')
	else:
		return open(arg, 'r')

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
