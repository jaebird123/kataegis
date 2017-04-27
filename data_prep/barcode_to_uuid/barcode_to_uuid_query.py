#     file: barcode_to_uuid.py
#   author: Jesse Eaton
#  created: April 26, 2017
# modified: April 26, 2017
#  purpose: changes a comma separated list of TCGA barcodes to a comma separated list of UUIDS


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

LOOKUP_FILE = 'barcode_to_uuid.tsv'
INPUT_DELIMINATOR = ','

QUERY_KEY_STRING = 'cases.case_id = '
QUERY_DELIMINATOR = ' OR '


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	to_uuid = get_uuid_dict(LOOKUP_FILE)

	uuids = []
	barcodes = []
	barcodes_not_found = []
	for line in args['input_file']:
		for barcode in line.split(INPUT_DELIMINATOR):
			if barcode in to_uuid:
				uuids.append(to_uuid[barcode])
				barcodes.append(barcode)
			else:
				barcodes_not_found.append(barcode)

	if args['map_output']:
		print_barcodes_uuids(barcodes, uuids, barcodes_not_found)
	else:
		print_uuids_query(uuids, args['limit_query_size'])

# returns dict. key are TCGA barcode and vals are UUIDs
def get_uuid_dict(fname):
	dic = {}
	with open(fname) as f:
		next(f) # skip the first line (it is header information)
		for line in f:
			cols = line.split('\t')
			barcode = cols[26] # TCGA barcodes are in the 26th column (zero indexed)
			uuid = cols[29]    # UUIDs are in the 29th column (zero indexed)
			dic[barcode] = uuid
	return dic

def print_barcodes_uuids(barcodes, uuids, barcodes_not_found):
	s = ''
	for i, barcode in enumerate(barcodes):
		s += barcode + '\t' + uuids[i] + '\n'
	print s
	print 'did not find ' + str(len(barcodes_not_found)) + ' barcode(s):'
	for barcode in barcodes_not_found:
		print barcode

def print_uuids_query(uuids, num_uuids_per_query):
	# add 'cases.case_id = ' infront of each UUID
	for i, uuid in enumerate(uuids):
		uuids[i] = QUERY_KEY_STRING + uuid

	# put all UUIDs in one query if no limit specified
	if not num_uuids_per_query:
		num_uuids_per_query = len(uuids)

	query_lists = chunkify_list(uuids, num_uuids_per_query)
	for uuids_in_query in query_lists:
		print QUERY_DELIMINATOR.join(uuids_in_query)

	print '\nYou must perform ' + str(len(query_lists)) + ' queries to get all ' + str(len(uuids)) + ' files.'

#  input: data (list of anytype)
#         chunk_size (int)
# output: chunks (2D list of input type) input list split into multiple lists
#           each of size chunk_size. last list may be any size <= chunk_size
def chunkify_list(data, chunk_size):
	return [ data[x : x+chunk_size] for x in xrange(0, len(data), chunk_size) ]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'barcode_to_uuid.py', description = "changes a comma separated list of TCGA barcodes to a comma separated list of UUIDS")
	parser.add_argument('input_file', help = 'input .txt file', type = lambda x: is_valid_file(parser, x))
	parser.add_argument('-m', '--map_output', help = 'changes output format to include original TCGA barcode with UUIDs', action = 'store_true')
	parser.add_argument('-l', '--limit_query_size', type = int, help = 'limit the number of UUIDs in each query. If this number is lower than the number of UUIDs multiple queries will be outputed')
	parsed_args = parser.parse_args(argv)
	if (parsed_args.limit_query_size) and (parsed_args.map_output):
		parser.error('should not use --map_output and --limit_query_size together')
	return vars(parsed_args)

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
