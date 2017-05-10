#     file: kataegis_splitter.py
#   author: Jesse Eaton and Jacob West-Roberts
#  created: April 29, 2017
# modified: April 29, 2017
#  purpose: For splitting RNA sequence data to kataegis positive and kataegis negative samples


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import numpy as np # for manipulating matricies


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

STR_LENGTH   = 20 # set maximum length of a string in a numpy array


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

#   input: rnas (np.array) [num_rna_samples+1, num_genes] RNA sequence data with header
#          kats (np.array) [num_kat_samples, 2] has TCGA barcodes and q-value
#  output: rnas_out (np.array) [num_same_barcode_samples, num_genes]
#          rnas_out (np.array) [num_same_barcode_samples, 2]
# putpose: keeps all elements in each array that have the same barcode
def keep_same_barcode(kats, rnas):
	rnas_header = rnas[0]
	rnas = np.delete(rnas, 0, 0) # remove header

	kats = kats[np.argsort(kats[:, 0])]
	rnas = rnas[np.argsort(rnas[:, 0])]

	m, n = len(kats), len(rnas)
	i, j = 0, 0
	keep_kats, keep_rnas = np.zeros(m, dtype = bool), np.zeros(n, dtype = bool)
	while True:
		if i >= m:
			while j < n:
				keep_rnas[j] = False
				j += 1
			break
		if j >= n:
			while i < m:
				keep_kats[i] = False
				i += 1
			break
		kat_bar, rna_bar = kats[i, 0], rnas[j, 0]
		if kat_bar == rna_bar:
			keep_kats[i], keep_rnas[j] = True, True
			i += 1
			j += 1
		elif kat_bar < rna_bar:
			keep_kats[i] = False
			i += 1
		else:
			keep_rnas[j] = False
			j += 1
	kats_out = kats[keep_kats]
	rnas_out = np.insert(rnas[keep_rnas], 0, rnas_header, axis = 0)
	return kats_out, rnas_out

#  input: rnas (np.array) [num_samples+1, num_genes] RNA sequence data with header
#         kats (np.array) [num_samples, 2] has TCGA barcodes and q-value
#         cutoff (float) q-value enrichment cutoff to separate kataegis pos from neg
# output: rnas_pos (np.array) [num_samples_kataegis_positive]
#         rnas_neg (np.array) [num_samples_kataegis_negative]
def kat_split(rnas, kats, cutoff):
	rnas = np.delete(rnas, 0, axis = 0) # remove first (header) row
	num_genes = rnas.shape[1]

	# dictionary. key: barcode, val: True/Fals for Kat pos/neg
	has_kat = get_has_kat_dic(kats, cutoff)
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
#         cutoff (float) q-value enrichment cutoff to separate kataegis pos from neg
# output: dic (dict) key is TCGA barcode. val is boolean True if q-value <= 0.05
#                                                        False if q-value > 00.05
def get_has_kat_dic(kats, cutoff):
	dic = {}
	for barcode, q_val in kats:
		if float(q_val) <= cutoff:
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


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
