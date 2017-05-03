#     file: plot_ego_graphs.py
#   author: Jesse Eaton
#  created: May 2, 2017
# modified: May 2, 2017
#  purpose: Plots ego graph (subgraph of node and all neighbors) for set of genes


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import numpy as np    # for manipulating matricies
import networkx as nx # for creating a graph
from operator import itemgetter # TEMPORARY for getting which hub in the graph to display
import matplotlib.pyplot as plt # for plotting and saving plots


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

GENE_NAMES = ['APOBEC1|339', 'APOBEC2|10930', 'APOBEC3A|200315', 'APOBEC3B|9582',
              'APOBEC3C|27350', 'APOBEC3D|140564', 'APOBEC3F|200316', 'APOBEC3G|60489',
              'APOBEC3H|164668', 'APOBEC4|403314']
TOP_HUBS_TO_DISPLAY = 20
THRESHES = np.linspace(0.5, 0.9, 20)


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	pairs = np.genfromtxt(args['input_file'], dtype = str, delimiter = '\t')
	
	# create graph
	graph = get_graph(pairs, args['edge_threshold'])

	for gene_name in GENE_NAMES:
		write_ego_graph(graph, gene_name)

	if args['kat_neg']:
		kat_pos_pairs = pairs
		kat_neg_pairs = np.genfromtxt(args['kat_neg'], dtype = str, delimiter = '\t')
		kat_pos_gene_neighbor_dic = get_gene_neighbor_dic(kat_pos_pairs, THRESHES)
		kat_neg_gene_neighbor_dic = get_gene_neighbor_dic(kat_neg_pairs, THRESHES)
		for gene_name, pos_vals in kat_pos_gene_neighbor_dic.iteritems():
			neg_vals = kat_neg_gene_neighbor_dic[gene_name]
			write_plot_neighbors_vs_edge_weights(pos_vals, neg_vals, THRESHES, gene_name)

	exit()

	# PLOT TOP HUBS
	node_and_degree = graph.degree()
	nodes = sorted(node_and_degree.items(), key = itemgetter(1))[::-1]
	counter = 0
	for gene_name, degree in nodes:
		edges = graph.edges(gene_name, data = True)
		weights = [ round(data['weight'], 3) for _, _, data in edges ]
		num_pos, num_neg = count_pos_neg(weights)
		counter += 1
		s = str(counter) + '.\t'
		s += gene_name + '\t' + str(degree) + '\t+' + str(num_pos) + '/-' + str(num_neg)
		print s
		if counter >= TOP_HUBS_TO_DISPLAY:
			break

def write_plot_neighbors_vs_edge_weights(num_pos_neighbors, num_neg_neighbors, edge_weight_thresholds, gene_name):
	line_pos, = plt.plot(edge_weight_thresholds, num_pos_neighbors, 'r', label = 'kataegis positive')
	line_neg, = plt.plot(edge_weight_thresholds, num_neg_neighbors, 'b', label = 'kataegis negative')
	plt.xlabel('edge weight threshold')
	plt.ylabel('number of neighbors')
	plt.title('number of neighbors with varying weight for gene ' + gene_name)
	plt.legend(handles = [line_pos, line_neg])
	plt.savefig('neighbors_' + gene_name.replace('|', '_') + '.png')
	plt.clf()

def get_gene_neighbor_dic(pairs, threshes):
	gene_neighbor_dic = {}
	for gene_name in GENE_NAMES:
		gene_neighbor_dic[gene_name] = []
	for thresh in threshes:
		graph = get_graph(pairs, thresh)
		for gene_name in GENE_NAMES:
			gene_neighbor_dic[gene_name].append(len(graph.edges(gene_name)))
	return gene_neighbor_dic

# input: vals (list of float)
def count_pos_neg(vals):
	num_pos, num_neg = 0, 0
	for val in vals:
		if val < 0:
			num_neg += 1
		else:
			num_pos += 1
	return num_pos, num_neg

# writes the ego graph of gene with gene_name to a file
def write_ego_graph(graph, gene_name):
	hub_ego = nx.ego_graph(graph, gene_name)
	pos = nx.spring_layout(hub_ego)
	nx.draw(hub_ego, pos, node_color='b', node_size=50, with_labels = True)
	nx.draw_networkx_nodes(hub_ego, pos, nodelist = [gene_name], node_size = 300, node_color='r')
	plt.savefig('ego_' + gene_name.replace('|', '_') + '.png')
	plt.clf()

# input: pairs (np.array) [num_pairs, 3] (string)
#        edge_thresh (float) values where if absolute value of edge weight is below we do not include edge
def get_graph(pairs, edge_thresh):
	G = nx.Graph()
	genes = np.unique(pairs[:, 0])
	G.add_nodes_from(genes)
	for geneA, geneB, edge_val in pairs:
		try:
			edge_val = float(edge_val)
		except:
			continue
		if geneA == geneB or abs(edge_val) < edge_thresh:
			continue
		G.add_edge(geneA, geneB, {'weight': edge_val})
	return G


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'barcode_and_q_values.py', description = "creates a tab separated value sheet with TCGA barcodes and q-value enrichment")
	parser.add_argument('input_file', help = 'input .pairs gene network file file', type = lambda x: is_valid_file(parser, x))
	parser.add_argument('-t', '--edge_threshold', type = lambda x: float_between(parser, x, 0.0, 1.0), required = True, help = 'cutoff where we will not display any edges with absolute weight lower than threshold. must be between 0.0 and 1.0 non inclusive')
	parser.add_argument('-n', '--kat_neg', type = lambda x: is_valid_file(parser, x), help = '.pairs for kataegis negative. input file should then be kataegis positive')
	return vars(parser.parse_args(argv))

def float_between(parser, arg, low, high):
	err_msg = 'must be float between ' + str(low) + ' and ' + str(high) + ' non-inclusive'
	try:
		arg = float(arg)
	except:
		parser.error(err_msg)
	if not (low < arg and arg < high):
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
