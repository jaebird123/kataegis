#     file: shared_network.py
#   author: Jesse Eaton and Jacob West-Roberts
#  created: May 9, 2017
# modified: May 10, 2017
#  purpose: Finds shared edges between inferred and known regulatory network. Plots
#             ego graph with shared edges for specified gene


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import numpy as np              # for manipulating matricies
import networkx as nx           # for creating a graph
import matplotlib.pyplot as plt # for plotting and saving plots


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

EGO_GRAPH_GENE = 'SPI1' # gene name for gene to plot on ego graph
# EGO_GRAPH_GENE = 'SLA2'


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)

	inf_gene_pairs = np.genfromtxt(args['reg_net_file'], dtype = str, delimiter = '\t')
	large_known_net = np.genfromtxt(args['known_network_file'], dtype = str, delimiter = '\t')
	
	# get set of gene names that are in both inferred and known networks
	inf_gene_set = get_inf_gene_set(inf_gene_pairs)
	kno_gene_set = set(large_known_net[:, (0, 2)].flat) # 0th an 2nd col are gene names
	gene_set = inf_gene_set.intersection(kno_gene_set)

	inf_graph = get_inf_graph(inf_gene_pairs, gene_set, 0.776)
	kno_graph = get_kno_graph(large_known_net, gene_set)

	# print number of nodes and edges in inferred and known regulatory networks
	print str(len(inf_graph.nodes())) + ' nodes in inferred graph'
	print str(len(kno_graph.nodes())) + ' nodes in known graph'
	print str(len(inf_graph.edges())) + ' edges in inferred graph'
	print str(len(kno_graph.edges())) + ' edges in known graph'

	# find edges exclusive to each graph and shared between graphs
	inf_graph_only = nx.difference(inf_graph, kno_graph)
	kno_graph_only = nx.difference(kno_graph, inf_graph)
	shared_graph   = nx.difference(inf_graph, inf_graph_only)

	# print edges shared between graphs
	print str(len(inf_graph_only.edges())) + ' of ' + str(len(inf_graph.edges())) + ' edges unique to inf_graph'
	print str(len(kno_graph_only.edges())) + ' of ' + str(len(kno_graph.edges())) + ' edges unique to kno_graph'
	print 'the ' + str(len(inf_graph.edges()) - len(inf_graph_only.edges())) + ' edge(s) in common are:'
	print shared_graph.edges()
	
	# plot ego graph for specific gene
	write_ego_graph(shared_graph, EGO_GRAPH_GENE)

#  input: known_net (np.array) [num_known_genes, 4]
#         gene_set (set of string) gene names. to be used as nodes
# output: G (nx.Graph) edges for each 
def get_kno_graph(known_net, gene_set):
	G = nx.Graph()
	G.add_nodes_from(gene_set)
	for geneA, _, geneB, _ in known_net:
		if geneA in gene_set and geneB in gene_set:
			G.add_edge(geneA, geneB)
	return G.to_undirected()

#  input: pairs (np.array) [num_pairs, 3] (string)
#         gene_set (set of string) gene names for genes to keep as nodes in graph
#         edge_thresh (float) values where if absolute value of edge weight is below we do not include edge
# output: G (nx.Graph) edges between each node pair in pairs if pairs p^2 value >= edge_threshold
def get_inf_graph(pairs, gene_set, edge_thresh):
	G = nx.Graph()
	G.add_nodes_from(gene_set)
	for geneA, geneB, edge_val in pairs:
		geneA = geneA.split('|')[0]
		geneB = geneB.split('|')[0]
		try:
			edge_val = float(edge_val)
		except:
			continue
		if geneA == geneB or abs(edge_val) < edge_thresh:
			continue
		if geneA in gene_set and geneB in gene_set:
			G.add_edge(geneA, geneB)
	return G.to_undirected()

#  input: pairs (np.array) [num_inferred_genes^2, 3] gene1, gene2, p^value. for all pairs
# output: gene_set (set of string) names of genes in pairs
def get_inf_gene_set(pairs):
	gene_list = []
	for gene1_name, gene2_name, _ in pairs:
		gene_list.append(gene1_name.split('|')[0])
		gene_list.append(gene2_name.split('|')[0])
	return set(gene_list)


# # # # # # # # #
#   P L O T S   #
# # # # # # # # #

# writes the ego graph of gene with gene_name to a file
def write_ego_graph(graph, gene_name):
	hub_ego = nx.ego_graph(graph, gene_name)
	pos = nx.spring_layout(hub_ego)
	nx.draw(hub_ego, pos, alpha = 0.5, node_color='b', node_size=1000, with_labels = True)
	nx.draw_networkx_nodes(hub_ego, pos, alpha = 0.5, nodelist = [gene_name], node_size = 1000, node_color='r', width = 5.0)
	plt.savefig('ego_' + gene_name.replace('|', '_') + '.png')
	plt.clf()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'barcode_and_q_values.py', description = "creates a tab separated value sheet with TCGA barcodes and q-value enrichment")
	parser.add_argument('-r', '--reg_net_file', type = lambda x: is_valid_file(parser, x), required = True, help = 'file containing all pairs of genes from RNA sequence data with corresponding p^2 values')
	parser.add_argument('-k', '--known_network_file', type = lambda x: is_valid_file(parser, x), required = True, help = 'file containing edges between genes in known regulatory network')
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
