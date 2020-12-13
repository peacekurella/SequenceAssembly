import os
import argparse
from Bio import SeqIO
import networkx as nx


def read_sequences(directory, k):
    """
    Reads in the sequences from the fastq files
    :param directory:
    :return:
    """
    files = os.listdir(directory)
    kmer_lists = []

    for file in files:
        kmers = []
        filename = os.path.join(directory, file)
        for record in SeqIO.parse(filename, "fastq"):
            sequence = ''.join(list(record.seq))
            kmers.extend(split_into_kmers(sequence, k))
        kmer_lists.append(kmers)
    return kmer_lists


def split_into_kmers(sequence, k):
    """
    Return a list of all kmers in the sequence
    :param sequence: sequence string
    :param k: size of k mers
    :return: a list of all kmers in the sequences
    """
    kmers = []
    for i in range(len(sequence)):
        if len(sequence[i: i + k]) < k:
            break
        else:
            kmers.append(sequence[i: i + k])

    return kmers

def filter_kmers(kmers, threshold) :

    count_map = {}
    for kmer in kmers :
        if kmer in count_map :
            count_map[kmer] += 1
        else :
            count_map[kmer] = 1

    return_list = []
    for kmer in kmers :
        if count_map[kmer] >= threshold :
            return_list.append(kmer)

    return return_list

def generate_de_bruijin_graph(kmers):
    # Set of all nodes in the DB Graph
    nodes = set()
    # Set of nodes having in-degrees
    not_starts = set()
    # List of all directed edges in the graph
    edges = []

    for kmer in kmers:
        # From k-mers, get k-1mers
        k1mer1 = kmer[:-1]
        k1mer2 = kmer[1:]
        # Add k-1 mers to the nodes list
        nodes.add(k1mer1)
        nodes.add(k1mer2)

        # Add an edge between the two k-1mers
        edges.append((k1mer1, k1mer2))

        # Add destination node to the set of nodes having in degrees
        not_starts.add(k1mer2)

    start_nodes = list(nodes - not_starts)

    # Return the nodes, edges and the starting nodes in the graph.
    return nodes, edges, start_nodes


def make_node_edge_map(edges):
    # Make a map of starting nodes to the adjacency list of that node
    node_edge_map = {}

    # Go through all edges
    for e in edges:
        n = e[0]
        # If start node exists, add destination node to adjacency list
        if n in node_edge_map:
            node_edge_map[n].append(e[1])
        # Add start node to map and initialize the adjacency list with the destination node
        else:
            node_edge_map[n] = [e[1]]
    return node_edge_map


def traverse_graph(graph):
    return []


def align_to_reference(assembled_seq, args):
    return []

def generate_de_bruijin_graph_alt(kmers):

    # to keep track of degree of nodes
    node_degree = {}

    # Set of all nodes in the DB Graph
    nodes = set()

    # Set of nodes having in-degrees
    not_starts = set()

    # List of all directed edges in the graph
    edges = []

    for kmer in kmers:

        # From k-mers, get k-1mers
        prefix = kmer[:-1]
        suffix = kmer[1:]

        # Add k-1 mers to the nodes list
        nodes.add(prefix)
        nodes.add(suffix)

        # Add an edge between the two k-1mers
        edges.append((prefix, suffix))

        # fetching and updating degrees
        p_count = node_degree.setdefault(prefix, 0)
        node_degree[prefix] = p_count+1
        s_count = node_degree.setdefault(suffix, 0)
        node_degree[suffix] = s_count-1

    degree_node = {v:k for k, v in node_degree.items() if v != 0}

    start = None

    if len(degree_node) > 0:
        start_stop = sorted(degree_node, reverse=True)
        start = degree_node[start_stop[0]]

    # Return the nodes, edges and the starting nodes in the graph.
    return nodes, edges, start


def DB_graph(nodes, edges) :
    graph = nx.MultiDiGraph()

    for node in nodes :
        graph.add_node(node)

    graph.add_edges_from(edges)

    return graph


if __name__ == "__main__":
    # add the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", default=30, help="k-mer size")
    parser.add_argument("--threshold", default=2, help="min count of k-mers to consider")
    parser.add_argument("--input", default="sars-cov-2-trimmed/SE-SW-4-15", help="input directory")
    parser.add_argument("--reference", default=None, help="reference genome")

    args = parser.parse_args()

    # each file is a seperate sequence
    kmer_lists = read_sequences(args.input, args.k)
    pruned_list = filter_kmers(kmer_lists, args.threshold)

    nodes, edges, start = generate_de_bruijin_graph_alt(kmer_lists[0])

    graph = DB_graph(nodes, edges)

    #print(nx.number_weakly_connected_components(graph))
    weak_components = nx.weakly_connected_components(graph)
    counter=0
    for c in weak_components:
        subgraph = graph.subgraph(list(c))
        if nx.has_eulerian_path(subgraph) == True:
            counter+=1
    print(counter)
