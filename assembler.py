import os
import argparse
from Bio import SeqIO


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
        if len(sequence[i: i + k + 1]) < k:
            break
        else:
            kmers.append(sequence[i: i + k + 1])

    return kmers


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
