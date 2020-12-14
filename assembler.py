import os
import json
import argparse
from Bio import SeqIO
from Bio import Align
import numpy as np
import networkx as nx
import pandas as pd


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


def filter_kmers(kmers, threshold):
    count_map = {}

    for kmer in kmers:
        if kmer in count_map:
            count_map[kmer] += 1
        else:
            count_map[kmer] = 1

    return_list = []
    for kmer in kmers:
        if count_map[kmer] >= threshold:
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
        node_degree[prefix] = p_count + 1
        s_count = node_degree.setdefault(suffix, 0)
        node_degree[suffix] = s_count - 1

    degree_node = {v: k for k, v in node_degree.items() if v != 0}

    start = None

    if len(degree_node) > 0:
        start_stop = sorted(degree_node, reverse=True)
        start = degree_node[start_stop[0]]

    # Return the nodes, edges and the starting nodes in the graph.
    return nodes, edges, start


def get_start_nodes(edges):
    node_degree = {}
    for edge in edges:
        prefix = edge[0]
        suffix = edge[1]
        p_count = node_degree.setdefault(prefix, 0)
        node_degree[prefix] = p_count + 1
        s_count = node_degree.setdefault(suffix, 0)
        node_degree[suffix] = s_count - 1

    degree_node = {v: k for k, v in node_degree.items() if v != 0}

    start = None

    if len(degree_node) > 0:
        start_stop = sorted(degree_node, reverse=True)
        start = degree_node[start_stop[0]]

    # Return the nodes, edges and the starting nodes in the graph.
    return start


def traverse_graph_alt(graph, start):
    # to collect all eulerian paths/cycles in a graph
    all_trails = list()

    # the eularian tour of the entire graph
    tour = []
    tour.append(start)

    skip_trail = True

    while (True):
        # start an eulerian trail
        trail = []

        curr_node = start

        # traverse a trail until we can't go further
        while (True):

            # terminate if can't go further
            if curr_node not in graph:
                break

            # pick a next node
            next_node = graph[curr_node].pop()

            # if the adjacency list becomes emtpy for the current node, delete
            if len(graph[curr_node]) == 0:
                del graph[curr_node]

            # append next node to trail
            trail.append(next_node)

            # if we circle back to start we have covered the trail
            if next_node == start:
                break;

            # if not move on to next node
            curr_node = next_node

        # we skip adding the first trail as it would reflect in the tour
        if skip_trail == False:
            # after finishing a trail, add it to all tours
            all_trails.append(list(trail))

        skip_trail = False

        # where to append the trail in the tour
        append_at = tour.index(start)

        # introducing the trail inbetween the tour
        tour = tour[:append_at + 1] + trail + tour[append_at + 1:]

        # done if no more nodes to explore
        if len(graph) == 0:
            break

        new_start_possible = False

        for node in tour:
            if node in graph:
                start = node
                new_start_possible = True
                break

        if not new_start_possible:
            print("error, tour exploration terminated with remaining graph:")
            # print(graph)
            break

    return tour, all_trails


def DB_graph(nodes, edges):
    graph = nx.MultiDiGraph()

    for node in nodes:
        graph.add_node(node)

    graph.add_edges_from(edges)

    return graph


def get_contig_from_path(path):
    contig = ''
    if path is None:
        return
    for p in path:
        if p is None:
            continue
        if contig == '':
            contig += p
        else:
            contig += p[-1]
    return contig


def get_alignment_scores(contigs, reference_file):
    # read the reference genome
    for record in SeqIO.parse(reference_file, "fasta"):
        reference = ''.join(list(record.seq))

    aligner = Align.PairwiseAligner()
    scores = []
    alignments = []
    for contig in contigs:
        contig = contig.strip('\n')
        if len(contig) > 0:
            alignment = aligner.align(reference, contig)[0]
            scores.append(alignment.score)
            alignments.append(alignment.aligned)
        else:
            scores.append(0)
            alignments.append([])

    scores = np.array(scores)
    best_contig = np.argmax(scores)

    return scores, scores[best_contig], alignments[best_contig], best_contig


def get_lens(contigs):
    num_contigs = len(contigs)

    lens = []
    for contig in contigs:
        lens.append(len(contig))

    return num_contigs, lens


if __name__ == "__main__":
    # add the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", default=30, help="k-mer size")
    parser.add_argument("--threshold", default=2, help="min count of k-mers to consider")
    parser.add_argument("--input", default="sars-cov-2-trimmed/PE-SW-4-15", help="input directory")
    parser.add_argument("--output", default="output", help="output folder")
    parser.add_argument("--ref", default="reference/reference.fna", help="reference genome file")

    args = parser.parse_args()

    results = pd.DataFrame({
        'k': [],
        'threshold': [],
        'al_scores': [],
        'lens': [],
        'n_cons': [],
        'best_al_score': [],
        'read': []
    })

    for args.threshold in range(2, 5):
        for args.k in range(30, 100, 5):
            # each file is a seperate sequence
            kmer_lists = read_sequences(args.input, args.k)

            for read_id in [0, 1]:

                print('Processing read:', read_id)

                output_directory = os.path.join(args.output, 'k_{}_th_{}'.format(args.k, args.threshold))
                if not os.path.isdir(output_directory):
                    os.makedirs(output_directory)

                pruned_list = filter_kmers(kmer_lists[read_id], args.threshold)

                nodes, edges, start = generate_de_bruijin_graph_alt(pruned_list)

                graph = DB_graph(nodes, edges)

                weak_components = nx.weakly_connected_components(graph)
                contigs = []
                paths = []
                for c in weak_components:
                    subgraph = graph.subgraph(list(c))
                    c_edges = subgraph.edges
                    c_map = make_node_edge_map(c_edges)
                    start_node = get_start_nodes(c_edges)
                    path, trail = traverse_graph_alt(c_map, start_node)
                    contig = get_contig_from_path(path)
                    if contig is not None:
                        paths.append(path)
                        contigs.append(contig)

                output_file = os.path.join(output_directory, 'read_{}.txt'.format(read_id))

                with open(output_file, 'w') as op:
                    op.writelines("%s\n" % contig for contig in contigs)

                scores, best_score, best_alignment, best_contig_idx = get_alignment_scores(contigs, args.ref)

                path_file = os.path.join(output_directory, "path_file_{}.json".format(read_id))
                with open(path_file, 'w') as op:
                    json.dump({
                        'path': paths[best_contig_idx]
                    }, op)

                best_alignment_file = os.path.join(output_directory, 'best_alignment_{}.json'.format(read_id))
                with open(best_alignment_file, 'w') as op:
                    json.dump({
                        'align': best_alignment
                    }, op)

                num_contigs, lens = get_lens(contigs)

                results = results.append({
                    'k': args.k,
                    'threshold': args.k,
                    'al_scores': scores,
                    'lens': lens,
                    'n_cons': num_contigs,
                    'best_al_score': best_score,
                    'read': read_id
                }, ignore_index=True)

                results.to_csv(os.path.join(args.output, 'results.csv'))
            break
        break
