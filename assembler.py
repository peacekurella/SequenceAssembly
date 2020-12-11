import os
import argparse
from Bio import SeqIO


def read_sequences(directory):
    """
    Reads in the sequences from the fastq files
    :param directory:
    :return:
    """
    files = os.listdir(directory)
    sequences = []

    for file in files:
        sequence = ""
        filename = os.path.join(directory, file)
        for record in SeqIO.parse(filename, "fastq"):
            sequence = sequence + ''.join(list(record.seq))
        sequences.append(sequence)
    return sequences


def split_into_kmers(sequence, k, threshold):
    kmers = {}

    for i in range(0, len(sequence), k):
        kmer = sequence[i: i + k + 1]
        if len(kmer) > k:
            if kmer in kmers.keys():
                kmers[kmer] = kmers[kmer] + 1
            else:
                kmers[kmer] = 1

    for key, count in kmers:
        if count <= threshold:
            kmers.pop(key)

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
    
    # reverse mapping to figure out start and stop nodes
    degree_node = {v:k for k, v in node_degree.items() if v != 0}
    
    start = None
    
    # reverse mapping would be empty if there is not eulerian path, but only cycles
    if len(degree_node) > 0:
        start_stop = sorted(degree_node, reverse=True)
        start = degree_node[start_stop[0]]
        
    # Return the nodes, edges and the starting nodes in the graph.
    return nodes, edges, start

def make_node_edge_map(edges):

    # Make a map of starting nodes to the adjacency list of that node
    node_edge_map = {}

    # Go through all edges
    for e in edges:
        n = e[0]
        # If start node exists, add destination node to adjacency list
        if n in node_edge_map:
            node_edge_map[n].add(e[1])
        # Add start node to map and initialize the adjacency list with the destination node
        else:
            adjacency_set = set() # using a set to support faster udpates
            adjacency_set.add(e[1])
            node_edge_map[n] = adjacency_set

    return node_edge_map

def traverse_graph(graph, start):
    # if there is no explicit start node
    if len(start) == 0:
        # pick any node as the start
        start = list(graph.keys())[0]
    else:
        start = start[0]

    # maintain a stack to store the nodes to visit
    path = [start]

    # accumulates the eulerian path
    eulerian_path = []

    # while stack is non-empty
    while path:
        # pick up the topmost node in the stack
        curr_node = str(path[-1])

        # if the current node is a key and has entries in its ajacency list
        if curr_node in graph and graph[curr_node]:

            # get the adjacency list
            adj_nodes = graph[curr_node]
            next_node = None

            # if only one entry, proceed to that node
            if len(adj_nodes) == 1:
                next_node = adj_nodes.pop()

            # if not, proceed to next node that would allow us to traverse
            # the rest of the graph
            else:
                # iterate through all adj_nodes
                for node in adj_nodes:
                    # if the node leads to other nodes, set that as next
                    if node in graph.keys():
                        next_node = node
                        break
                # clear the edge from the current node
                graph[curr_node].remove(next_node)
            # update the path we want to explore
            path.append(next_node)
        else:
            # if we can go any further, append the node to the path
            eulerian_path.append(path.pop())
    # reverse the accumulator as the path was populated in reverse
    eulerian_path.reverse()

    return eulerian_path

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
        tour = tour[:append_at+1] + trail + tour[append_at+1:]
        
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
            print(graph)
            break
            
    return tour, all_trails

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
    sequences = read_sequences(args.input)

    for sequence in sequences:
        kmers = split_into_kmers(sequence, args.k, args.threshold)
        k1mers = split_into_kmers(sequence, args.k - 1, 0)
        nodes, edges, start = generate_de_bruijin_graph_alt(kmers)
        # start would be None if there is not eulerian path, but only cycles
        if start == None:
            nodes = list(nodes)
            start = nodes[0]
        graph = make_node_edge_map(edges)
        tour, all_trails = traverse_graph(graph, start)
    
    #     graph = generate_de_bruijin_graph(kmers, k1mers, args)  # save the graph for plotting
    #     assembled_seq = traverse_graph(graph)
    #     mismatches = align_to_reference(assembled_seq, args)  # save the mismatches for plotting
