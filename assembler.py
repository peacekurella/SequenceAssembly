import argparse


def read_sequences(directory):
    return []


def split_into_kmers(sequence, args):
    return []


def generate_de_bruijin_graph(kmers, args):
    return []


def traverse_graph(graph):
    return []


def align_to_reference(assembled_seq, args):
    return []


if __name__ == "__main__":
    # add the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", default=30, help="k-mer size")
    parser.add_argument("--threshold", default=2, help="min count of k-mers to consider")
    parser.add_argument("--input", default="sars-cov-2-raw-data", help="input directory")
    parser.add_argument("--reference", default=None, help="reference genome")

    args = parser.parse_args()

    # each file is a seperate sequence
    sequences = read_sequences(args.input)

    for sequence in sequences:
        kmers = split_into_kmers(sequence, args)
        graph = generate_de_bruijin_graph(kmers, args)  # save the graph for plotting
        assembled_seq = traverse_graph(graph)
        mismatches = align_to_reference(assembled_seq, args)  # save the mismatches for plotting
