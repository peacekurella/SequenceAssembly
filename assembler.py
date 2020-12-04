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


def generate_de_bruijin_graph(kmers, k1mers, args):
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
    parser.add_argument("--input", default="sars-cov-2-trimmed/SE-SW-4-15", help="input directory")
    parser.add_argument("--reference", default=None, help="reference genome")

    args = parser.parse_args()

    # each file is a seperate sequence
    sequences = read_sequences(args.input)

    for sequence in sequences:
        kmers = split_into_kmers(sequence, args.k, args.threshold)
        k1mers = split_into_kmers(sequence, args.k - 1, 0)
    #     graph = generate_de_bruijin_graph(kmers, k1mers, args)  # save the graph for plotting
    #     assembled_seq = traverse_graph(graph)
    #     mismatches = align_to_reference(assembled_seq, args)  # save the mismatches for plotting
