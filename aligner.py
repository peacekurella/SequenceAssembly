from Bio import Align
from Bio import SeqIO
import numpy as np
import json

with open('contigs1.txt') as contigs:
    contigs = contigs.readlines()

for record in SeqIO.parse("reference/reference.fna", "fasta"):
    reference = ''.join(list(record.seq))


aligner = Align.PairwiseAligner()
scores = []
for contig in contigs:
    contig = contig.strip('\n')
    alignment = aligner.align(reference, contig)[0]
    print(alignment.score)
    scores.append(alignment.score)

scores = np.array(scores)
best_contig = np.argmax(scores)

best_contig = contigs[best_contig].strip('\n')
alignment = aligner.align(reference, best_contig)[0]

with open('alignment_plot.json', 'w') as op:
    json.dump({
        'align': alignment.aligned
    }, op)