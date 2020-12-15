# Sequence Assembly Using deBruijn Graphs

## Description
De-Bruijn graph traversal is a popular method to assemble contigs
from genomes. De-Bruijn graphs are constructed by splitting the 
genome into K-mers and using them as the edges of the graph while 
K-1 mers become the nodes. This essentially turns the assembly
problem into a graph traversal problem. While De-Bruijn graphs 
have advantages, like offering compact representation of genomes,
they also have drawbacks. The major drawback is that the technique
assumes perfect sequencing, which is quite difficult to achieve in
practice. Other drawbacks include repeat edges which lead to bloated
graphs and the existence of multiple assemblies which can be difficult
to resolve. In this project, we explored how De-Bruijn graphs can be used
for assembling the SARS-Cov2 genome, while addressing some of the 
shortcomings and taking note of research that provides possible 
solutions to the others. We also give a description of the steps and 
methodology used in this project, like the data cleaning techniques, 
which will then enable other researchers to replicate our results or 
expand on our work. We also perform experiments with a range of 
parameters as a way to conduct analysis and provide neat and interactive 
visualizations of our results. 

## Visualization Demo
<img src="https://github.com/peacekurella/SequenceAssembly/blob/main/readme_res/demo.gif" width="800" height="600" />

