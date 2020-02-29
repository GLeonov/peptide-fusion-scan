# Peptide Fusion Scan App
Pre-processing of all de novo peptide options. Eliminates any peptide scan IDs
if they have been found in a databse matching the human proteome.

Using a set of experimentally determined peptides (by MS/MS), this application
partitions each peptide into sub-fragments and determines their potential
origin from a reference library of amino acid sequences (FASTA format).

Generates outputs segmenting by relevant groups of interest.

Author: Dr German Leonov, University of York, leonov89@gmail.com
