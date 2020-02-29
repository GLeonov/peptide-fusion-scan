# Peptide Fusion Scan App
Pre-processing of all de novo peptide options. Eliminates any peptide scan IDs
if they have been found in a databse matching the human proteome.

Using a set of experimentally determined peptides (by MS/MS), this application
partitions each peptide into sub-fragments and determines their potential
origin from a reference library of amino acid sequences (FASTA format).

Generates an Excel output file providing the origins of the subfractions of
each experimentally determined peptide, if possible.

Author:         Dr German Leonov
                University of York
                leonov89@gmail.com

Created:        29/02/2020
Last modified:  29/02/2020
