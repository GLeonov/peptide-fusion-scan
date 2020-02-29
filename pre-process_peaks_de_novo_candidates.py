"""
Pre-processing of all de novo peptide options. Eliminates any peptide scan IDs
if they have been found in a databse matching the human proteome. An output is
produced as an input for the fusion software.

Author:         German Leonov
                University of York
                german.leonov@york.ac.uk
                
Created:        04/07/2017
Last modified:  04/07/2017
"""

import re
import os
import xlsxwriter


class Peptide:
    def __init__(self):
        self.experiment = None
        self.name = None
        self.scanID = None
        self.score = None


def read_in_peptides(de_novo=False, matched_DB_IDs={}):
    """
    Reads in peptides info and filters.
    """
    data_files = os.listdir('data')
    
    if de_novo:
        data_files = [f for f in data_files if 'all_de_novo' in f]
        peptides = []   # list of peptides that need to be considered
    else:
        data_files = [f for f in data_files if 'DB-matched_peptides' in f]
    
    for data_file in data_files:
        print data_file
        dataset_name = data_file.split('_')[0]
        matched_DB_IDs[dataset_name] = matched_DB_IDs.get(dataset_name, [])
        
        for line in open('data/%s' % data_file).read().split('\n')[1:]:
            line = line.split(',')
            
            if len(line) > 1:
                
                if de_novo:
                    peptide = Peptide()
                    peptide.experiment = dataset_name
                    peptide.name = re.sub(r'\([^()]*\)', '', line[3])
                    peptide.scanID = int(line[4])
                    peptide.score = int(line[6])
                    
                    # discard peptides that are not 8-12mers
                    if len(peptide.name) < 7 or len(peptide.name) > 25:
                        continue
                    
                    # discard peptides that have scan IDs matching the human proteome
                    if peptide.scanID in matched_DB_IDs[dataset_name]:
                        continue
                    
                    # discard peptides with an ALC score below 50%
                    if peptide.score < 50:
                        continue
                    
                    peptides.append(peptide)
                else:
                    matched_DB_IDs[dataset_name].append(int(line[9]))
    
    if de_novo:
        return peptides
    else:
        return matched_DB_IDs


def write_output(peptides):
    """
    Writes out peptides for splicing.
    """
    f = open('data/HLA-Alleles-DeNovoSearch.csv', 'w')
    f.write('Experiment,Peptide,ScanID,ScoreALC\n')
    for peptide in peptides:
        f.write('%s,%s,%d,%d\n' % (peptide.experiment, peptide.name, peptide.scanID, peptide.score))
    f.close()

if __name__ == '__main__':
    matched_DB_IDs = read_in_peptides(de_novo=False)
    peptides = read_in_peptides(de_novo=True, matched_DB_IDs=matched_DB_IDs)
    write_output(peptides)
