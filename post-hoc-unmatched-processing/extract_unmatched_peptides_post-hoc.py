"""
Extracts peptides post-hoc peptide fusion scan.

Examples:
python3 extract_unmatched_peptides_post-hoc.py -u --experiment_name N1_Normal --all_denovo_file data/N1_Normal_all_de_novo.csv --db_matched_file data/N1_Normal_DB-matched_peptides.csv --results_file data/N1_Normal_CONTROL_full_proteome_scan.xls --reference_proteome_file ../data/SPHu.fasta
python3 extract_unmatched_peptides_post-hoc.py -u --experiment_name C6_chRCC --all_denovo_file data/C6_chRCC_all_de_novo.csv --db_matched_file data/C6_chRCC_DB-matched_peptides.csv --results_file data/C6_chRCC_CONTROL_full_proteome_scan.xls --reference_proteome_file ../data/SPHu.fasta

2-3h per 100k peptides.
"""

import argparse
import itertools
import os
import pandas
import pathlib
import re
import tqdm

from Bio import SeqIO


def peptide_in_reference_proteome(ref_db, peptide):
    """
    Method to run matching to proteome for unmatched_prior_to_scan - return and store if no match.
    """
    for protein_name, sequence in ref_db.items():
        sequence = sequence.strip('X*')
        if 'X' in sequence:
            splitter = 'X'
        else:
            splitter = '*'
        for ORF in sequence.split(splitter):
            for _ in re.finditer(peptide, ORF):
                return True
    return False


def _create_combinations_with_isoleucines(peptide):
    """
    Performs all possible combinations of string reassignment of L to I.
    """
    count = peptide.count('L')
    product = [''.join(seq) for seq in itertools.product("01", repeat=count)]
    positions = [pos for pos, char in enumerate(peptide) if char == 'L']
    peptides = []
    for combo in product:
        partitioned_peptide = list(peptide)
        for pos, item in enumerate(positions):
            if combo[pos] == '1':
                partitioned_peptide[item] = 'I'
        peptides.append(''.join(partitioned_peptide))
    return peptides


def matched_genome(ref_db, peptide):
    for il_variant in _create_combinations_with_isoleucines(peptide):
        if peptide_in_reference_proteome(ref_db, il_variant):
            return False
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='python3 main.py')
    parser.add_argument('--experiment_name', type=str, help='name of the experiment; also sets the folder name for results', required=True)
    parser.add_argument('--all_denovo_file', type=str, help='all de novo data file from Peaks', required=True)
    parser.add_argument('--db_matched_file', type=str, help='DB-matched data file from Peaks', required=True)
    parser.add_argument('--results_file', type=str, help='Fused pepetide scan results file', required=True)
    parser.add_argument('--reference_proteome_file', type=str, help='FASTA file to use as the reference proteome', required=True)
    parameters = parser.parse_args()

    results_directory = 'results'
    results_path = results_directory + os.sep + parameters.experiment_name
    pathlib.Path(results_path).mkdir(parents=True, exist_ok=True)

    # get processed data
    df_processed = pandas.read_excel(parameters.results_file)

    # process raw data
    df_db_matched = pandas.read_csv(parameters.db_matched_file, dtype={'Scan': int})
    df_all_denovo = pandas.read_csv(parameters.all_denovo_file, dtype={'ALC (\%)': int, 'Scan': int})
    df_db_matched['Peptide'] = df_db_matched['Peptide'].apply(lambda x: re.sub(r'\([^()]*\)', '', x))
    df_all_denovo['Peptide'] = df_all_denovo['Peptide'].apply(lambda x: re.sub(r'\([^()]*\)', '', x))

    ref_db = {seq.name: str(seq.seq) for seq in SeqIO.parse(parameters.reference_proteome_file, 'fasta')}

    # selection of de novo (oroginal HLA) peptides for fusion candidate scan
    scan_ids_not_db_matched = set(df_all_denovo['Scan'].values) - set(df_db_matched['Scan'].values)
    df_all_denovo_unmatched = df_all_denovo[df_all_denovo['Scan'].isin(scan_ids_not_db_matched)]

    unmatched = set(df_all_denovo_unmatched.Scan.values) - set(df_processed.ScanID.values)
    unmatched_prior_to_scan = df_all_denovo_unmatched[df_all_denovo_unmatched['Scan'].isin(unmatched)]
    unmatched_prior_to_scan = unmatched_prior_to_scan[['Peptide', 'Scan', 'ALC (%)']].drop_duplicates()

    unmatched = []
    for values in tqdm.tqdm([data[1] for data in unmatched_prior_to_scan.iterrows()]):
        peptide = values['Peptide']
        unmatched.append(matched_genome(ref_db, peptide))

    unmatched_prior_to_scan[unmatched].to_csv(results_path + os.sep + 'Unmatched_denovo_peptides.csv', index=False)
