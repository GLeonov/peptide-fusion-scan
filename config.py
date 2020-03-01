import argparse
import os
import pathlib

parser = argparse.ArgumentParser(prog='python3 main.py')
parser.add_argument('--experiment_name', type=str, help='name of the experiment; also sets the folder name for results', required=True)
parser.add_argument('--all_denovo_file', type=str, help='all de novo data file from Peaks', required=True)
parser.add_argument('--db_matched_file', type=str, help='DB-matched data file from Peaks', required=True)
parser.add_argument('--reference_proteome_file', type=str, help='FASTA file to use as the reference proteome', required=True)
parser.add_argument('--min_peptide_length', type=int, help='smallest peptide length to scan', required=True)
parser.add_argument('--max_peptide_length', type=int, help='largest peptide length to scan', required=True)
parser.add_argument('--alc_score_cutoff', type=int, help='ALC score cutoff', required=True)
parser.add_argument('--max_fusion_distance', type=int, help='maximum distance to scan for peptide fusions', required=True)
parser.add_argument('--include_PTMs', type=bool, help='if True, converts C(+57.02)FKHSGTGM(+15.99)VHR to CFKHSGTGMVHR; else, skips the row', required=True)
parameters = parser.parse_args()

# create working directories if it does not exist; all output will be stored here
results_directory = 'results'
results_path = results_directory + os.sep + parameters.experiment_name
pathlib.Path(results_path).mkdir(parents=True, exist_ok=True)
temp_files_path = results_directory + os.sep + parameters.experiment_name + os.sep + 'multithreaded_temp_files' + os.sep
