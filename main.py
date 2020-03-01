"""
Peptide fusion scanning application that takes parameters passed from a bash script.
This allows to setup multiple runs to be run in a sequence.
"""

# Python 3.7 built-in
import datetime
import multiprocessing
import os
import pandas
import pathlib
import re
import shutil
import subprocess
import time

# additional dependencies
from Bio import SeqIO  # pip3 install biopython
import tqdm  # pip3 install tqdm

# local application
import config
import fused_peptides


def running_python_jobs_count():
    return len([i for i in subprocess.check_output(['ps', 'uax']).splitlines() if 'python' in str(i) and 'main.py' in str(i)])


def get_input_data():
    print('\n%s || %s || Reading in data' % (datetime.datetime.now(), config.parameters.experiment_name))
    # read data
    df_db_matched = pandas.read_csv(config.parameters.db_matched_file, dtype={'Scan': int})
    df_all_denovo = pandas.read_csv(config.parameters.all_denovo_file, dtype={'ALC (\%)': int, 'Scan': int})
    ref_db = {seq.name: str(seq.seq) for seq in SeqIO.parse(config.parameters.reference_proteome_file, 'fasta')}

    # process PTMs
    if config.parameters.include_PTMs:
        df_db_matched['Peptide'] = df_db_matched['Peptide'].apply(lambda x: re.sub(r'\([^()]*\)', '', x))
        df_all_denovo['Peptide'] = df_all_denovo['Peptide'].apply(lambda x: re.sub(r'\([^()]*\)', '', x))
    else:
        df_db_matched = df_db_matched[~df_db_matched['Peptide'].str.contains('\+')]
        df_all_denovo = df_all_denovo[~df_all_denovo['Peptide'].str.contains('\+')]

    print('%s || %s || Storing initial peptide files' % (datetime.datetime.now(), config.parameters.experiment_name))
    # save DB-matched file with razor peptides removed
    df_db_matched_no_razor_peptides = df_db_matched[df_db_matched['Accession'].apply(lambda x: str(x).count(':')) == 0]
    df_db_matched_no_razor_peptides.to_csv(config.results_path + os.sep + 'DB-matched_no_razor_peptides.csv',
                                           index=False)

    # save the unique protein IDs with razor peptides removed
    # TODO: clarify which entries to keep - i.e. do we just keep the first occurence of each duplicate protein?
    df_db_matched_no_razor_peptides_unique = df_db_matched_no_razor_peptides.drop_duplicates(subset='Accession',
                                                                                             keep='first')
    df_db_matched_no_razor_peptides_unique.to_csv(
        config.results_path + os.sep + 'DB-matched_no_razor_peptides_unique_proteinIDs.csv', index=False)

    # get all denovo scan IDs that are not found in DB-matched data
    scan_ids_not_db_matched = set(df_all_denovo['Scan'].values) - set(df_db_matched['Scan'].values)
    df_all_denovo_unmatched = df_all_denovo[df_all_denovo['Scan'].isin(scan_ids_not_db_matched)]
    df_all_denovo_unmatched.to_csv(config.results_path + os.sep + 'Fused_peptide_candidates.csv', index=False)
    df_all_denovo_unmatched = df_all_denovo_unmatched[['Peptide', 'Scan', 'ALC (%)']].drop_duplicates()
    print('%s || %s || Reading in data finished' % (datetime.datetime.now(), config.parameters.experiment_name))
    return df_all_denovo_unmatched, ref_db


def run(fused_peptide_candidates, ref_db):
    print('%s || %s || Setting up runs' % (datetime.datetime.now(), config.parameters.experiment_name))
    fused_peptide_candidates = [fused_peptides.Peptide(*data[1].values) for data in fused_peptide_candidates.iterrows()]
    number_of_cores = multiprocessing.cpu_count()

    # perform peptide fusion; if recent files exist, prevent override unless multithreaded_temp_files are deleted
    pathlib.Path(config.temp_files_path).mkdir(parents=True, exist_ok=True)

    # setup jobs
    jobs = []
    for n_peptide, peptide in enumerate(fused_peptide_candidates[:10]):
        job = multiprocessing.Process(target=fused_peptides.find_peptide_fusions,
                                      args=(config.parameters.experiment_name,
                                            peptide,
                                            n_peptide,
                                            config.temp_files_path,
                                            ref_db,
                                            config.parameters.max_fusion_distance))
        jobs.append(job)

    # run jobs
    print('%s || %s || Scanning for peptide fusions' % (datetime.datetime.now(), config.parameters.experiment_name))
    progress_bar = tqdm.tqdm(total=len(jobs))
    while len(jobs) > 0:
        if running_python_jobs_count() <= number_of_cores:
            job = jobs.pop()
            job.start()
            progress_bar.update(1)
    while running_python_jobs_count() > 1:  # wait for jobs to complete
        time.sleep(1)
    progress_bar.close()
    print('%s || %s || Scanning completed' % (datetime.datetime.now(), config.parameters.experiment_name))


def output_final_results():
    print('%s || %s || Aggregating results files' % (datetime.datetime.now(), config.parameters.experiment_name))
    # assemble all results into a single file
    header = ['Peptide', 'Length', 'ProteinID',
              'Fragment 1', 'Match 1 Start', 'Match 1 Stop', 'Cis Distance',
              'Fragment 2', 'Match 2 Start', 'Match 2 Stop',
              'Scan', 'ALCScore', 'Origin', 'ParentID']

    df = pandas.DataFrame(columns=header)
    for f in os.listdir(config.temp_files_path):
        df = pandas.concat([df, pandas.read_csv(config.temp_files_path + os.sep + f, names=header)], axis=0)

    try:
        shutil.rmtree(config.temp_files_path)
    except OSError:
        print('Directory already removed or no permissions to delete.')

    df = df.fillna('N/A')
    df['Direction'] = df.apply(lambda row: fused_peptides.get_direction(row['Match 1 Start'], row['Match 2 Start'],
                                                                        row['Origin']), axis=1)
    df['Fused Peptide'] = df.apply(
        lambda row: fused_peptides.get_fused_fragments(row['Direction'], row['Fragment 1'], row['Fragment 2']), axis=1)

    # file outputting
    df.to_csv(config.results_path + os.sep + 'Matched_all.csv', index=False)
    df[df['Direction'] == 'Forward'].to_csv(config.results_path + os.sep +
                                            'Matched_denovo_Fused_Forward.csv', index=False)
    df[df['Direction'] == 'Reverse'].to_csv(config.results_path + os.sep +
                                            'Matched_denovo_Fused_Reverse.csv', index=False)
    df[df['Origin'] == 'de novo DB-matched'].to_csv(config.results_path + os.sep +
                                                    'Matched_denovo_DB-matched.csv', index=False)
    df[df['Origin'] == 'de novo non-matched'].to_csv(config.results_path + os.sep + 'Unmatched_denovo.csv',
                                                     index=False)
    print('%s || %s || Fusion peptide analysis completed\n' % (datetime.datetime.now(), config.parameters.experiment_name))


if __name__ == '__main__':
    run(*get_input_data())
    output_final_results()
