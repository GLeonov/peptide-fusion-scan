"""
Using a set of experimentally determined peptides (by MS/MS), this application
partitions each peptide into sub-fragments and determines their potential
origin from a reference library of amino acid sequences (FASTA format).

Generates an Excel output file providing the origins of the subfractions of
each experimentally determined peptide, if possible.


Dependencies:   xlsxwriter  (python-xlsxwriter)

Author:         German Leonov
                University of York
                german.leonov@york.ac.uk
                
Created:        17/11/2016
Last modified:  10/10/2017
"""

import itertools
import multiprocessing
import time, datetime
import re
import os
from operator import itemgetter
import shutil
import random

import xlsxwriter  


class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


class Peptide:
    def __init__(self, name, scanID, score):
        self.name = name
        self.scanID = scanID
        self.score = score
        self.fragments = []


class FragmentPair:
    def __init__(self, pair):
        self.pair = pair
        self.matchPairs = []


class Match:
    def __init__(self, name, start, stop, number):
        self.name = name
        self.start = start
        self.stop = stop
        self.ORF = number


class FusedPeptide:
    def __init__(self):
        self.name = None
        
        self.fragment_1 = None
        self.fragment_1_start = None
        self.fragment_1_stop = None
        self.match_1_name = None
        self.match_1_start = None
        self.match_1_stop = None
        
        self.fragment_2 = None
        self.fragment_2_start = None
        self.fragment_2_stop = None
        self.match_2_name = None
        self.match_2_start = None
        self.match_2_stop = None
        
        self.cis_distance_12 = 'N/A'


def time_stamp():
    """
    Returns the date and time.
    """
    return '%s%s%s' % (color.RED, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), color.END)
    

def read_in_peptides(peptides_file, min_peptide_length, max_peptide_length, strain):
    """
    Reads in the peptides data.
    """
    print time_stamp(), 'Reading in %s ...' % peptides_file
    data = {}
    
    f = open(peptides_file)
    f.readline()        # discard header
    c = 0
    for line in f:
        c += 1
        if c % 10000 == 0:
            print c
        line = line.strip('\n').split(',')
        if len(line) > 1:
            experiment = line[0]
            name = line[1]
            scanID = int(line[2])
            score = int(line[3])
            peptide = Peptide(name, scanID, score)
            data[experiment] = data.get(experiment, []) + [peptide]
    
    f.close()
    print time_stamp(), 'Completed!\n'
    return data


def read_in_references(references_file):
    """
    Reads in the reference proteins.
    
    Returns: {protein_name -> amino_acid_sequence}
    """
    references = {}
    print time_stamp(), 'Reading in %s ...' % references_file
    f = open(references_file)
    first = True
    for line in f:
        if '>' in line[0]:
            if not first:
                references[name] = sequence
                #print [name, sequence]
            else:
                first = False
            name = line[1:].split()[0].strip('\r\n\\')
            sequence = ''
        else:
            sequence += line.strip('\r\n\\')
    print time_stamp(), 'Completed!\n'

    return references


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


def _fragment_peptide(peptide):
    """
    Fragments the provided peptide according to the length requirement.
    
    Returns a list of fragment pairs, in order, sorted by full length peptide fragment first.
    """
    fragments = []
    
    # generate I/L variants
    IL_variants = _create_combinations_with_isoleucines(peptide.name)
    for fragment in IL_variants:
        fragments.append(FragmentPair(['', fragment]))
    
    # generate spliced peptide fragments
    for IL_variant in IL_variants:
        for position in range(2, len(IL_variant)-1):
            fragments.append(FragmentPair([IL_variant[:position], IL_variant[position:]]))
    
    # sort fragments, longest fragments (full length) first - necessary for excluding them from splice searches if they match database in full
    fragments.sort(key=lambda x: len(x.pair[-1]), reverse=True)
    
    return fragments


def _search_for_fragment_in_protein_by_splitting_ORFs(peptide, references_data, max_search_distance=40):
    """
    Creates the reference search space for the sequence. Updates the Peptide
    object's Fragment objects with any Match objects it finds.
    
    Returns the updated Peptide.
    """
    full_length_match_found = False
    for fragment in peptide.fragments:
        # stop looking for matches if full length fragment match found
        if full_length_match_found:
            break
        
        for protein_name, sequence in references_data.items():
            sequence = sequence.strip('X*')
            
            distance = 0
            if 'X' in sequence:
                splitter = 'X'
            else:
                splitter = '*'
            for number, ORF in enumerate(sequence.split(splitter)):
                
                # determine the side which has the big fragment (defined in fragment_peptide method)
                if len(fragment.pair[1]) < len(fragment.pair[0]):
                    large_fragment = fragment.pair[0]
                    small_fragment = fragment.pair[1]
                else:
                    large_fragment = fragment.pair[1]
                    small_fragment = fragment.pair[0]
                
                for match_data_large in re.finditer(large_fragment, ORF):
                    
                    match_large = Match(protein_name, match_data_large.start()+distance, match_data_large.end()-1+distance, number)
                    
                    # determine the search space start for small fragment matching
                    starting_position = match_data_large.start()-(max_search_distance+len(small_fragment))
                    if starting_position < 0:
                        starting_position = 0
                    
                    # full length matches
                    if len(small_fragment) == 0:
                        full_length_match_found = True
                        fragment.matchPairs.append(['N/A', match_large])
                        continue
                    
                    for match_data_small in re.finditer(small_fragment, ORF[starting_position:match_data_large.end()+1+max_search_distance+len(small_fragment)]):
                        
                        match_small = Match(protein_name, starting_position+match_data_small.start()+distance, starting_position+match_data_small.end()-1+distance, number)
                        
                        # discard results if large fragment match overlaps with small fragment match
                        if match_small.stop >= match_large.start and match_small.start <= match_large.stop:
                            continue
                        
                        if len(fragment.pair[1]) < len(fragment.pair[0]):
                            fragment.matchPairs.append([match_large, match_small])
                        else:
                            fragment.matchPairs.append([match_small, match_large])
                    
                distance += len(ORF)+1
    return peptide


def create_fused_peptide(f, fragment_1, fragment_2, match_1, match_2, peptide):
    """
    Writes out fused peptide info.
    """
    # create fused peptide info
    fused_peptide = FusedPeptide()
    fused_peptide.name = peptide.name
    
    # fused peptide match
    if fragment_1:
        fused_peptide.fragment_1 = fragment_1
        fused_peptide.match_1_name = match_1.name
        fused_peptide.match_1_start = match_1.start
        fused_peptide.match_1_stop = match_1.stop
        
        fused_peptide.fragment_2 = fragment_2
        fused_peptide.match_2_name = match_2.name
        fused_peptide.match_2_start = match_2.start
        fused_peptide.match_2_stop = match_2.stop
        
        if match_1.stop < match_2.start:
            if match_1.stop+1 == match_2.start:
                return 0    # skip output if fusion is not sensible
            fused_peptide.cis_distance_12 = match_2.start - match_1.stop - 1
        else:
            fused_peptide.cis_distance_12 = match_1.start - match_2.stop - 1
    
    # full length peptide match
    else:
        fused_peptide.fragment_1 = fragment_2
        fused_peptide.match_1_name = match_2.name
        fused_peptide.match_1_start = match_2.start
        fused_peptide.match_1_stop = match_2.stop
        fused_peptide.cis_distance_12 = 'N/A'
        
    write_temp_file(f, fused_peptide, peptide)


def find_peptide_fusions(strain, experiment, peptide, n_peptide, directory, references_data, max_search_distance):
    """
    Finds peptides that could be fused and outputs their matches.
    """
    peptide.fragments = _fragment_peptide(peptide)
    
    # find fragments in experimental data
    peptide = _search_for_fragment_in_protein_by_splitting_ORFs(peptide, references_data, max_search_distance)
    
    f_name = '%s%s_%d.txt' % (directory, experiment, n_peptide)
    f = open(f_name, 'w')
    
    # outputting the fusions
    for fragment in peptide.fragments:
        for match_pair in fragment.matchPairs:
            create_fused_peptide(f, fragment.pair[0], fragment.pair[1], match_pair[0], match_pair[1], peptide)                             
    f.close()


def write_temp_file(f, fused_peptide, peptide):
    """
    Write a temporary text file for later.
    """
    if fused_peptide.cis_distance_12 != 'N/A':
        f.write('%s,%d,%s,%s,%d,%d,%s,%s,%d,%d,%d,%d\n' %
                (fused_peptide.name,
                len(fused_peptide.name),
                fused_peptide.match_1_name,
                fused_peptide.fragment_1,
                fused_peptide.match_1_start+1,
                fused_peptide.match_1_stop+1,
                str(fused_peptide.cis_distance_12),
                fused_peptide.fragment_2,
                fused_peptide.match_2_start+1,
                fused_peptide.match_2_stop+1,
                peptide.scanID,
                peptide.score))
    else:
        f.write('%s,%d,%s,%s,%d,%d,%s,%s,%s,%s,%d,%d\n' %
                (fused_peptide.name,
                len(fused_peptide.name),
                fused_peptide.match_1_name,
                fused_peptide.fragment_1,
                fused_peptide.match_1_start+1,
                fused_peptide.match_1_stop+1,
                'N/A',
                'N/A',
                'N/A',
                'N/A',
                peptide.scanID,
                peptide.score))


def run_jobs(jobs_to_do, max_CPUs):
    """
    Runs jobs that need to be submitted for processing.
    """
    running_jobs = {}
    total_jobs = len(jobs_to_do)

    print time_stamp(), '%sTotal jobs to do: %s%d%s' % (color.BOLD, color.BLUE,  total_jobs, color.END)
    
    while len(jobs_to_do) > 0:
        finished_jobs = []

        # identify finished jobs
        for job in running_jobs:
            alive = job.is_alive()
            if not alive:
                finished_jobs.append(job)
                

        # remove finished jobs from running jobs list
        for job in finished_jobs:
            running_jobs.pop(job)

        # try start more jobs
        jobs_to_run = max_CPUs - len(running_jobs)
        if jobs_to_run > 0:
            for job in range(jobs_to_run):
                if len(jobs_to_do) > 0:
                    new_job = jobs_to_do.pop(0)
                    new_job.start()
                    running_jobs[new_job] = 0              
    
    while len(running_jobs) > 0:
        finished_jobs = []
        for job in running_jobs:
            alive = job.is_alive()
            if not alive:
                finished_jobs.append(job)
        for job in finished_jobs:
            running_jobs.pop(job)
    
    print time_stamp(), '%s%s%d%%%s%s of jobs completed.%s' % (color.BOLD, color.BLUE, int((total_jobs - len(jobs_to_do)) / float(total_jobs) * 100), color.END, color.BOLD, color.END)


def multiprocess_peptide_search(strain, peptides_data, max_CPUs, references_data, max_search_distance):
    """
    Submits peptide search for multiprocessing
    """
    try:
        directory = os.getcwd() + os.sep + 'multithreaded_temp_files' + os.sep
        # clean up first if possible
        try:
            shutil.rmtree(directory)
        except OSError:
            pass
        os.makedirs(directory)
    except OSError:
        pass
    
    jobs_to_do = []
    for experiment, peptides in sorted(peptides_data.items(), key=itemgetter(0)):
        for n_peptide, peptide in enumerate(peptides):
            job = multiprocessing.Process(target=find_peptide_fusions, args=(strain, experiment, peptide, n_peptide, directory, references_data, max_search_distance))
            jobs_to_do.append(job)
    
    if len(jobs_to_do) > 0:
        run_jobs(jobs_to_do, max_CPUs)
    else:
        print '\nDid not find any jobs to do! Stopping the application.'


def _write_header(sheet, bold):
    """
    Writed a header in an Excel spread sheet.
    """
    sheet.autofilter('A1:P1')
    sheet.set_column(0, 25, 10)
    sheet.write(0, 0, 'ExperimentID', bold)
    sheet.write(0, 1, 'Peptide', bold)
    sheet.write(0, 2, 'Length', bold)
    sheet.write(0, 3, 'ProteinID', bold)
    sheet.write(0, 4, 'Fragment 1', bold)
    sheet.write(0, 5, 'Match 1 Start', bold)
    sheet.write(0, 6, 'Match 1 Stop', bold)
    sheet.write(0, 7, 'Cis Distance', bold)
    sheet.write(0, 8, 'Fragment 2', bold)
    sheet.write(0, 9, 'Match 2 Start', bold)
    sheet.write(0, 10, 'Match 2 Stop', bold)
    sheet.write(0, 11, 'In Control', bold)
    sheet.write(0, 12, 'Organism', bold)
    sheet.write(0, 13, 'Direction', bold)
    sheet.write(0, 14, 'ScanID', bold)
    sheet.write(0, 15, 'ALC Score', bold)


def get_formats(book):
    """
    Sets and returns formatting for Excel file.
    """
    # Excel formatting
    bold = book.add_format({'bold': True, 'align': 'right', 'font_name': 'Times New Roman'})
    italic = book.add_format({'italic': True, 'align': 'right', 'font_name': 'Times New Roman'})
    blue = book.add_format({'font_color': 'blue', 'align': 'right', 'font_name': 'Times New Roman'})
    red = book.add_format({'font_color': 'red', 'align': 'right', 'font_name': 'Times New Roman'})
    green = book.add_format({'font_color': 'green', 'align': 'right', 'font_name': 'Times New Roman'})
    regular = book.add_format({'align': 'right', 'font_name': 'Times New Roman'})
    
    lbold = book.add_format({'bold': True, 'align': 'left', 'font_name': 'Times New Roman'})
    lblue = book.add_format({'font_color': 'blue', 'align': 'left', 'font_name': 'Times New Roman'})
    lred = book.add_format({'font_color': 'red', 'align': 'left', 'font_name': 'Times New Roman'})
    lgreen = book.add_format({'font_color': 'green', 'align': 'left', 'font_name': 'Times New Roman'})
    return bold, italic, blue, red, green, regular, lbold, lblue, lred, lgreen


def combine_output(strain, directory, peptides_data, control_peptides):
    """
    Writes out an Excel file.
    """
    print
    try:
        directory_out = os.getcwd() + os.sep + 'analysed' + os.sep
        os.makedirs(directory_out)
    except OSError:
        pass
    
    # determine experiments to write out and group them
    peptide_files_by_experiment = {}
    for f in sorted(os.listdir(directory)):
        experiment = f.split('_')[0]
        peptide_files_by_experiment[experiment] = peptide_files_by_experiment.get(experiment, []) + [f]
    
    print time_stamp(), 'Writing out to file:', color.BOLD, '%s_full_proteome_scan.xls' % strain, color.END
    book = xlsxwriter.Workbook('%s%s_full_proteome_scan.xls' % (directory_out, strain))
    bold, italic, blue, red, green, regular, lbold, lblue, lred, lgreen = get_formats(book)
    sheet_n = 1
    sheet = book.add_worksheet('Hits')
    _write_header(sheet, lbold)
    row = 1
    
    unique_peptides_per_experiment = {}
    
    for experiment, files in peptide_files_by_experiment.items():
        unique_peptides_per_experiment[experiment] = unique_peptides_per_experiment.get(experiment, set())
        
        for f_name in files:
            for line in open(directory + f_name):
                raw_data = line.strip('\n').split(',')
                
                # convert fields to integers where possible
                data = []
                for val in raw_data:
                   try:
                       d = int(val)
                   except ValueError:
                       d = val
                   data.append(d)
                
                # skip outputting full length match entries
                if 'N/A' == data[6]:
                    continue
                
                # prevent Excel from complaining
                if row > 60000:
                    sheet_n += 1
                    sheet = book.add_worksheet('Hits' + '_' + str(sheet_n))
                    _write_header(sheet, bold)
                    row = 1
                
                sheet.write(row, 0, experiment, regular)  #Experiment
                sheet.write(row, 1, data[0], lbold)       #Peptide
                sheet.write(row, 2, data[1], regular)     #Peptide Length
                sheet.write(row, 3, data[2], regular)     #Match 1
                sheet.write(row, 4, data[3], lred)        #Fragment 1
                sheet.write(row, 5, data[4], red)         #Match 1 Start
                sheet.write(row, 6, data[5], red)         #Match 1 Stop
                sheet.write(row, 7, data[6], italic)      #Cis Distance 1
                sheet.write(row, 8, data[7], lgreen)      #Fragment 2
                sheet.write(row, 9, data[8], green)       #Match 2 Start
                sheet.write(row, 10, data[9], green)      #Match 2 Stop
                
                if not control_peptides.has_key(experiment):
                    sheet.write(row, 11, 'N/A', regular)  #Present in control
                elif data[0] in control_peptides[experiment]:
                    sheet.write(row, 11, 'Yes', regular)  #Present in control
                else:
                    sheet.write(row, 11, 'No', bold)      #Present in control
                
                if 'NL4.3' == data[2][:5] or 'IIIB' == data[2][:4]:
                    sheet.write(row, 12, 'HIV', bold)  #Present in HIV
                else:
                    sheet.write(row, 12, 'Host/Other', regular)  #Present in host
                
                if data[4] < data[8]:
                    sheet.write(row, 13, 'Normal', bold)
                else:
                    sheet.write(row, 13, 'Reverse', regular)
                
                sheet.write(row, 14, data[10], regular)
                sheet.write(row, 15, data[11], regular)
                
                row += 1
                unique_peptides_per_experiment[experiment].add(data[0])
    book.close()
    
    
    
    print time_stamp(), 'Cleaning up temporary files directory ...'
    shutil.rmtree(directory)
    print time_stamp(), 'Completed!\n'
    
    print '%sData points explained:%s' % (color.BOLD, color.END)
    for experiment, values in unique_peptides_per_experiment.items():
        print '%s: %4d / %4d explained peptides' % (experiment, len(values), len(peptides_data[experiment]))


if __name__ == '__main__':
    peptides_file = os.path.abspath('data/HLA-Alleles-DeNovoSearch.csv')
    control_peptides_file = os.path.abspath('data/Control-HLA-Alleles-DeNovoSearch.csv')
    
    # parameters for processing
    min_peptide_length = 7
    max_peptide_length = 25
    max_fusion_search_distance = 40
    
    # multithreading parameters
    max_peptides_per_group = 1
    max_CPUs = 48
    
    
    ### Analyse data ###
    strain = 'CONTROL'     # must set strain parameter manually
    reference_file = os.path.abspath('data/SPHu.fasta')
    
    # get data
    control_peptides = read_in_peptides(control_peptides_file, min_peptide_length, max_peptide_length, strain)
    peptides_data = read_in_peptides(peptides_file, min_peptide_length, max_peptide_length, strain)
    references_data = read_in_references(reference_file)
    
    # get matches and establish which peptides can be fused
    multiprocess_peptide_search(strain, peptides_data, max_CPUs, references_data, max_fusion_search_distance)
    
    # combine results
    combine_output(strain, os.getcwd() + os.sep + 'multithreaded_temp_files' + os.sep, peptides_data, control_peptides)
