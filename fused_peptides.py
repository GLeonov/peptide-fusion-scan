"""
Classes and methods for peptide fusion scanning and output.
"""

import itertools
import re
import uuid


class Peptide:
    def __init__(self, name, scan_id, score):
        self.name = name
        self.scanID = scan_id
        self.score = score
        self.unique_identifier = str(uuid.uuid4())
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


# methods for peptide fusion analysis and running
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
    il_variants = _create_combinations_with_isoleucines(peptide.name)
    for fragment in il_variants:
        fragments.append(FragmentPair(['', fragment]))

    # generate spliced peptide fragments
    for IL_variant in il_variants:
        for position in range(2, len(IL_variant) - 1):
            fragments.append(FragmentPair([IL_variant[:position], IL_variant[position:]]))

    # sort fragments, longest fragments (full length) first,
    # necessary for excluding them from splice searches if they match database in full
    fragments.sort(key=lambda x: len(x.pair[-1]), reverse=True)

    return fragments


def search_for_fragment_in_protein_by_splitting_orfs(peptide, references_data, max_search_distance=40):
    """
    Creates the reference search space for the sequence. Updates the Peptide
    object's Fragment objects with any Match objects it finds.

    First fragment to be searched is full length peptide itself - which means if matched,
    the search will stop before next fragment set gets a chance to be matched.

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

                    match_large = Match(protein_name, match_data_large.start() + distance,
                                        match_data_large.end() - 1 + distance, number)

                    # determine the search space start for small fragment matching
                    starting_position = match_data_large.start() - (max_search_distance + len(small_fragment))
                    if starting_position < 0:
                        starting_position = 0

                    # full length matches
                    if len(small_fragment) == 0:
                        full_length_match_found = True
                        fragment.matchPairs.append(['N/A', match_large])
                        continue

                    for match_data_small in re.finditer(small_fragment, ORF[starting_position:match_data_large.end() + 1 + max_search_distance + len(small_fragment)]):

                        match_small = Match(protein_name, starting_position + match_data_small.start() + distance,
                                            starting_position + match_data_small.end() - 1 + distance, number)

                        # discard results if large fragment match overlaps with small fragment match
                        if match_small.stop >= match_large.start and match_small.start <= match_large.stop:
                            continue

                        if len(fragment.pair[1]) < len(fragment.pair[0]):
                            fragment.matchPairs.append([match_large, match_small])
                        else:
                            fragment.matchPairs.append([match_small, match_large])

                distance += len(ORF) + 1
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
            if match_1.stop + 1 == match_2.start:
                return 0  # skip output if fusion is not sensible
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


def write_temp_file(f, fused_peptide, peptide, unmatched=False):
    """
    Write a temporary text file for later.
    """
    if unmatched:
        f.write('%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%s,%s\n' %
                (fused_peptide,
                 len(fused_peptide),
                 'N/A',
                 'N/A',
                 'N/A',
                 'N/A',
                 'N/A',
                 'N/A',
                 'N/A',
                 'N/A',
                 peptide.scanID,
                 peptide.score,
                 'de novo non-matched',
                 peptide.unique_identifier))
        return

    if fused_peptide.cis_distance_12 != 'N/A':
        f.write('%s,%d,%s,%s,%d,%d,%s,%s,%d,%d,%d,%d,%s,%s\n' %
                (fused_peptide.name,
                 len(fused_peptide.name),
                 fused_peptide.match_1_name,
                 fused_peptide.fragment_1,
                 fused_peptide.match_1_start + 1,
                 fused_peptide.match_1_stop + 1,
                 str(fused_peptide.cis_distance_12),
                 fused_peptide.fragment_2,
                 fused_peptide.match_2_start + 1,
                 fused_peptide.match_2_stop + 1,
                 peptide.scanID,
                 peptide.score,
                 'de novo matched and fused',
                 peptide.unique_identifier))
    else:
        f.write('%s,%d,%s,%s,%d,%d,%s,%s,%s,%s,%d,%d,%s,%s\n' %
                (fused_peptide.name,
                 len(fused_peptide.name),
                 fused_peptide.match_1_name,
                 fused_peptide.fragment_1,
                 fused_peptide.match_1_start + 1,
                 fused_peptide.match_1_stop + 1,
                 'N/A',
                 'N/A',
                 'N/A',
                 'N/A',
                 peptide.scanID,
                 peptide.score,
                 'de novo DB-matched',
                 peptide.unique_identifier))


def find_peptide_fusions(experiment, peptide, n_peptide, directory, references_data, max_search_distance):
    """
    Finds peptides that could be fused and outputs their matches.
    """
    peptide.fragments = _fragment_peptide(peptide)

    # find fragments in experimental data
    peptide = search_for_fragment_in_protein_by_splitting_orfs(peptide, references_data, max_search_distance)

    f_name = '%s%s_%d.txt' % (directory, experiment, n_peptide)
    f = open(f_name, 'w')

    # outputting the fusions
    matched = False
    for fragment in peptide.fragments:
        for match_pair in fragment.matchPairs:
            matched = True
            create_fused_peptide(f, fragment.pair[0], fragment.pair[1], match_pair[0], match_pair[1], peptide)

    # store the peptides that were not matched as fused or to the reference proteome
    if not matched:
        write_temp_file(f, peptide.name, peptide, unmatched=True)
    f.close()


def get_direction(start, end, origin):
    if origin == 'de novo matched and fused':
        if start < end:
            return 'Forward'
        else:
            return 'Reverse'
    else:
        return 'N/A'


def get_fused_fragments(direction, fragment_1, fragment_2):
    if direction != 'N/A':
        return fragment_1 + fragment_2
    else:
        return 'N/A'
