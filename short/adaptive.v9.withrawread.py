import sys

import pysam
import pandas as pd
import argparse

def checkNoRepeats(read_name, target_chr, target_start, r1_to_check, suboptimalThreshold):
    for r_aln in r1_to_check[read_name]:
        if r_aln.reference_name == target_chr and r_aln.reference_start + 1 == int(target_start):
            if r_aln.has_tag("AS") and r_aln.has_tag("XS"):
                AS = r_aln.get_tag("AS")
                XS = r_aln.get_tag("XS")
                if AS >= XS:  # only if the best hit alignment score is better than subotimal one, go ahead, else this read is not valid
                    # find delta: as - (xs/as*100)
                    delta = 100.0 - (XS / float(AS) * 100)
                    # compare delta with threshold
                    if delta >= float(suboptimalThreshold):
                        return True
                return False
            else:
                return True
    print("WTF: ", read_name)
    return False
def getListFromCigar2(cigar):
    '''
    Generate a list of tuples with integer value and relative cigar operation character
    @param cigar: The CIGAR string
    @return: A list of tuples of the type (int, op)
    '''
    pos = 0
    cigar_listf = []
    for i in range(len(cigar)):
        if cigar[i].isalpha():
            cigar_listf.append((int(cigar[pos:i]), cigar[i]))
            pos = i + 1
    return cigar_listf


def getListFromCigar(cigar):
    '''
    Generate a list of tuples with integer value and relative cigar operation character
    @param cigar: The CIGAR string
    @return: A list of tuples of the type (int, op)
    '''
    pos = 0
    cigar_listf = []
    for i in range(len(cigar)):
        if cigar[i].isalpha():
            if cigar[i] != 'I' and cigar[i] != 'D':
                cigar_listf.append((int(cigar[pos:i]), cigar[i]))
            elif cigar[i] == 'D':
                cigar_listf.append((int(cigar[pos:i]), 'M'))
            pos = i + 1
    cleaned_cigar_list = []
    cleaned_cigar_list.append(cigar_listf[0])
    j = 0
    for i in range(1, len(cigar_listf)):
        if cigar_listf[i][1] == cigar_listf[i - 1][1]:
            matches = cleaned_cigar_list[j][0] + cigar_listf[i][0]
            cleaned_cigar_list[j] = (matches, cigar_listf[i][1])
        else:
            cleaned_cigar_list.append(cigar_listf[i])
            j += 1

    return cleaned_cigar_list

def getListFromCigarOnQuery(cigar):
    '''
    Generate a list of tuples with integer value and relative cigar operation character.
    The only difference between this and getListFromCigar is that here I becomes a M and not D.
    This because as explained in the SAM format I consumes query but not reference, for D is the opposite.
    This is useful for gap evaluation and alignment_on_query
    @param cigar: The CIGAR string
    @return: A list of tuples of the type (int, op)
    '''
    pos = 0
    cigar_listf = []
    for i in range(len(cigar)):
        if cigar[i].isalpha():
            if cigar[i] != 'I' and cigar[i] != 'D':
                cigar_listf.append((int(cigar[pos:i]), cigar[i]))
            elif cigar[i] == 'I':
                cigar_listf.append((int(cigar[pos:i]), 'M'))
            pos = i + 1
    cleaned_cigar_list = []
    cleaned_cigar_list.append(cigar_listf[0])
    j = 0
    for i in range(1, len(cigar_listf)):
        if cigar_listf[i][1] == cigar_listf[i - 1][1]:
            matches = cleaned_cigar_list[j][0] + cigar_listf[i][0]
            cleaned_cigar_list[j] = (matches, cigar_listf[i][1])
        else:
            cleaned_cigar_list.append(cigar_listf[i])
            j += 1

    return cleaned_cigar_list


def getStartingMismatch(cigar_listf):
    '''
    Find the number of the starting Mismatch (including insertion and deletion) and position of the
    first Match in a cigar string starting from a cigar list of the type [(int, op), ...]
    A tuple containing number of mismatch and position of the first match is returned.
    We can use the position of the first match to slice the cigar list e.g.:
    >>> cigar_list = getListFromCigar("32S21M20S")
    >>> print(cigar_list)
    [(32, 'S'), (21, 'M'), (20, 'S')]
    >>> no_match, pos = getStartingMismatch(cigar_list)
    >>> print("Cigar tuples before first match: " cigar_list[:pos])
    Cigar tuples before first match: [(32, 'S')]
    >>> print("First match: ", cigar_list[pos])
    First match: (21, 'M')

    @param cigar_listf: list of cigar tuples as returned from getListFromCigar
    @return: a tuple of number of mismatch and position of the first Match
    '''
    no_match = 0
    for i in range(len(cigar_listf)):
        bases, cigar_type = cigar_listf[i]
        if cigar_type != 'M':
            no_match += bases
        else:
            return (no_match, i)


def getClosingMismatch(cigar_listf):
    '''
    Find the number of Mismatch (including insertion and deletion) and the position of the first Mismatch
    after the last Match at the end of a cigar string starting from a cigar list of the type [(int, op), ...].
    This is simply obtained by passing the reversed list to getStartingMismatch.
    We can use the position of the first mismatch after the last match to slice the cigar list e.g.:
    >>> cigar_list = getListFromCigar("32S21M20S")
    >>> print(cigar_list)
    [(32, 'S'), (21, 'M'), (20, 'S')]
    >>> no_match, pos = getClosingMismatch(cigar_list)
    >>> print("Cigar tuples until last match: " cigar_list[:pos])
    Cigar tuples until last match: [(32, 'S'), (21, 'M')]
    >>> print("Cigar tuples after last match: ", cigar_list[pos:])
    Cigar tuples after last match: [(20, 'S')]
    >>> print("Last Match: ", cigar_list[pos-1])
    Last Match: (21, 'M')

    @param cigar_listf: list of cigar tuples as returned from getListFromCigar
    @return: a tuple of number of mismatch at the end of the cigar string
             and position of the first Mismatch after the last match
    '''
    no_match, pos = getStartingMismatch(cigar_listf[::-1])
    return no_match, len(cigar_listf) - pos


def checkChrPos(dict_list, chr_to_check, pos_to_check):
    for i_chr, i_pos in dict_list:
        if i_chr == chr_to_check and -3000 < pos_to_check - i_pos < 3000:
            return True
    return False


def getMatchesFromCigarString(cigar_str):
    cigar_list = getListFromCigar(cigar_str)
    matches = [int(x[0]) for x in cigar_list if x[1] == 'M']
    return sum(matches)

def getMatchesFromCigarStringOnQuery(cigar_str):
    cigar_list = getListFromCigarOnQuery(cigar_str)
    matches = [int(x[0]) for x in cigar_list if x[1] == 'M']
    return sum(matches)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input bam file')
parser.add_argument('-o', '--output', help='path/basename for outputs')
parser.add_argument('-F', '--input_F4', help='input F4 bam file')
parser.add_argument('-l', '--lower_gap', help='maximum negative gap accepted')
parser.add_argument('-g', '--higher_gap', help='maximum positive gap accepted')
# parser.add_arxgument('-R', '--input_F320', help='input F320 bam file')
args = parser.parse_args()
file = args.input
file_F4 = args.input_F4
# file_F320 = args.input_F320
output = args.output
max_gap = int(args.higher_gap)
min_gap = -int(args.lower_gap)

"""
file = "/Users/cipriani.carlo/Desktop/old_AAV-Short/GROMPE/v5/analysis/LTR51.LC30.sorted.md.rel.pg.iss.bam"
file_F4 = "/Users/cipriani.carlo/Desktop/old_AAV-Short/GROMPE/v5/analysis/LTR51.LC30.F4.sorted.md.bam"
output = "/Users/cipriani.carlo/Desktop/old_AAV-Short/GROMPE/v5/analysis/prove.newfilter."
max_gap = int(50)
min_gap = -int(50)
"""
aav_gap_threshold = 50
suboptimalThreshold = 40
# file = "/home/carlo/Scrivania/Master_Degree/Tesi/adaptative/bams/cleaned.26.AAV53.LC50.sorted.md.rel.pg.iss.bam"
dict_r2 = {}
r1_step3 = {}
r2_step3 = {}
r2_samfile = pysam.AlignmentFile(file_F4, "rb")
for read in r2_samfile.fetch():
    if read.is_read2:
        if read.query_name in dict_r2:
            dict_r2[read.query_name].append([read.reference_name, read.reference_start + 1])
            r2_step3[read.query_name].append(read)
        else:
            dict_r2[read.query_name] = [[read.reference_name, read.reference_start + 1]]
            r2_step3[read.query_name] = [read]
    else:
        if read.query_name in r1_step3:
            r1_step3[read.query_name].append(read)
        else:
            r1_step3[read.query_name] = [read]

samfile = pysam.AlignmentFile(file, "rb")
list_results = []
list_only_rear = []
list_no_integration = [[], []]
list_results_no_integration = []
list_repeats = []
r1_counter = 0
r1_chimeric = 0
r1_noresult = 0  # reads which map not on chrV and are not chimeric
r1_no_target_counter = 0
r1_onlyaav_counter = 0
r1_except_counter = 0
# pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)
for read in samfile.fetch():
    if read.is_read1:
        r1_counter += 1
        try:
            if read.has_tag("SA"):
                    """
                    If read has secondary alignment we are interested in it.
                    We create a list of all the alignments: SA + primary.
                    Each alignment is represented by a list:
                    [chr, pos, strand, cigarstring, map_quality, NM (not used)]
                    So we will have a list of lists.
                    We then split this list in chrV alignments and others (every chr except chrV)
                    """
                    sa = read.get_tag("SA")
                    all_alignments = [nl.split(",") for nl in sa.split(";")[:-1]]
                    strand = "-" if read.is_reverse else "+"
                    all_alignments.append(
                        [read.reference_name, read.reference_start + 1, strand, read.cigarstring, read.mapq, 0])
                    chrV_alignments = []
                    others_alignments = []
                    for i in range(len(all_alignments)):
                        if all_alignments[i][0] == "chrV":
                            chrV_alignments.append(all_alignments[i])
                        else:
                            others_alignments.append(all_alignments[i])
                    results = {}
                    for x_oth in others_alignments:
                        x_oth.append(getMatchesFromCigarString(x_oth[3]))
                    # from now everything becomes unreadable, must add doc and maybe change list to dict
                    """
                    For each alignment on the vector genome chrV we look where on the read it begins so
                    we can sort the arrangements if any.
                    """
                    # starts_chrV is a list of tuples of the type:
                    # [((starting_mismatches, index_of the_first_match), index in the chrV alignments), ...]
                    # it is used to sort vector alignments to find the right junction point and the number of rear.

                    starts_chrV = []
                    chrV_cigars_list = []
                    for i in range(len(chrV_alignments)):
                        cigar_list = getListFromCigar(chrV_alignments[i][3])
                        chrV_cigars_list.append(cigar_list)
                        # We are interested in the position of the FIRST MATCH.
                        # If the read is reverse strand we use the getClosingMismatch.
                        # It returns the first mismatch after the last match.
                        # Using it minus one we obtain the first match
                        # e.g with a cigar string = 32S21M and the relative cigar list [(32, 'S'), (21, 'M')]
                        # it will return (0,2) this means 0 mismatch at the end and the first mismatch after
                        # the last match is in position 2 -> since the list has only two elements this means
                        # that the string ends with a match. To obtain the index position of the last match
                        # we can simply do pos - 1. For the negative strand the closing mismatches are equal
                        # to the starting mismatches.
                        c = getStartingMismatch(cigar_list) if chrV_alignments[i][2] == '+' else getClosingMismatch(
                            cigar_list)
                        if chrV_alignments[i][2] == '-':
                            c = (c[0], c[1] - 1)  # adjusting the pos if reverse strand as explained above
                        starts_chrV.append((c, i))

                    starts_chrV.sort(key=lambda tup: tup[0][0])
                    # each element of the list is a tuple like this: ((start, index), elem) as explained above
                    # start is the start of the alignment in the read.
                    # Index is the index of the first MATCH in the cigar list
                    # elem is the position of this chrV alignment in the list of the vector alignments
                    # stressing out on this since is not very intuitive at the beginning, maybe another structure
                    # would be better.
                    # Example of a read with two chrV alignments e.g. (reporting only useful data)
                    # chrV_alignments = [[chrV,3209,"42S90M","+", ...],["chrV","2018","90S42M","-",...]
                    # we will have two cigar_lists : [(42,S),(90,M)] and [(90,S),(42,M)]
                    # for the first one since is on the positive strand we will have:
                    # start = 42  (since we have 42 softclipped at the beginning of the read)
                    # index = 1   (since we have a softclip before the match)
                    # or the second we will have (it is on the negative strand):
                    # start = 0  (since we have a match in the last position)
                    # index = 1   (since the last match is in the last position)
                    # after processing and sorting
                    # starts_chrV = [((0,1),1),((42,1),0)]
                    start, index = starts_chrV[0][0]
                    elem = starts_chrV[0][1]
                    end = start + getMatchesFromCigarStringOnQuery(chrV_alignments[elem][3])
                    other_chrV_aln = len(chrV_alignments) - 1
                    results['name'] = read.query_name
                    results['n_aav_aln'] = 1
                    results['start_aav'] = start
                    results['aav_matches'] = chrV_cigars_list[elem][index][0]
                    aav_string = ",".join([str(x) for x in chrV_alignments[elem]])
                    #aav_string += ";"
                    last_chrV_alignment = chrV_alignments[elem]
                    last_matches = chrV_cigars_list[elem][index][0]
                    new_aav_string = ",".join(str(x) for x in ["chrV", chrV_alignments[elem][1],
                                                               int(chrV_alignments[elem][1]) + int(last_matches) - 1,
                                                               chrV_alignments[elem][2]])
                    alignment_on_query = ",".join(str(x) for x in ["chrV",start,end,chrV_alignments[elem][2]])

                    for i in range(1, other_chrV_aln + 1):
                        start, index = starts_chrV[i][0]
                        elem = starts_chrV[i][1]
                        aav_gap = start - end
                        if abs(aav_gap) <= aav_gap_threshold:
                            end = start + getMatchesFromCigarStringOnQuery(chrV_alignments[elem][3])
                            results['n_aav_aln'] += 1
                            aav_string += ";" + ",".join([str(x) for x in chrV_alignments[elem]])

                            last_chrV_alignment = chrV_alignments[elem]
                            last_matches = chrV_cigars_list[elem][index][0]
                            new_aav_string += ";" + ",".join(str(x) for x in ["chrV", chrV_alignments[elem][1],
                                                                              int(chrV_alignments[elem][1]) + int(
                                                                                  last_matches) - 1,
                                                                              chrV_alignments[elem][2]])

                            results['aav_matches'] += chrV_cigars_list[elem][index][0]
                            alignment_on_query += ";"+ ",".join(str(x) for x in ["chrV", start, end,chrV_alignments[elem][2]])
                            if aav_gap < 0:
                                results['aav_matches'] += aav_gap

                    results['input_file'] = file_F4
                    results['aav_alignments'] = aav_string
                    results['aav_alignments_start_end'] = new_aav_string
                    results['end_aav'] = end
                    results['aav_last_start'] = int(last_chrV_alignment[1])
                    results['aav_last_end'] = int(last_chrV_alignment[1]) + int(last_matches) - 1
                    results['aav_last_strand'] = last_chrV_alignment[2]
                    if last_chrV_alignment[2] == '+':
                        results['junction_locus'] = results['aav_last_end']
                    else:
                        results['junction_locus'] = results['aav_last_start']
                    starts_others = []
                    others_cigars_list = []
                    results['gap'] = None
                    results['other'] = None
                    others_alignments.sort(key=lambda x: int(x[-1]), reverse=True)
                    for i in range(len(others_alignments)):
                        cigar_list = getListFromCigar(others_alignments[i][3])
                        others_cigars_list.append(cigar_list)
                        start, index = getStartingMismatch(cigar_list) if others_alignments[i][
                                                                              2] == '+' else getClosingMismatch(cigar_list)
                        others_alignments[i].append(start - end)

                    # others_alignments.sort(key=lambda x: abs(int(x[-1])), reverse=False)

                    dict_r2_chr_list = [c[0] for c in dict_r2[read.query_name]]
                    for i in range(len(others_alignments)):
                        # cigar_list = getListFromCigar(others_alignments[i][3])
                        # others_cigars_list.append(cigar_list)
                        # start, index = getStartingMismatch(cigar_list) if others_alignments[i][2] == '+' else getClosingMismatch(cigar_list)
                        gap = others_alignments[i][-1]
                        if max_gap >= gap >= min_gap and others_alignments[i][0] in dict_r2_chr_list and checkChrPos(
                                dict_r2[read.query_name], others_alignments[i][0], int(others_alignments[i][1])):
                            if checkNoRepeats(read.query_name, others_alignments[i][0], others_alignments[i][1], r1_step3, suboptimalThreshold):
                                results['gap'] = gap
                                results['other'] = others_alignments[i]
                                results['target_chr'] = others_alignments[i][0]
                                results['target_start'] = others_alignments[i][1]
                                cigar_list = getListFromCigar(others_alignments[i][3])
                                start, index = getStartingMismatch(cigar_list) if others_alignments[i][
                                                                                      2] == '+' else getClosingMismatch(
                                    cigar_list)
                                if others_alignments[i][2] == '-':
                                    index -= 1
                                results['target_end'] = int(others_alignments[i][1]) + int(cigar_list[index][0]) - 1
                                results['target_strand'] = others_alignments[i][2]
                                results['target_cigar'] = others_alignments[i][3]
                                results['target_matches'] = int(cigar_list[index][0])
                                results['total_matches'] = results['aav_matches'] + results['target_matches']
                                results['from_junc_to_plusN'] = None
                                results['from_junc_to_plus20'] = None
                                results['from_minusN_to_IS'] = None
                                original_sequence = read.get_forward_sequence()
                                if others_alignments[i][2] == '+':
                                    results['integration_locus'] = results['target_start']
                                else:
                                    results['integration_locus'] = results['target_end']
                                if results['gap'] >= 0:
                                    results['seq_gap'] = original_sequence[
                                                         results['end_aav']:results['end_aav'] + results['gap']]
                                    results['from_junc_to_plusN'] = original_sequence[
                                                                    results['end_aav']:results['end_aav'] + 12]
                                    results['from_junc_to_plus20'] = original_sequence[
                                                                     results['end_aav']:results['end_aav'] + 20]
                                    loc_x = results['end_aav'] + results['gap']
                                    results['from_minusN_to_IS'] = original_sequence[loc_x - 12:loc_x]
                                else:
                                    results['total_matches'] += results['gap']
                                    results['seq_gap'] = original_sequence[
                                                         results['end_aav'] + results['gap']:results['end_aav']]
                                    loc_x = results['end_aav'] + results['gap']
                                    results['from_junc_to_plusN'] = original_sequence[loc_x:loc_x + 12]
                                    loc_x = results['end_aav'] + results['gap']
                                    results['from_minusN_to_IS'] = original_sequence[results['end_aav'] - 12:results['end_aav']]
                                # results['pool_tag'] = read.get_tag('RG')
                                results['integration'] = ":".join(str(x) for x in others_alignments[i][0:3])
                                alignment_on_query += ";" + ",".join(str(x) for x in [others_alignments[i][0], start, start+getMatchesFromCigarStringOnQuery(others_alignments[i][3]), others_alignments[i][2] ])
                                break
                            else:
                                list_repeats.append(read.query_name)
                    results['alignment_on_query'] = alignment_on_query
                    results['read_length'] = read.infer_read_length()
                    results['raw_read'] = read.get_forward_sequence()
                    if results['gap'] is not None:
                        r1_v2 = results['name']
                        list_reads_v2 = r2_step3[r1_v2]
                        r1_integration_v2 = results['integration'].split(":")
                        r1_integration_chr_v2 = r1_integration_v2[0]
                        r1_integration_pos_v2 = int(r1_integration_v2[1])
                        r1_integration_strand_v2 = r1_integration_v2[2]

                        for read_v2 in list_reads_v2:
                            strand_v2 = "-" if read_v2.is_reverse else "+"
                            if read_v2.reference_name == r1_integration_chr_v2 and (
                                    -3000 < read_v2.reference_start + 1 - r1_integration_pos_v2 < 3000):
                                target_alignment = [read_v2.reference_name, read_v2.pos, strand_v2, read_v2.cigarstring,
                                                    read_v2.mapq, 0]
                                results['cigar_r2'] = read_v2.cigarstring
                                results['pos_r2'] = read_v2.reference_start + 1
                                results['strand_r2'] = strand_v2
                                results['target_matches_r2'] = getMatchesFromCigarString(read_v2.cigarstring)
                                break

                        list_results.append(results)
                        r1_chimeric += 1
                    elif results['n_aav_aln'] > 1:
                        list_only_rear.append(results)
                        r1_onlyaav_counter += 1
                        if len(others_alignments) == 0:
                            r1_no_target_counter += 1
                    else:

                        list_no_integration[0].append(read.query_name)
                        list_no_integration[1].append(sa)
                        if read.reference_name == "chrV":
                            r1_onlyaav_counter += 1
                            if len(others_alignments) == 0:
                                r1_no_target_counter += 1

                        results = {}
                        strand = "-" if read.is_reverse else "+"
                        results['name'] = read.query_name
                        results['chr'] = read.reference_name
                        # results['cigar'] = read.cigarstring
                        results['aav_last_start'] = read.reference_start + 1
                        results['aav_last_end'] = read.reference_end
                        results['aav_last_strand'] = "-" if read.is_reverse else "+"
                        results['aav_alignments'] = ",".join([str(i) for i in
                                                              [read.reference_name, read.reference_start + 1, strand,
                                                               read.cigarstring, read.mapq, 0]])
                        results['aav_alignments_start_end'] = ",".join([str(results[res_var]) for res_var in
                                                                        ['chr', 'aav_last_start', 'aav_last_end',
                                                                         'aav_last_strand']])
                        results['input_file'] = file_F4
                        results['alignment_on_query'] = alignment_on_query
                        results['read_length'] = read.infer_read_length()
                        list_no_integration.append(read.query_name)
                        if read.reference_name == 'chrV':
                            results['start_aav'] = getStartingMismatch(getListFromCigar(read.cigarstring))[
                                0] if strand == "+" else getClosingMismatch(getListFromCigar(read.cigarstring))[0]
                            list_only_rear.append(results)
                        else:
                            results['start_aav'] = None
                            list_results_no_integration.append(results)
                            r1_noresult += 1

            else:

                    if read.reference_name == "chrV":
                        r1_onlyaav_counter += 1
                        r1_no_target_counter += 1

                    strand = "-" if read.is_reverse else "+"
                    results = {}
                    results['read_length'] = read.infer_read_length()
                    results['name'] = read.query_name
                    results['chr'] = read.reference_name
                    # results['cigar'] = read.cigarstring
                    results['aav_last_start'] = read.reference_start + 1
                    results['aav_last_end'] = read.reference_end
                    results['aav_last_strand'] = "-" if read.is_reverse else "+"
                    results['aav_alignments'] = ";".join([str(i) for i in
                                                          [read.reference_name, read.reference_start + 1, strand,
                                                           read.cigarstring, read.mapq, 0]])
                    matches =  results['aav_last_end'] - results['aav_last_start']
                    results['aav_alignments_start_end'] = ",".join(
                        [str(results[res_var]) for res_var in ['chr', 'aav_last_start', 'aav_last_end', 'aav_last_strand']])
                    results['input_file'] = file_F4
                    if read.reference_name == 'chrV':
                        results['start_aav'] = getStartingMismatch(getListFromCigar(read.cigarstring))[
                            0] if strand == "+" else getClosingMismatch(getListFromCigar(read.cigarstring))[0]
                        alignment_on_query = ",".join(str(x) for x in ["chrV", results['start_aav'], results['start_aav']+matches, results['aav_last_strand']])
                        results['alignment_on_query'] = alignment_on_query
                        list_only_rear.append(results)

                    else:
                        results['start_aav'] = None
                        list_results_no_integration.append(results)
                        r1_noresult += 1
                    list_no_integration.append(read.query_name)
        except Exception as ex:
            r1_except_counter += 1
            print()
            print(output, "EXCEPTION: ",read.query_name, ex )
            print()
df_r1 = pd.DataFrame(list_results)
try:
    df_r1['GCperc'] = df_r1.seq_gap.apply(lambda x: (x.count('G') + x.count('C')) / len(x) if len(x) > 3 else None)
except:
    pass
if not df_r1.empty:
    df_r1.to_csv(output + ".R1.tsv", sep="\t", index=False)

df_only_aav = pd.DataFrame(list_only_rear)
if not df_only_aav.empty:
    df_only_aav[
        ['name', 'aav_last_start', 'aav_last_end', 'aav_last_strand', 'aav_alignments', 'aav_alignments_start_end',
         'input_file', 'start_aav','read_length','alignment_on_query']].to_csv(output + "_only_aav.tsv", sep="\t", index=False)
df_noresults = pd.DataFrame(list_results_no_integration)
if not df_noresults.empty:
    df_noresults.to_csv(output + "_noresults.tsv", sep="\t", index=False)
with open(output + "_stats.txt", "w") as stats_file:
    stats_file.write(str(r1_counter) + "\t" + str(r1_chimeric) + "\t" + str(r1_onlyaav_counter) + "\t" + str(
        r1_noresult) + "\t" + str(r1_except_counter) + "\n")
print(output, "R1 parsed")
# sys.exit()
##### Start R2 processing ######
# samfile = pysam.AlignmentFile(file_F320, "rb")

# df_r1 = pd.read_csv(df_file,sep="\t")
try:
    r1_reads = df_r1['name'].unique()
    r1_dict = {}
    for line in list_results:
        key = line.pop('name')
        r1_dict[key] = line
except:
    from sys import exit

    print(output, "NO R1 integration found. Exiting")
    exit()
list_results = []
list_only_rear = []
list_no_integration = [[], []]
for r1 in r1_reads:
    results = {}
    list_reads = r2_step3[r1]
    r1_integration = r1_dict[r1]['integration'].split(":")
    r1_integration_chr = r1_integration[0]
    r1_integration_pos = int(r1_integration[1])
    r1_integration_strand = r1_integration[2]
    for read in list_reads:
        strand = "-" if read.is_reverse else "+"
        if read.reference_name == r1_integration_chr and (-3000 < read.reference_start + 1 - r1_integration_pos < 3000):
            target_alignment = [read.reference_name, read.pos, strand, read.cigarstring, read.mapq, 0]

            results['n_aav_aln'] = 0
            results['name'] = read.query_name
            results['start_aav'] = 0
            results['aav_aln'] = ";"
            results['end_aav'] = 0
            results['junction'] = None
            results['gap'] = 0
            results['other'] = target_alignment
            results['seq_gap'] = ""
            results['integration'] = ":".join(str(x) for x in target_alignment[0:3])
            list_results.append(results)
            break
"""
for read in samfile.fetch():
    results ={}
    if read.query_name in r1_reads:
        r1_integration = df_r1[df_r1['name'] == read.query_name].reset_index()['integration'][0]
        r1_integration_chr = r1_integration.split(":")[0]
        r1_integration_strand = r1_integration.split(":")[2]
        strand = "-" if read.is_reverse else "+"

        if read.has_tag("SA"):
            sa = read.get_tag("SA")
            all_alignments = [nl.split(",") for nl in sa.split(";")[:-1]]
            strand = "-" if read.is_reverse else "+"
            primary_alignment = [read.reference_name, read.pos, strand, read.cigarstring, read.mapq, 0]
            all_alignments.append(primary_alignment)
            chrV_alignments = []
            others_alignments = []
            for i in range(len(all_alignments)):
                if all_alignments[i][0] == "chrV":
                    chrV_alignments.append(all_alignments[i])
                else:
                    others_alignments.append(all_alignments[i])
            found = False
            for aln in others_alignments:
                if aln[0]==r1_integration_chr:
                    target_alignment = aln

                    results['n_aav_aln'] = 0
                    results['name'] = read.query_name
                    results['start_aav'] = 0
                    results['aav_aln'] = ";"
                    results['end_aav'] = 0
                    results['junction'] = None
                    results['gap'] = 0
                    results['other'] = target_alignment
                    results['seq_gap'] = ""
                    results['integration'] = ":".join(str(x) for x in target_alignment[0:3])
                    list_results.append(results)
                    break

        elif read.reference_name!="chrV" and read.reference_name==r1_integration_chr:
            strand = "-" if read.is_reverse else "+"
            primary_alignment = [read.reference_name, read.pos, strand, read.cigarstring, read.mapq, 0]
            results['n_aav_aln'] = 0
            results['name'] = read.query_name
            results['start_aav'] = 0
            results['aav_aln'] = ";"
            results['end_aav'] = 0
            results['junction'] = None
            results['gap'] = 0
            results['other'] = primary_alignment
            results['seq_gap'] = ""
            results['integration'] = ":".join(str(x) for x in primary_alignment[0:3])

            list_results.append(results)
        else:
            list_no_integration.append(read.query_name)
"""
df_r2 = pd.DataFrame(list_results)
try:
    df_r2['GCperc'] = df_r2.seq_gap.apply(lambda x: (x.count('G') + x.count('C')) / len(x) if len(x) > 3 else None)
except:
    pass
df_r2.to_csv(output + ".R2.tsv", sep="\t", index=False)
r1_reads = set(df_r1['name'].unique())
try:
    r2_reads = set(df_r2['name'].unique())
except:
    print(output, "R2 integrations not found")
    from sys import exit

    exit()

r2_not_found = r1_reads - r2_reads
with open(output + "noR2.log", "w") as logs:
    for r in r2_not_found:
        logs.write(r + "\n")
print(output, "R2 parsed")
"""out_samfile = pysam.AlignmentFile(output + "results.bam", "wb", template=samfile)
sam_file = pysam.AlignmentFile(file_F4, "rb")
for index, row in df_r1.iterrows():
    integration = row['integration'].split(":")
    chr = integration[0]
    pos = int(integration[1])
    strand = integration[2]
    for read in r1_step3[row['name']]:
        if read.reference_start + 1 <= pos + 1 and read.reference_start + 1 >= pos - 1 and read.reference_name == chr:

            read_to_write = read
            # read_to_write.query_alignment_start = read.query_alignment_start
            # read_to_write.query_alignment_end = read.query_alignment_end            #read_to_write.query_qualities = row['other'][4]
            read_to_write.is_read1 = True
            read_to_write.is_proper_pair = True
            if strand == '+':
                read_to_write.is_reverse = False
                read_to_write.mate_is_reverse = True
                read_to_write.flag = 99
            else:
                read_to_write.is_reverse = True
                read_to_write.mate_is_reverse = False
                read_to_write.flag = 83
            out_samfile.write(read_to_write)
            break

for index, row in df_r2.iterrows():
    integration = row['integration'].split(":")
    chr = integration[0]
    pos = int(integration[1])
    strand = integration[2]
    for read in r2_step3[row['name']]:
        if read.reference_start + 1 <= pos + 1 and read.reference_start + 1 >= pos - 1 and read.reference_name == chr:
            read_to_write = read
            read_to_write.is_read2 = True
            read_to_write.is_proper_pair = True
            if strand == '+':
                read_to_write.is_reverse = False
                read_to_write.mate_is_reverse = True
                read_to_write.flag = 163
            else:
                read_to_write.is_reverse = True
                read_to_write.mate_is_reverse = False
                read_to_write.flag = 147
            out_samfile.write(read_to_write)
            break
out_samfile.close()"""
print(output, " FINISHED")





