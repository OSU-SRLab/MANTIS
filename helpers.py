# @file helpers.py
# @author Esko Kautto (esko.kautto@osumc.edu)
# @updated 2016-06-16

import time

"""
Since the functionality of .iteritems() has changed between 
Python 2 and 3, this function attempts to bridge the gap and
provide similar functionality regardless of Python version.
"""
def iteritems(d):
    if hasattr(d, 'iteritems'):
        return d.iteritems()
    else:
        return list(d.items())
    # end .iteritems()

"""
Returns the current time in a m/d/y H:M:S format.
"""
def timestamp():
    return time.strftime('%m/%d/%y %H:%M:%S')
    # end .timestamp()

"""
Prints the input string's line(s) with prepended timestamp(s).
"""
def tprint(string):
    output = string.split('\t')
    for line in output:
        print('[{0}] {1}'.format(timestamp(), line))
    # end .tprint()

"""
Checks to make sure the BED file has the expected 6-column
format. The fourth column is expected to be a formatted
kmer repeat unit and count string. The fifth and sixth
column are not currently utilized, but are required to
maintain the BED 6-col format.
e.g.
chr1  10357206    10357223    (T)17   0   +
"""
def check_bedfile_format(filepath):
    with open(filepath, 'Ur') as f:
        line_number = 0
        for line in f:
            line = line.strip()
            line_number += 1
            if len(line):
                line = line.split('\t')
                if len(line) != 6:
                    print('Error: MANTIS expects a 6-column BED file with' + 
                        ' the 4th column containing the kmer repeat' + 
                        ' sequence and count (e.g. (T)15 or (CAC)5 ).')
                    print('\nOffending line (line {0}) has:'.format(line_number))
                    for n, value in enumerate(line):
                        print('[{0}]\t{1}'.format(n, value))

                    return False
                else:
                    # Check the k-mer column
                    kmer_format_ok = True
                    kmer = line[3]
                    if '(' in kmer and ')' in kmer:
                        kmer = kmer.split(')')
                        if '(' in kmer[0]:
                            kmer_unit = kmer[0].split('(')[1]
                            kmer_count = kmer[1]
                            if not kmer_count.isdigit():
                                kmer_format_ok = False
                            # Make sure only valid ATCG characters were used.
                            for c in list('ATCG'):
                                kmer_unit = kmer_unit.replace(c, '')
                            if len(kmer_unit) > 0:
                                kmer_format_ok = False
                        else:
                            kmer_format_ok = False
                    else:
                        kmer_format_ok = False


                    if not kmer_format_ok:
                        print('Error: MANTIS expects the kmer column in the' + 
                            ' BED file to follow a (XX)NN type format, with' +
                            ' the repeat unit wrapped in parentheses,' +
                            ' followed by the expected number of repeats.')
                        return False
        return True
    # end .check_bedfile_format()


"""
Determines how many times the k-mer repeats in a row
within the given sequence.
"""
def kmer_repeat_count(kmer, sequence, offset = 0):
    repeats = 0
    klen = len(kmer)
    while True:
        if sequence[offset:offset+klen] != kmer:
            # End of repeat sequence
            break
        repeats += 1
        offset += klen
        if offset > len(sequence):
            # End of sequence
            break
    return repeats
# end kmer_repeat_count()

# Checks to make sure required modules are present in environment
def required_modules_present(modules):
    missing = []
    for module in modules:
        try:
            __import__('imp').find_module(module.lower())
            # Everything is fine; module is available
        except ImportError:
            # Module not found
            missing.append(module)

    if len(missing):
        for module in missing:            
            print('Error: You must have {0} available in your environment!'.format(module))
        print('Please check your $PYTHONPATH to make sure you have properly ' +
            'included required moudles/libraries in it.')
        exit(1)
    return True
    # end required_modules_present()
