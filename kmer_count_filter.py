# @file kmer_count_filter.py
# @author Esko Kautto (esko.kautto@osumc.edu)
# @updated 2018-05-23

import os
import argparse
import math
import numpy
from datetime import datetime, timedelta
from helpers import iteritems, timestamp, tprint


class MANTIS_Filter(object):
    class Locus(object):
        def __init__(self):
            self.locus = ''
            self.k = set()
            self.n = {}
            self.t = {}

        def add(self, line):
            if type(line) is str:
                line = line.strip().split('\t')
            if self.locus == '':
                self.locus = line[0]
            elif line[0] != self.locus:
                tprint('MANTIS_Filter.Locus error: Fed data for ' + 
                    '{0} to locus {1}!'.format(line[0], self.locus))
                return False
            k = int(line[1])
            n = int(line[2])
            t = int(line[3])
            if k in self.k:
                tprint('MANTIS_Filter.Locus error: Duplicate entry for ' + 
                    'repeat count {0}!'.format(k))
                return False

            self.k.add(k)
            self.n[k] = n
            self.t[k] = t
            return True
            # end .add()

        """
        Generates output data so locus can be written into output file.
        """
        def generate_output(self):
            output = []
            for k in sorted(self.k):
                n = 0
                t = 0
                if k in self.n:
                    n = self.n[k]
                if k in self.t:
                    t = self.t[k]
                if (n + t) > 0:
                    line = [self.locus, str(k), str(n), str(t)]
                    output.append(line)
            
            if len(output) == 0:
                # Nothing to output.
                return False
            return output
            # end .generate_output()


        """
        Checks if both the normal and tumor data have 
        the minimum required coverage.
        """
        def has_min_coverage(self, min_coverage):
            if (sum(self.n.values()) < min_coverage) \
                or (sum(self.t.values()) < min_coverage):
                return False
            return True
            # end .has_min_coverage()


        """
        Generates a list of values from the repeat counts and
        support for each repeat count.
        """
        def generate_value_list(self,data):
            values = []

            # Iterate through the data and only record values
            # for k-repeats that have some support
            for k, count in iteritems(data):
                if count > 0:
                    values += [k] * count
            return sorted(values)
            # end .generate_value_list()

        """
        Filters the data for outliers. The normal and tumor
        data are treated as separate sets for the outlier 
        calculations, since each will likely have a different
        mean if the locus is unstable.
        """
        def outlier_filter(self, sds):
            # Process each subset of data separately.
            self.n = self.subset_outlier_filter(self.n, sds)
            self.t = self.subset_outlier_filter(self.t, sds)
            self.remove_unsupported_ks()
            self.regenerate_k_set()
            return True
            # end .outlier_filter()

        """
        Zeroes out any k-value repeat counts that don't meet the 
        minimum amount of supported reads; mainly gets rid of
        one-off values that are likely to only exist due to
        sequencing or alignment errors in the reads.
        """
        def support_filter(self, min_support):
            for k in self.k:
                n_support = False
                t_support = False
                if k in self.n and self.n[k] >= min_support:
                    n_support = True
                if k in self.t and self.t[k] >= min_support:
                    t_support = True

                if not n_support and not t_support:
                    # Neither subset had required support.
                    self.n[k] = 0                    
                    self.t[k] = 0
            self.remove_unsupported_ks()
            self.regenerate_k_set()
            # end .support_filter()

        """
        Regenerates the set of k-repeat values based on the
        existing normal and tumor data.
        """
        def regenerate_k_set(self):
            self.k = set()
            for k in self.n.keys():
                self.k.add(k)
            for k in self.t.keys():
                self.k.add(k)
            return True
            # end .regenerate_k_set()

        """
        Removes any k-repeat values that no longer have any support.
        """
        def remove_unsupported_ks(self):
            unsupported = set()
            for k in self.k:
                n = 0
                t = 0
                if k in self.n:
                    n = self.n[k]
                if k in self.t:
                    t = self.t[k]
                if (n + t) == 0:
                    # Has no support
                    unsupported.add(k)
            if len(unsupported):
                for k in unsupported:
                    self.k.discard(k)
            return True
            # end .remove_unsupported_ks()

        """
        Keeps only the values that fall within the acceptable
        window, removing outliers.
        """
        def subset_outlier_filter(self, data, sds):
            values = self.generate_value_list(data)
            output = {}
            if len(values):
                mean = numpy.mean(values)
                std = numpy.std(values)
                # Round the min and max for the window to allow for
                # some leniency in the filter.
                min_k = int(math.floor(mean - (sds * std)))
                max_k = int(math.ceil(mean + (sds * std)))

                for k, count in iteritems(data):
                    if min_k <= k <= max_k:
                        # Acceptable
                        output[k] = count
            return output
            # end .subset_outlier_filter()

        # end MANTIS_Filter.Locus subclass definition





    def __init__(self, config):
        self.min_locus_coverage = config['min_locus_coverage']
        self.min_repeat_reads = config['min_repeat_reads']
        self.standard_deviations = config['standard_deviations']
        # end .__init__()

    """
    Main filtering function that calls all other methods
    required by the filter.
    """
    def filter(self, input_filepath, output_filepath):
        self.load_loci(input_filepath)
        self.filtered_loci = {}
        for locus_name, locus in iteritems(self.loci):
            locus = self.filter_locus(locus)
            if locus is not False:
                self.filtered_loci[locus_name] = locus
        self.write_output(output_filepath)
        # end .filter()

    """
    Runs the different locus filters on the locus instance.
    """
    def filter_locus(self, locus):
        locus.outlier_filter(self.standard_deviations)
        locus.support_filter(self.min_repeat_reads)
        if not locus.has_min_coverage(self.min_locus_coverage):
            return False
        return locus
        # end .filter_locus()


    """
    Writes output data into specified file, iterating
    through each (filtered) locus, generating the output
    data for it and writing it out.
    """
    def write_output(self, output_filepath):
        fileout = open(output_filepath, 'w')
        header = ['Locus', 'Repeats', 'Normal', 'Tumor']
        fileout.write('\t'.join(header) + '\n')
        for l in self.ordered_loci:
            if l in self.filtered_loci:
                locus = self.filtered_loci[l]
                output = locus.generate_output()
                for line in output:
                    fileout.write('\t'.join(line) + '\n')
        fileout.close()
        return True
        # end .write_output()


    """
    Loads the data from the input file and stores the per-locus
    data for filtering.
    """
    def load_loci(self, input_filepath):
        self.loci = {}
        self.ordered_loci = []
        with open(input_filepath, 'Ur') as f:
            header = True
            for line in f:
                if header:
                    # Make sure this is a header line so we don't
                    # discard valid data.
                    if self.line_is_header(line):
                        # Line is header line, skip
                        continue
                line = line.strip().split('\t')
                if line[0] not in self.loci:
                    self.loci[line[0]] = MANTIS_Filter.Locus()
                    self.ordered_loci.append(line[0])
                self.loci[line[0]].add(line)
        # end .load_loci()

    """
    Tries to determine if line is a header line based on column values.
    """
    def line_is_header(self, line):
        line = line.upper().strip().split('\t')
        if len(line) != 4:
            # Malformed line, definitely not data.
            return True
        header = ['LOCUS', 'REPEATS', 'NORMAL', 'TUMOR']
        is_header = True
        for i in range(0, 4):
            if line[i] != header[i]:
                # Unexpected token.
                is_header = False
                break
        return is_header
        # end .line_is_header()


    # end MANTIS_Filter class definition

"""
Attempts to filter out undesirable noise in the filter, in
an attempt to compare results with more confidence. The filter
will try to get rid of reads that are obvious outliers, loci that
don't meet a minimum coverage depth requirement, and repeat counts
that don't have enough supporting reads.
"""
if __name__ == "__main__":
    prog_name = 'MSI Locus Kmer Counter Filter'
    tprint(prog_name)

    parser = argparse.ArgumentParser(description=prog_name)

    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
        help='Input file (.kmer_counts)')

    parser.add_argument('-o', '--output', dest='output', type=str, required=True,
        help='Output filename.')

    parser.add_argument('-mlc', '--min-locus-coverage', dest='mlc', type=int,
        default=20, help='Minimum coverage required for each of the normal ' +
        'and tumor results.')

    parser.add_argument('-mrr', '--min-repeat-reads', dest='mrr', type=int,
        default=5, help='Minimum reads supporting a specific repeat count.')

    parser.add_argument('-sd', '--standard-deviations', dest='sd', type=float,
        default=3.0, help='Standard deviations from mean before repeat count is ' +
        'considered an outlier')

    args = parser.parse_args()



        
    input_filepath = os.path.abspath(args.input)
    if not os.path.isfile(input_filepath):
        tprint('Error: Input file {0} does not exist!'.format(input_filepath))
        exit(1)

    output_filepath = os.path.abspath(args.output)

    if args.mlc < 0:
        tprint('Error: Minimum locus coverage cannot be below 0.')
        exit(1)

    if args.mrr < 0:
        tprint('Error: Minimum read count cannot be below 0.')
        exit(1)        

    if args.sd < 0.0:
        tprint('Error: Standard deviation count cannot be below 0.')
        exit(1)        


    config = {
        'min_locus_coverage': int(args.mlc),
        'min_repeat_reads': int(args.mrr),
        'standard_deviations': float(args.sd)
    }

    locus_filter = MANTIS_Filter(config)
    locus_filter.filter(input_filepath, output_filepath)
    # Done.
    exit(0)
