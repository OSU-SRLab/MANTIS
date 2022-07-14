#!/usr/bin/env python

# @file kmer_repeat_counter.py
# @author Esko Kautto (esko.kautto@osumc.edu)
# @updated 2017-12-18

import os
import argparse
import functools
import subprocess
import time
from datetime import datetime, timedelta
from multiprocessing import Process, Queue, BoundedSemaphore, Semaphore
import pysam
import numpy
from structures import SAMRead, CIGAR, Locus, MSILocus
from offset_finder import OffsetFinder
from helpers import iteritems, tprint, timestamp
from defaults import Parameter, load_settings
import re
import sys


class MSILocusLoader:
    def __init__(self, config):
        self.genome_path = None

        if 'genome' in config:
            self.genome_path = config['genome']
        else:
            tprint('MSILocusLoader warning: You need to specify a path' +
                ' to the reference genome!')
        # end .__init__()


    """
    Loads BED file into a list of Locus objects.
    """
    def load_loci(self, bedfile):
        loci = []
        bedfile = os.path.abspath(bedfile)
        if not os.path.isfile(bedfile):
            tprint('MSILocusLoader error: File {0}'.format(bedfile) + 
                ' does not exist!')
        else:
            with open(bedfile, 'r') as filein:
                for line in filein.readlines():
                    if line[0] != '@':
                        locus = MSILocus(line)

                        # Correct any off-by-one errors that may occur because of 
                        # unstandardized open- and closed-endedness of bed file coordinates.
                        self.correct_off_by_one_errors(locus)                    
                        loci.append(locus)
            filein.close()
        return loci
        # end .load_loci()

    strip_coord_re = re.compile(r'^[\d\-]+')
    """
    Fixes any off-by-one errors for the locus that could result from
    different interpretations of the open and closed ends of the BED
    format regions.
    """
    def correct_off_by_one_errors(self, locus):
        # Generate the locus position by adding a 1 basepair padding around the locus
        # coordinates, which accounts for the off-by-one errors.
        position = '{0}:{1}-{2}'.format(locus.chromosome, 
            locus.start - 1, 
            locus.end + 1)
        
        raw_sequence = self.get_sequence(position)
        #accommodate newer versions of PySam, which causes this to return chromosome and position
        raw_sequence = raw_sequence.split(":")[-1]
        sequence = MSILocusLoader.strip_coord_re.sub("", raw_sequence)
        if sequence[1:1+locus.kmer_length] != locus.kmer:
            # Sequence doesn't start where expected; shift accordingly.
            if sequence[0:locus.kmer_length] == locus.kmer:
                # Shift back by one
                locus.start -= 1
                locus.end -= 1
            elif sequence[2:2+locus.kmer_length] == locus.kmer:
                # Shift forward by one
                locus.start += 1
                locus.end += 1
            else:
                tprint('Error: Specified locus does not appear '
                    + ' to be the starting point for kmer {kmer}.'.format(
                        kmer=locus.kmer))
        # end .correct_off_by_one_errors()

    """
    Helper method that uses SAM tools to fetch the sequence of the 
    locus from the reference genome.
    """
    def get_sequence(self, locus):
        if self.genome_path is None:
            tprint('MSILocusLoader error: Can not use .get_sequence() without' 
                + ' specifying the path to the reference genome!')
            exit(1)
    
        sequence = []
        for subseq in pysam.faidx(self.genome_path, locus)[1:]:
            sequence.append(subseq.strip())
        return ''.join(sequence).upper()
        # end .get_sequence()    

    # end MSILocusLoader class definition.







class KmerRepeatCounter:
    def __init__(self, config):
        self.min_read_quality = config['min_read_quality']
        self.min_read_length = config['min_read_length']
        self.min_locus_quality = config['min_locus_quality']
        self.debug_output = config['debug_output']
        self.__reset()
        # end .__init__()

    def __reset(self):
        self.__consumers = list()
        self.__producer = None
        # end .__reset()


    """
    Estimates if the read is likely to be truncated within the target
    locus, which would indicate the repeat count will be too low.
    """
    def read_has_truncation(self, read, locus, offset, repeat_count):
        truncated = False
        if read.length - (offset + (repeat_count * locus.kmer_length)) < locus.kmer_length:
            # Possible 3' truncation
            diff = read.length - (offset + (repeat_count * locus.kmer_length))
            if read.seq[-1 * diff:] == locus.kmer[0:diff]:
                # Partial copy of the kmer; almost certainly truncated
                truncated = True

        elif offset < locus.kmer_length:
            # Possible 5' truncation
            if (offset is 0) or read.seq[0:offset] == locus.kmer[offset:]:
                # No offset OR offset matches fragment of kmer; likely truncation
                truncated = True
        return truncated
        # end .read_has_truncation()


    """
    Checks if the read passes the minimum QC requirements; currently
    consists of a minimum average read quality and read length.
    """
    def passes_qc_filter(self, read):
        passes = True
        if read.quality < self.min_read_quality:
            # Read quality too low.
            passes = False
        elif read.length < self.min_read_length:
            # Read sequence too short.
            passes = False
        return passes
        # end .passes_qc_filter()

    """
    Calculates the mean quality for the subset of the read.
    """
    @staticmethod
    def subset_mean_quality(read, start, end):
        scores = []
        for i in range(start, end):
            score = KmerRepeatCounter.symbol_to_quality_score(read.qual[i])
            scores.append(score)
        if len(scores):
            return numpy.mean(scores)
        return 0.0
        # end .subset_mean_quality()

    """
    Converts the symbolic quality score to a numeric (int) value.
    """
    @staticmethod
    def symbol_to_quality_score(symbol):
        symbols = '!"#$%&' + "'" + '()*+,-./0123456789:;<=>?@ABCDEFGHIJK'
        score = 0
        if symbol in symbols:
            score = symbols.find(symbol)
        return score
        # End .symbol_to_quality_score()


    """
    Checks if the repeating locus (as opposed to the entire read)
    has sufficient quality to pass the QC filter.
    """
    def passes_locus_qc_filter(self, read, offset_start, offset_end):
        locus_mean_quality = KmerRepeatCounter.subset_mean_quality(read, 
            offset_start, offset_end)
        if locus_mean_quality < self.min_locus_quality:
            # failed locus QC filter
            return False 
        return True
        # end .passes_locus_qc_filter()

    """
    Checks if any of the threads are still live.
    """
    def has_live_threads(self):
        live = False
        if self.has_live_producer():
            # The producer thread is still alive.
            live = True
        elif self.has_live_consumers():
            # Has at least one live consumer thread.
            live = True
        return live
        # end .has_live_threads()

    """
    Checks if the producer thread is still live.
    """
    def has_live_producer(self):
        return self.__producer.is_alive()
        # end .has_live_producer()

    """
    Checks if any of the consumer threads are still live.
    """
    def has_live_consumers(self):
        live = False
        for c in self.__consumers:
            if c.is_alive():
                live = True
                break
        return live
        # end .has_live_consumers()

    """
    Converts the numeric quality scores to the respective
    phred scale character values.
    """
    @staticmethod
    def quality_scores_to_symbols(scores):
        if type(scores) is int:
            return KmerRepeatCounter.score_to_symbol(scores)
        elif type(scores) is str:
            return KmerRepeatCounter.score_to_symbol(int(scores))
        else:
            symbols = []
            for score in scores:
                symbols.append(KmerRepeatCounter.score_to_symbol(int(score)))
            return ''.join(symbols)
        # end .quality_scores_to_symbols

    """
    Returns the character corresponding to the quality score.
    """
    @staticmethod
    def score_to_symbol(score):
        values = '!"#$%&' + "'" + '()*+,-./0123456789:;<=>?@ABCDEFGHIJK'
        if score >= len(values):
            symbol = values[-1]
        else:
            symbol = values[score]
        return symbol
        # end .score_to_symbol().

    """
    Determines how many times the k-mer repeats in a row
    within the given sequence.
    """
    @staticmethod
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
    # end .kmer_repeat_count()

    """
    Calculates the repeat count of the targeted k-mer (defined by the locus)
    within the read.
    """
    def locus_repeat_count(self, read, locus):
        repeat_count = -1
        offset = 0

        if read.start <= locus.start and read.end >= locus.end:
            offset = OffsetFinder.find_offset(read, locus)
            if offset < read.length:
                # Read covers locus; get repeat count
                repeat_count = KmerRepeatCounter.kmer_repeat_count(locus.kmer, read.seq, offset)

                # Check if the read might be truncated either on the 5' or 3'
                # end by comparing the offset and length of kmer to the length
                # of the read sequence.
                if self.read_has_truncation(read, locus, offset, repeat_count):
                    # Exclude this read
                    repeat_count = -1
        return (repeat_count, offset,)
        # end .locus_repeat_count()

    """
    Helper method utilizing PySam to determine if reads in the BAM file
    use a 'chr' prefix as part of the chromosomes or not.
    """
    @staticmethod
    def bam_uses_chr_prefixes(filename, loci):
        found_with_chr = False
        found_without_chr = False
        # Iterate through the loci and check for reads in each locus,
        # with and without the 'chr' prefix, until we get a match for one.
        # Using a 5-base buffer to allow better matching due to potentially
        # misspecified genomic coordinates.

        source = pysam.AlignmentFile(filename, 'rb')  

        for locus in loci:
            with_chr = locus.chromosome
            if with_chr[0:3] != 'chr':
                with_chr = 'chr{0}'.format(with_chr)
            without_chr = with_chr[3:]


            # Make sure the start coordinate isn't below 1
            start_pos = locus.start - 5
            if start_pos < 1:
                start_pos = 1
            end_pos = locus.end + 5

            # First check if we can find one WITH the 'chr' prefix.
            try:
                reads = source.fetch(with_chr, start_pos, end_pos)
                # If we proceed beyond this line, we found matches
                found_with_chr = True
                break
            except ValueError:
                found_with_chr = False

            if not found_with_chr:
                # Check if we can find reads WITHOUT the 'chr' prefix.
                try:
                    reads = source.fetch(without_chr, start_pos, end_pos)
                    # If we proceed beyond this line, we found matches
                    found_without_chr = True
                    break
                except ValueError:
                    found_without_chr = False

            if found_with_chr or found_without_chr:
                break

        if found_with_chr:
            return True
        elif found_without_chr:
            return False
        else:
            # No matches for anything!
            return None
        # end .bam_uses_chr_prefixes()


    """
    Producer thread that reads data from the BAM file and populates the queue
    with matched reads using pysam; uses semaphores to prevent loss of 
    synchronization state or overwhelming the read buffer.
    """
    def extract_reads(self, filename, loci, full, empty, mutex, queue, consumers):
        use_chr_prefix = KmerRepeatCounter.bam_uses_chr_prefixes(filename, loci)
        if use_chr_prefix is None:
            tprint('Fatal error! Could not find any matched reads.')
            exit(1)

        if self.debug_output:
            tprint('Extractor> Thread starting for {0}'.format(filename))
        source = pysam.AlignmentFile(filename, 'rb')   

        available_chromosomes = set([hash(str(x)) for x in source.references])
        for locus in loci:

            # Format the chromosome to be compatible with how the reads are
            # stored in the BAM file (i.e. with or without 'chr' prefix).
            chromosome = locus.chromosome
            if use_chr_prefix and chromosome[0:3] != 'chr':
                # Prepend the 'chr' prefix.
                chromosome = 'chr{0}'.format(chromosome)
            elif not use_chr_prefix and chromosome[0:3] == 'chr':
                # Remove the 'chr' prefix.
                chromosome = chromosome[3:]

            if hash(chromosome) in available_chromosomes:

                # Make sure the start coordinate isn't below 1
                start_pos = locus.start - 5
                if start_pos < 1:
                    start_pos = 1
                end_pos = locus.end + 5


                for read in source.fetch(chromosome, start_pos, end_pos):
                    # Use AlignedSegment object to create a list, which is
                    # then used in the creation of the SAMRead object,
                    # since the AlignedSegment objects are C-structs and cannot
                    # be passed to the consumer threads.
                    data = [
                        read.query_name,
                        read.flag,
                        chromosome,
                        read.reference_start,
                        read.mapping_quality,
                        read.cigarstring,
                        '', 
                        '',
                        '',
                        read.query_sequence,
                        KmerRepeatCounter.quality_scores_to_symbols(read.query_qualities)
                    ]
                    # CIGAR of None means it was likely an asterisk (*), so the read
                    # will get ignored since something was wrong with the alignment.
                    if read.cigarstring is not None:
                        read = SAMRead('\t'.join([str(x) for x in data]))
                        item = [locus, read]

                        # Use semaphores to handle proper writing into the queue.
                        empty.acquire()
                        mutex.acquire()
                        queue.put(item)
                        mutex.release()
                        full.release()

        source.close()

        if self.debug_output:
            tprint('Extractor> Extracted all reads for all target loci.')
        # Add set amount of end signals to queue to end consumers
        for i in range(0, consumers * 2):
            empty.acquire()
            mutex.acquire()
            queue.put([False,False])
            mutex.release()
            full.release()
        if self.debug_output:
            tprint('Extractor> Queued up {0} termination signals.'.format(consumers))
        return True
        # end .extract_reads()



    """
    Consumer thread that reads data put on the queue by the BAM reading thread; 
    processes each read to get the repeat count for the repeat unit (k-mer).
    """
    def read_analyzer(self, queue_in, queue_out, full, empty, mutex_in, mutex_out, queue_full):
        query_delay = 0.100 # in seconds
        n = 0
        if self.debug_output:
            tprint('Analyzer> Thread {0} started.'.format(os.getpid()))
        while True:
            n += 1
            if self.debug_output and (n % 10000 is 0):
                tprint('Analyzer> Thread {0}'.format(os.getpid()) +
                    ' still ALIVE, loop {0}'.format(n))


            full.acquire()
            mutex_in.acquire()
            item = queue_in.get()
            locus = item[0]
            read = item[1]
            mutex_in.release()
            empty.release()   

            if not locus:
                # Element was set to FALSE; signals thread termination.
                if self.debug_output:
                    tprint('Analyzer> Thread {0}'.format(os.getpid()) +
                    ' received termination signal.')
                break

            if self.passes_qc_filter(read):
                # Get the repeat count for the repeat unit, along with the
                # internal offset (how many bases into the read the locus
                # actually starts at).
                repeat_count, repeat_start = self.locus_repeat_count(read, locus)
                    
                # A repeat_count of -1 would mean the read did not contain
                # the target locus, or the locus position within the read
                # could not be determined.
                if repeat_count >= 0:
                    # Figure out where the repeat region ends (offset-wise).
                    repeat_end = repeat_start + (repeat_count * locus.kmer_length)

                    # Make sure the locus itself has sufficient quality.
                    if self.passes_locus_qc_filter(read, repeat_start, repeat_end):
                        # Passed locus QC filter; read is processed. Acquire semaphore
                        # to return repeat count as output.
                        queue_full.acquire()
                        mutex_out.acquire()
                        queue_out.put([locus.locus(), repeat_count])
                        mutex_out.release()

        if self.debug_output:
            tprint('Analyzer> Ending thread {0}'.format(os.getpid()))
        return True   
        # end .read_analyzer()

    """
    Performs a status check on the producer (extractor), consumer (analyzer),
    and shared queue structures to make sure the program hasn't lost proper
    multiprocessing functionality.
    """
    def status_check(self, qsize):
        ok = '\033[92m'
        reset = '\033[0m'
        fail = '\033[91m'
        if qsize is 0:
            queue_status = 'EMPTY'
        else:
            queue_status = '{0} ITEMS'.format(qsize)

        if self.has_live_producer() and not self.has_live_consumers():
                tprint(fail + 'Main> Analyzer(s) LIVE, Extractor DEAD, Queue ' + 
                    queue_status + '. ' + reset)
                tprint(fail + 'Main> Teriminating process due to ' + 
                    'multiprocessing failure.' + reset)
                exit(1)

        # end .status_check()    

    """
    Inspects the reads, finds corresponding MSI loci,
    and returns data on the locus size/length distributions.
    """
    def process(self, input_filepath, msi_loci, config):
        self.__reset()

        # Generate dictionary read counts for loci
        counts = {}
        loci = []
        for locus in msi_loci:
            counts[locus.locus()] = {}
            loci.append(locus)

        # Generate input and output queues and (mutex) semaphores
        # for each.
        queue_out = Queue()
        queue_in = Queue()
        queue_full = BoundedSemaphore(100)
        full = Semaphore(0)
        empty = BoundedSemaphore(40)
        mutex_out = Semaphore(1)
        mutex_in = Semaphore(1)

        # Set amount of consumer threads; minimum one.
        consumer_threads = config['threads'] - 1
        if consumer_threads < 1:
            consumer_threads = 1

        # Create producer thread; currently only using single thread
        # since I/O is more of the limiter than CPU bound processes.
        self.__producer = Process(target=self.extract_reads, args=(
            input_filepath, 
            msi_loci, 
            full, 
            empty, 
            mutex_out, 
            queue_in, 
            consumer_threads))
        self.__producer.start()


        # Spawn the set amount of threads/processes
        if self.debug_output:
            tprint('Main> Generating {0} analyzer process(es).'.format(consumer_threads))
        for i in range(0, consumer_threads):
            p = Process(target=self.read_analyzer, args=(
                queue_in, 
                queue_out, 
                full, 
                empty, 
                mutex_in, 
                mutex_out,
                queue_full))
            self.__consumers.append(p)
            self.__consumers[-1].start()

        # Iterate through the loci, fetching any reads and pushing them to 
        # the pool of threads, collecting the output as they process it.
        query_delay = 0.050 # In seconds
       
        loop_counter = 0
        proc_check_interval = 100
        while (not queue_out.empty() or self.has_live_threads()):
            # Sleep for the set amount of time so the queue isn't constantly
            # getting hammered with queries
            time.sleep(query_delay)
            loop_counter += 1
            if loop_counter % proc_check_interval is 0:
                # Time to check that the consumers
                # didn't die while the producer is still producing
                mutex_out.acquire()
                self.status_check(queue_out.qsize())
                mutex_out.release()

            while not queue_out.empty():
                # There is data on the queue to be processed;
                # the return from the queue should be a tuple
                # with (locus, repeat_count)
                mutex_out.acquire()
                result = queue_out.get()
                locus = result[0]
                repeat_count = result[1]
                if repeat_count >= 0:
                    if locus not in counts:
                        counts[locus] = {}
                    if repeat_count not in counts[locus]:
                        counts[locus][repeat_count] = 0
                    counts[locus][repeat_count] += 1
                mutex_out.release()
                queue_full.release()

            if not self.has_live_threads():
                # All processes should have terminated.
                if self.debug_output:
                    tprint('Main> All processes complete.')
                break
        # end while loop

        return counts
        # end .process()
    
    # end KmerRepeatCounter class definition.

 




"""
Use arguments from the command line parameters (input: args)
to generate a config setting dict.
"""
def generate_config(args):
    config = {}

    config['genome'] = os.path.abspath(args.genome)
    if not os.path.isfile(config['genome']):
        tprint('Error: {0} does not exist!'.format(config['genome']))
        exit(1)

    config['threads'] = int(args.threads)
    if config['threads'] < 1:
        tprint('Error: Cannot specify less than one thread. '
            + '(Provided {0}).'.format(config['threads']))
        exit(1)

    config['bedfile'] = os.path.abspath(args.bedfile)
    if args.bedfile is None:
        tprint('Error: BED file not provided!')
        exit(1)
    else:
        bedfile = os.path.abspath(args.bedfile)
        if not os.path.isfile(bedfile):
            tprint('Error: {0} does not exist!'.format(bedfile))
            exit(1)

    if args.normal is None:
        tprint('Error: Normal BAM/SAM file not provided!')
        exit(1)
    else:
        config['normal_filepath'] = os.path.abspath(args.normal)
        if not os.path.isfile(config['normal_filepath']):
            tprint('Error: {0} does not exist!'.format(config['normal_filepath']))
            exit(1)
    
    if args.tumor is None:
        tprint('Error: Tumor BAM/SAM file not provided!')
        exit(1)
    else:
        config['tumor_filepath'] = os.path.abspath(args.tumor)
        if not os.path.isfile(config['tumor_filepath']):
            tprint('Error: {0} does not exist!'.format(config['tumor_filepath']))
            exit(1)


    if args.output is None:
        tprint('Error: Output filepath must be specified!')
    else:
        config['output_filepath'] = os.path.abspath(args.output)
        # Make sure output folder exists
        output_dir = os.path.dirname(config['output_filepath'])
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)


    config['min_read_quality'] = float(args.mrq)
    config['min_read_length'] = int(args.mrl)
    config['min_locus_quality'] = float(args.mlq)
    config['debug_output'] = args.debug_output
    return config
    # end .generate_config()

strip_chr_re = re.compile(r'^chr')
#Helper method for ordering loci
def cmp_loci(x, y):
    def parse_locus(l):
        pieces = l.split(':')
        chr = strip_chr_re.sub("", pieces[0])
        pieces = pieces[1].split('-')
        start = int(pieces[0])
        end = int(pieces[1])
        return (chr, start, end)
    x_chr, x_start, x_end = parse_locus(x)
    y_chr, y_start, y_end = parse_locus(y)
    if (x_chr.isdigit() and not y_chr.isdigit()):
        return -1
    elif (y_chr.isdigit() and not x_chr.isdigit()):
        return 1
    elif x_chr.isdigit():
        x_chr = int(x_chr)
        y_chr = int(y_chr)
    if (x_chr < y_chr):
        return -1
    elif (x_chr > y_chr):
        return 1
    if (x_start < y_start):
        return -1
    elif (x_start > y_start):
        return 1
    elif (x_end < y_end):
        return -1
    elif (x_end > y_end):
        return 1
    else:
        return 0

"""
Write the output of the program into the target file.
The output consists of tab-separated values, with each line having:
- locus
- k-value
- supporting reads for k-value in normal BAM
- supporting reads for k-value in tumor BAM
"""
def write_output(filepath, normal, tumor):
    fileout = open(filepath, 'w')

    loci = set(normal.keys())
    loci = loci.union(tumor.keys())
    if (sys.version_info > (3, 0)):
        loci = sorted(loci, key=functools.cmp_to_key(cmp_loci))
    else:
        loci = sorted(loci, cmp=cmp_loci)
    header = '\t'.join(['Locus','Repeats','Normal','Tumor'])
    fileout.write(header + '\n')
    for locus in loci:
        locus_output = generate_locus_output(normal[locus], tumor[locus])
        for line in locus_output:
            # Prepend locus
            line = [locus] + line
            line = '\t'.join([str(x) for x in line])
            fileout.write(line + '\n')
    fileout.close()
    # end .write_output()


"""
Helper method for generating the output lines to be written.
"""
def generate_locus_output(normal, tumor):
    output = []
    # Get repeat lengths to use as keys for accessing
    # both the normal and tumor data.
    lengths = set(normal.keys())
    for length in tumor.keys():
        lengths.add(length)
    lengths = sorted(lengths)
    for length in lengths:
        n_repeats = 0
        t_repeats = 0
        if length in normal:
            n_repeats = normal[length]
        if length in tumor:
            t_repeats = tumor[length]
        output.append([length, n_repeats, t_repeats])
    return output
    # end .generate_locus_output()

"""
Uses PySAM to generate index for BAM file if one doesn't exist.
"""
def generate_index_if_needed(filepath):
    index_file = os.path.abspath(filepath) + '.bai'
    if not os.path.isfile(index_file) and not os.path.isfile(os.path.abspath(filepath)[:-4] + '.bai'):
        # Index file doesn't exist; generate it
        pysam.index(filepath, index_file)
    return True
    # end .generate_index_if_needed()

if __name__ == "__main__":
    prog_name = 'MANTIS K-Mer Repeat Counter'
    tprint(prog_name)

    parser = argparse.ArgumentParser(description=prog_name)

    parser.add_argument('-n', '--normal', dest='normal', type=str, 
        required=True, help='Normal input (SAM/BAM) file.')

    parser.add_argument('-t', '--tumor', dest='tumor', type=str, 
        required=True, help='Tumor input (SAM/BAM) file.')

    parser.add_argument('-b', '--bedfile', dest='bedfile', type=str,
        required=True, help='Input BED file.')

    parser.add_argument('-o', '--output', dest='output', type=str, 
        help='Output BED filename.')

    parser.add_argument('-mrq', '--min-read-quality', dest='mrq',
        type=float, default=10.0, help='Minimum average read quality.')

    parser.add_argument('-mlq', '--min-locus-quality', dest='mlq',
        type=float, default=20.0, help='Minimum average locus quality.')

    parser.add_argument('-mrl', '--min-read-length', dest='mrl',
        type=int, default=35, 
        help='Minimum (unclipped/unmasked) read length.')

    parser.add_argument('--threads', dest='threads', type=int,
        default=1, help='Number of threads/proceseses to use.')

    parser.add_argument('--genome', dest='genome', type=str,
        help='Path to reference genome (FASTA).')

    parser.add_argument('--debug', dest='debug_output', action='store_true',
        help='Print debug output from multithreading.')

    args = parser.parse_args()

    # Use helper method for checking args and creating config dict.
    config = generate_config(args)

    # Load MSI loci from the BED line.
    tprint('Loading target MSI loci from BED file ...')
    mll = MSILocusLoader(config)
    msi_loci = mll.load_loci(config['bedfile'])
    tprint('Loaded {0} loci.'.format(len(msi_loci)))

    # Generate instance of the Kmer Repeat Counter
    krc = KmerRepeatCounter(config)

    tprint('Processing normal input file with ' +
        '{0} thread(s) ...'.format(config['threads']))
    generate_index_if_needed(config['normal_filepath'])
    normal = krc.process(config['normal_filepath'], msi_loci, config)
    tprint('Done processing normal.')

    krc = KmerRepeatCounter(config)
    tprint('Processing tumor input file with ' +
        '{0} thread(s) ...'.format(config['threads']))
    generate_index_if_needed(config['tumor_filepath'])    
    tumor = krc.process(config['tumor_filepath'], msi_loci, config)  
    tprint('Done processing tumor.')

    write_output(config['output_filepath'], normal, tumor)
    print('Saved output in {0}'.format(config['output_filepath']))
    # Done.
    exit(0)
