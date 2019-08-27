#!/usr/bin/env python
# -*- coding: utf-8 -*-

# @file structures.py
# @author Esko Kautto (esko.kautto@osumc.edu)
# @updated 2016-06-16

from copy import deepcopy
import time

"""
Container class that holds information from SAM/BAM reads in a format
that is useful for MANTIS.
"""
class SAMRead(object):
    def __init__(self,line):
        if type(line) is str:
            line = line.strip().split('\t')
            self.qname = line[0]
            self.chromosome = line[2]
            self.pos = int(line[3])
            self.mapq = int(line[4])
            self.cigar = line[5]
            self.seq = line[9]
            self.qual = line[10]
            self.length = self.calculate_length()
            self.start = self.pos
            self.end = self.start + CIGAR.length(self.cigar)
            self.preprocess_read()
            self.quality = self.calculate_quality()
        # end .__init__()

    def locus(self):
        return '{0}:{1}-{2}'.format(self.chromosome, self.start, self.end)
        # end .locus()

    def preprocess_read(self):
        # Trims away soft clipping
        seq_length = 0
        gen_length = 0
        cigar = CIGAR.to_array(self.cigar)
        trimmed_cigar = []
        trimmed_seq = []
        trimmed_qual = []
        index = 0
        for i in range(0, len(cigar)):
            token = cigar[i][0]
            size = cigar[i][1]
            if token not in ['S', 'H']:
                trimmed_seq.append(self.seq[index:index+size])
                trimmed_qual.append(self.qual[index:index+size])
                seq_length += size
                gen_length += size
                trimmed_cigar.append([token, size])

            if token != 'H':
                # Hard-clipped segments are not included in the sequence;
                # only increment index for non-hard clipped section.
                index += size
        # Update self with updated values
        self.length = seq_length
        self.seq = ''.join(trimmed_seq)
        self.qual = ''.join(trimmed_qual)

        # Reassemble CIGAR
        c = []
        for i in range(0, len(trimmed_cigar)):
            c.append(str(trimmed_cigar[i][1]) + trimmed_cigar[i][0])
        self.cigar = ''.join(c)

        # Update end postion based on estimated genomic length of read
        self.end = self.start + gen_length
    # end .preprocess_read()



    def calculate_quality(self, offset=0, length=0):
        if length < 1:
            length = len(self.qual)

        if length is 0:
            return 0.0


        total_score = 0
        for i in range(offset, offset+length):
            score = SAMRead.base_score(self.qual[i])
            total_score += score
        average = (1.0 * total_score) / length
        return average
        # end .calculate_quality()

    def calculate_length(self):
        cigar = CIGAR.to_array(self.cigar)
        if cigar == '*':
            return 0
        length = len(self.seq)
        while cigar:
            segment = cigar.pop(0)
            token = segment[0]
            size = segment[1]
            if token is 'S':
                length -= size
        return length
        # end .calculate_length()

    @staticmethod
    def base_score(base):
        values = '!"#$%&' + "'" + '()*+,-./0123456789:;<=>?@ABCDEFGHIJK'
        score = values.find(base)
        return score
        # end .base_score()

    # end SAMRead class definition.


"""
Helper class for handling CIGAR strings.
"""
class CIGAR(object):
    # Converts the CIGAR string to a token+length array
    @staticmethod
    def to_array(string):
        if string == '*':
            return []
        tokens = ['M','S','N','I','D','H']
        for token in tokens:
            string = string.replace(token, token + '.')
        string = filter(None, string.split('.'))

        array = []

        for substring in string:
            # Last character is the token, all the other
            # characters are the length.
            token = substring[-1]
            value = int(substring[:-1])
            array.append([token, value])
        return array
        # end .to_array()

    # Estimates the length of the read from the CIGAR.
    @staticmethod
    def length(string, count_insertions = True):
        if type(string) is str:
            array = CIGAR.to_array(string)

        length = 0
        for segment in array:
            if segment[0] in ['M', 'D', 'S']:
                length += segment[1]
            elif segment[0] == 'I' and count_insertions:
                length += segment[1]

        return length
        # end of .length()

    @staticmethod
    def show_sequence(cigar, sequence):
        if type(cigar) is not list:
            cigar = CIGAR.to_array(cigar)

        output = []
        index = 0
        for segment in cigar:
            size = segment[1]
            if segment[0] != 'D':
                subsequence = sequence[index:index+size]
                index += size
            else:
                subsequence = '?' * size
            print('{0}{1}\t{2}'.format(segment[0], size, subsequence))
        # end .show_sequence()
    # end CIGAR class definition.

"""
Class that represents a genomic locus.
"""
class Locus(object):
    locus_number = 1
    def __init__(self, line):
        line = line.strip().split('\t')
        self.chromosome = line[0]
        self.start = int(line[1]) + 1
        self.end = int(line[2])
        self.name = 'Locus #{0}'.format(self.locus_number)
        Locus.locus_number += 1
        # end .__init__()

    def __lt__(self, other):
        if self.chromosome == other.chromosome:
            if self.start != other.start:
                return self.start < other.start
            return self.end < other.end

        return self.chromosome < other.chromosome
        # end .__lt__()

    def locus(self):
        locus = '{0}:{1}-{2}'.format(
            self.chromosome,
            self.start,
            self.end)
        return locus
        # end .locus()

    def get_hash(self):
        h = []
        h.append(self.locus())
        if self.kmer is not None:
            k = '({0}){1}'.format(self.kmer, self.repeats)
            h.append(k)
        return ':'.join(h)
        # end .get_hash()

    def __str__(self):
        out = '{0} covering {1}'.format(
            self.__class__.__name__,
            self.locus())
        return out
        # end .__str__()
    # end Locus class definition.

"""
Class that expands the base Locus class to include
other data needed by an MSI locus.
"""
class MSILocus(Locus):
    def __init__(self, line):
        line = line.strip().split('\t')
        self.chromosome = line[0]
        self.start = int(line[1]) + 1
        self.end = int(line[2])
        if len(line) > 3:
            if '(' in line[3]:
                k_r = line[3].split(')',2)
                self.kmer = k_r[0].split('(',2)[1]
                self.kmer_length = len(self.kmer)
                self.repeats = int(k_r[1])
            self.name = line[3]
            self.strand = line[5]
        # end .__init__()

    # end MSILocus class definition.
