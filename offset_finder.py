# @file offset_finder.py
# @author Esko Kautto (esko.kautto@osumc.edu)
# @updated 2016-06-16

from copy import deepcopy
from structures import CIGAR
import time

class OffsetFinder:
    max_iterations = 10

    @staticmethod
    def message(msg):
        if type(msg) is str:
            msg = msg.split('\n')

        for line in msg:
            timestamp = time.strftime('%m/%d/%y %H:%M:%S')        
            print('[{0}] OffsetFinder: {1}'.format(timestamp, line))
    # end .message()

    """
    Determine the offset for starting to read the read's
    sequence, based on the read's start position, the
    CIGAR string, and the msi locus start. Since there
    could be potential insertions or deletions before
    the MSI locus, the offset has to take those into account
    from the CIGAR.
    """
    @staticmethod
    def find_offset(read, locus):
        # First, get the best-effort estimate for the offset,
        # based on the read's expected genomic start position,
        # the locus start, and the read's CIGAR string.
        offset = OffsetFinder.estimate_offset(read, locus)

        # Check the boundary to make sure the entire read isn't
        # slightly offset; this can happen if the aligner shifts
        # the read.pos when an indel is very close to the edge of
        # the read sequence.
        offset_change = OffsetFinder.adjust_for_shifts(
            read,
            locus,
            offset)

        offset += offset_change
        return offset
        # end .find_offset()

    """
    Estimate the offset of the read sequence in relation to the 
    target locus start. The returned offset will be an integer in
    the range of [0, |sequence|]. An offset of 0 will indicate the
    read starts right at the locus boundary; an offset of |sequence|
    indicates the read doesn't cover the locus.
    """
    @staticmethod 
    def estimate_offset(read, locus):
        offset = 0
        genomic_pos = read.pos
        cigar = CIGAR.to_array(read.cigar)

        if cigar[0][0] == 'S':
            # Soft clipped. Pop the element off the list so that the soft
            # clipped region isn't double-processed later in the method. 
            # Ideally the SAMRead preprocessing will have removed any
            # preceding soft clipping so this won't be needed.
            offset += cigar[0][1]
            cigar.pop(0)                            

        if genomic_pos == locus.start:
            # Genomic position of the start of the read is the start
            # of the MSI locus. Alignment might have resulted in position
            # being slightly off; make sure the sequence starts with the 
            # target kmmer. If it doesn't, shift to the right (increment
            # offset) until the right kmer is encountered.
            # No (further) offset required.
            if not OffsetFinder.at_target_kmer(read, locus, offset):
                offset += OffsetFinder.adjust_reading_frame(read, locus, offset)
            return offset

        # Iterate through the cigar string until we find the
        # start of the locus
        while cigar:
            token = cigar[0][0]
            size = cigar[0][1]
            cigar.pop(0) # Remove this CIGAR segment.

            if token == 'S':
                # Soft clipped segment. This should not happen in the middle 
                # of a CIGAR string; since the SAMRead preprocessor gets rid
                # of prepended soft clippings, and there is a previous check
                # for it within this method, the only time this should be 
                # encountered is at the end of a sequence, at which point the
                # read likely did not cover the MSI locus.
                genomic_pos += size
                offset += size

            elif token == 'D':
                # A deletion in the sequence. The genomic position will be 
                # incremented by the size of the deletion, but the offset
                # will not be affected since the deleted bases are not in
                # the read sequence.
                genomic_pos += size

            elif token == 'I':
                # An insertion of one or more bases into the sequence. 
                # Since the insertion is not part of the original (reference)
                # sequence, the genomic position will not be incremented, but
                # the offset will be to essentially skip over the insertion.
                offset += size

            elif token == 'M':
                if genomic_pos + size > locus.start:
                    # Matched sequence crosses locus boundary; increment
                    # the offset only be the difference so that bases are
                    # not accidentally skipped.
                    offset += locus.start - genomic_pos
                    genomic_pos = locus.start
                
                elif genomic_pos + size == locus.start:
                    # Matched sequence reaches the exact boundary of the
                    # MSI locus. If the next CIGAR string segment is an 
                    # indel, run an additional calculation to check if
                    # further offset is necessary.
                    offset += size
                    genomic_pos = locus.start

                    if len(cigar) and cigar[0][0] in ['D', 'I']:
                        offset += OffsetFinder.boundary_indel_offset(
                            read, locus, cigar, offset, genomic_pos)


                elif genomic_pos + size == locus.start - 1:
                    # Right before the locus boundary. As with the previous
                    # case, make sure an indel in the next CIGAR segment
                    # is accounted for if necessary.
                    offset += size
                    genomic_pos = locus.start - 1
                    if len(cigar) and cigar[0][0] in ['D', 'I']:
                        offset += OffsetFinder.boundary_indel_offset(
                            read, locus, cigar, offset, genomic_pos)
                else:
                    # Not at the locus boundary yet; increment offset and
                    # genomic position normally.
                    offset += size

                    genomic_pos += size


            if genomic_pos >= locus.start:
                # MSI locus reached
                break


        # Return the estimate offset. 
        return offset

        # end .estimate_offset()


    """
    Handles cases where an indel exists right at the locus boundary.
    This is relatively likely if the aligner accounted for any repeat 
    count changes by placing the added or removed k-mers at the start
    of the MSI locus.
    """
    @staticmethod
    def boundary_indel_offset(read, locus, cigar, offset, genomic_pos):
        offset_change = 0
        if genomic_pos == locus.start and False:
            # Only checks insertions at the locus boundary.
            insertion_seq = read.seq[seq_pos:seq_pos+insert_size]
            if insertion_seq[0:locus.kmer_length] != locus.kmer:
                # The insertion appers to be a random insertion, not
                # part of the repeat count change for the MSI locus.
                # Shift the offset accordingly.
                offset += insert_size
        return offset_change
        # end .boundary_indel_offset()

    """
    Checks the boundary position to make sure the entire read before the MSI locus
    isn't slightly shifted; this could happen if the aligner shifts the 'start' 
    position of the entire read when the locus is close to the edge of the read
    sequence (e.g. the start position for caTTTTT.. might get decremented by N
    if there are N extra T's in the locus, despite the read starting with 'ca'.)
    """
    @staticmethod
    def adjust_for_shifts(read, locus, offset):
        # First try to get sequence into the proper reading frame (if it
        # isn't already) by left- or right-shifting it.

        offset_change = OffsetFinder.adjust_reading_frame(read, locus, offset)
        
        if offset_change < 0:
            # Encountered problem adjusting reading frame.
            return 0

        # Check if we can adjust the offset to an earlier position to
        # capture an earlier (continuous) start for the repeat sequence.
        offset_change += OffsetFinder.shift_start_left(read, locus, offset + offset_change)
        return offset_change
        # end .adjust_for_shifts()

    """
    Tries to adjust the read sequence into the proper reading frame.
    """
    @staticmethod
    def adjust_reading_frame(read, locus, offset):
        if OffsetFinder.at_target_kmer(read, locus, offset):
            # Already in proper reading frame, no adjusting required.
            return 0

        # Try to adjust the reading frame until the read is within
        # the correct frame.
        offset_change = 0
        in_frame = False
        for i in range(1, OffsetFinder.max_iterations + 1):
            if OffsetFinder.at_target_kmer(read, locus, offset - i):
                # Left-shifting read moved it into the proper reading frame.
                in_frame = True
                offset_change = -1 * i
            elif OffsetFinder.at_target_kmer(read, locus, offset + i):
                # Right-shifting read moved it into the proper reading frame.
                in_frame = True
                offset_change = i
            if in_frame:
                # Read is now in the proper reading frame.
                break

        if not in_frame:
            # Couldn't get the read into the right reading frame.
            offset_change = -1
        return offset_change
        # end of .adjust_reading_frame()


    """
    Tries to find an 'earlier' start for the repeat sequence
    by checking left (before) of the current offset. Assumes
    the read is already in the proper reading frame.
    """
    @staticmethod
    def shift_start_left(read, locus, offset = 0):
        # Start by assuming no change will occur.
        offset_change = 0

        # Kmer length acts as the shift unit (-/neg since we're shifting left).
        shift_unit = -1 * locus.kmer_length

        for i in range(1, OffsetFinder.max_iterations + 1):
            start = offset + (i * shift_unit)
            stop = offset + (i * shift_unit) + locus.kmer_length
            if start < 0 or stop > read.length:
                # Going out of boundaries on the read; stop.
                break

            if read.seq[start:stop] == locus.kmer:
                # Found the kmer at the earlier offset position
                offset_change = (i * shift_unit)
            else:
                # No longer within repeat sequence
                break

        return offset_change
        # end of .shift_start_left()


    """
    Basic helper method that checks if the current offset in the read
    sequence matches the targeted k-mer sequence.
    """
    @staticmethod
    def at_target_kmer(read, locus, offset = 0):
        return (read.seq[offset:offset+locus.kmer_length] == locus.kmer)
        # end .at_target_kmer()

    # end OffsetFinder class definition.

