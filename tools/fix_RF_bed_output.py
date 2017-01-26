# @file fix_RF_bed_output.py
# @author Esko Kautto (esko.kautto@osumc.edu)
# @updated 2017-01-26
#
# The script reads in the 6-column BED file output generated
# by the RepeatFinder program and parses through it to fix
# any incorrect annotations that might have been created due
# to unexpected formats of data in the reference genome.
# An example are the alternate haplotype reference ("chromosome")
# entries from HG38, which can result in the reference spanning
# two columns, causing the 6-column output to incorrectly
# become 7 columns and cause MANTIS to reject the input BED file.

import os
import sys
import argparse

"""
Attempts to determine if the input string follows the correct
format for the motif repeat column in the bedfile, i.e.
(A|C|T|G)XX
"""
def valid_motif_string(string):
	valid = True
	if string[0] == '(' and ')' in string:
		check = string[1:].split(')')
		if len(check) != 2:
			# Split improperly on parenthesis.
			valid = False
		else:
			if len(check[0]) < 1:
				# No motif specified
				valid = False
			else:
				for c in check[0].upper():
					if c not in ['A', 'T', 'C', 'G']:
						# Invalid base/character
						valid = False
						break
				if not check[1].isdigit():
					# Improper repeat count.
					valid = False
	else:
		valid = False
	return valid
	# end valid_motif_string()


"""
Checks if the input line appears to follow the expected
6-column BED format being used. If the input appears to
match the format, True is returned; if not, False.
"""
def check_line_validity(line):
	valid = True
	if len(line) != 6:
		# Wrong number of columns
		valid = False
	else:
		if len(line[0].strip()) is 0:
			# Reference/chromosome appears empty.
			valid = False
		elif not line[1].isdigit() or not line[2].isdigit():
			# The start and/or stop coordinates appear incorrect.
			valid = False
		else:
			# We only need to check the 4th column, since the
			# values 5th and 6th columns are ignored at this time.
			if not valid_motif_string(line[3]):
				valid = False

	return valid
	# end check_line_validity()

"""
Attempts to fix the line values. If successful,
returns the line as a list of the corrected values. If 
not, returns False to signify failure.
"""
def fix_line(line):
	output = []
	if len(line) is 6:
		# Not attempting to fix, return unchanged.
		output = line
	elif len(line) < 6:
		# Insufficient fields.
		print('Cannot fix line "{0}": insufficient columns'.format(
			'\t'.join(line)))
		output = False
	else:
		proceed = True
		if line[-1] in ['+', '-']:
			output.append(line[-1])
		else:
			proceed = False

		if proceed and line[-2].isdigit():
			output.append(line[-2])
		else:
			proceed = False

		if proceed and valid_motif_string(line[-3]):
			output.append(line[-3])
		else:
			proceed = False

		if proceed and line[-4].isdigit():
			output.append(line[-4])
		else:
			proceed = False

		if proceed and line[-5].isdigit():
			output.append(line[-5])
		else:
			proceed = False

		if proceed:
			output.append(line[0])

		if not proceed:
			output = False
		else:
			output = output[::-1]
			print('Fixed line "{0}" => "{1}"'.format(
				'\t'.join(line),
				'\t'.join(output),
				))

	return output
	# end fix_line()

if __name__ == '__main__':
	program_name = 'MANTIS RepeatFinder BED file fixer'
	parser = argparse.ArgumentParser(description=program_name)
	parser.add_argument('-i', '--input', dest='input', type=str, 
	    required=True, help='Input filepath.')
	parser.add_argument('-o', '--output', dest='output', type=str, 
	    required=True, help='Output filepath.')
	parser.add_argument('--overwrite', dest='overwrite', action='store_true')
	args = parser.parse_args()

	input_filepath = os.path.abspath(args.input)
	if not os.path.isfile(input_filepath):
		print('Error: Input filepath {0} does not exist!'.format(input_filepath))
		exit(1)

	output_filepath = os.path.abspath(args.output)
	if os.path.isfile(output_filepath):
		if args.overwrite:
			print('Warning: Output file is being overwritten.')
		else:
			print('Error: Output file {0} already exists!'.format(output_filepath))
			print('To overwrite the output, use the --overwrite argument.')
			exit(1)
	else:
		output_folder = os.path.dirname(output_filepath)
		if not os.path.isdir(output_folder):
			os.makedirs(output_folder)

	with open(output_filepath, 'w') as fileout:
		with open(input_filepath, 'Ur') as filein:
			for line in filein:
				line = line.strip()
				if len(line):
					line = line.split('\t')
					if check_line_validity(line):
						# Nothing to fix.
						pass
					else:
						# Attempts to fix line data.
						line = fix_line(line)

					if line:
						fileout.write('{0}\n'.format('\t'.join(line)))
			filein.close()
		fileout.close()
	# done with __main__


