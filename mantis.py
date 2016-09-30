# @file mantis.py
# @author Esko Kautto (esko.kautto@osumc.edu)
# @updated 2016-07-05

import os
import sys
import argparse
import operator
import subprocess
from helpers import iteritems, check_bedfile_format, required_modules_present
from defaults import load_settings, Parameter



if __name__ == "__main__":
    prog_name = 'Microsatellite Analysis for Normal-Tumor InStability (v1.0.3)'
    print(prog_name)

    # Make sure Pysam and NumPy are available in environment
    required_modules_present(['Pysam', 'NumPy'])


    parser = argparse.ArgumentParser(description=prog_name)

    parser.add_argument('-cfg', '--config', dest='cfg', type=str,
        help='Configuration file.')

    parser.add_argument('-n', '--normal', dest='normal', type=str, required=True,
        help='Normal input BAM file.')

    parser.add_argument('-t', '--tumor', dest='tumor', type=str, required=True,
        help='Tumor input BAM file.')

    parser.add_argument('--threads', dest='threads', type=int,
        help='How many threads (processes) to use.')

    parser.add_argument('-b', '--bedfile', dest='bedfile', type=str, 
        required=True, help='Input BED file.')

    parser.add_argument('-o', '--output', dest='output', type=str, 
        required=True, help='Output filename.')

    parser.add_argument('-mrq', '--min-read-quality', dest='mrq', type=float, 
        help='Minimum average read quality.')

    parser.add_argument('-mlq', '--min-locus-quality', dest='mlq', type=float,
        help='Minimum average locus quality.')

    parser.add_argument('-mrl', '--min-read-length', dest='mrl', type=int, 
        help='Minimum (unclipped/unmasked) read length.')

    parser.add_argument('-mlc', '--min-locus-coverage', dest='mlc', type=int, 
        help='Minimum coverage required for each of the normal ' +
        'and tumor results.')

    parser.add_argument('-mrr', '--min-repeat-reads', dest='mrr', type=int, 
        help='Minimum reads supporting a specific repeat count.')

    parser.add_argument('-sd', '--standard-deviations', dest='sd', type=float,
        help='Standard deviations from mean before repeat count is ' +
        'considered an outlier')

    parser.add_argument('--genome', dest='genome', type=str,
        help='Path to reference genome (FASTA).')

    parser.add_argument('--difference-threshold', dest='dif_threshold', type=float,
        help='Default difference threshold value for calling a sample unstable.')

    parser.add_argument('--distance-threshold', dest='euc_threshold', type=float,
        help='Default distance threshold value for calling a sample unstable.')

    parser.add_argument('--dissimilarity-threshold', dest='cos_threshold', type=float,
        help='Default dissimilarity threshold value for calling a sample unstable.')

    args = parser.parse_args()


    if args.cfg is not None:
        config_file_path = os.path.abspath(args.cfg)
        if not os.path.isfile(config_file_path):
            print('Error! Config file {0} not found!'.format(config_file_path))
            exit()
    else:
        config_file_path = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), 'mantis_config.cfg')

    # List of configuration parameters/arguments used by the three-tier
    # config priority loader.
    params = [
        Parameter(key='bedfile', required=True),
        Parameter(key='normal_filepath', arg_key='normal', required=True),
        Parameter(key='tumor_filepath', arg_key='tumor', required=True),
        Parameter(key='output_filepath', arg_key='output', required=True),
        Parameter(key='genome'),
        Parameter(key='mrq', default=25.0),
        Parameter(key='mlq', default=30.0),
        Parameter(key='mrl', default=35),
        Parameter(key='mlc', default=30),
        Parameter(key='mrr', default=3),
        Parameter(key='sd', default=3.0),
        Parameter(key='threads', default=1),
        Parameter(key='dif_threshold', default=0.4),
        Parameter(key='euc_threshold', default=0.187),
        Parameter(key='cos_threshold', default=0.07),
    ]

    # Use config/setting loader.
    config = load_settings(params, config_file_path, args)

    # Convert any variables that should be filepaths into the absolute
    # filepaths of the values.
    filepath_vars = [
        'bedfile', 
        'normal_filepath', 
        'tumor_filepath', 
        'output_filepath', 
        'genome',
        ]

    for key in filepath_vars:
        if key in config and config[key] is not None:
            config[key] = os.path.abspath(config[key])


        
    if 'bedfile' not in config:
        print('Error: BED file not provided!')
        exit(1)
    else:
        bedfile = os.path.abspath(config['bedfile'])
        if not os.path.isfile(bedfile):
            print('Error: {0} does not exist!'.format(bedfile))
            exit(1)
        else:
            config['bedfile'] = bedfile
            if not check_bedfile_format(config['bedfile']):
                exit(1)

    if 'normal_filepath' not in config:
        print('Error: Normal BAM file not provided!')
        exit()
    else:
        normal_filepath = os.path.abspath(config['normal_filepath'])
        if not os.path.isfile(normal_filepath):
            print('Error: {0} does not exist!'.format(normal_filepath))
            exit(1)
        else:
            # Make sure a corresponding .BAI index file exists
            if not os.path.isfile('{0}.bai'.format(normal_filepath)):
                print('Error: {0} needs corresponding .bai index file!'.format(normal_filepath))
                exit(1)
            config['normal_filepath'] = normal_filepath
    
    if 'tumor_filepath' not in config:
        print('Error: Tumor BAM file not provided!')
        exit(1)
    else:
        tumor_filepath = os.path.abspath(config['tumor_filepath'])
        if not os.path.isfile(tumor_filepath):
            print('Error: {0} does not exist!'.format(tumor_filepath))
            exit(1)
        else:
            # Make sure a corresponding .BAI index file exists
            if not os.path.isfile('{0}.bai'.format(tumor_filepath)):
                print('Error: {0} needs corresponding .bai index file!'.format(tumor_filepath))
                exit(1)            
            config['tumor_filepath'] = tumor_filepath

    if args.output is None:
        print('Error: Output filename not specified!')
        exit(1)
    else:
        config['output_filepath'] = os.path.abspath(args.output)
        # Create output directory if it doesn't exist.
        output_dir = os.path.dirname(config['output_filepath'])
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)         

    output_filepath = os.path.abspath(args.output)


    """
    Nothing below this line should need to be edited, as the code below
    handles the execution of the MANTIS program.
    """

    mantis_folder = os.path.dirname(os.path.realpath(__file__))
    kmer_repeat_counter = os.path.join(mantis_folder, 'kmer_repeat_counter.py')
    kmer_count_filter = os.path.join(mantis_folder, 'kmer_count_filter.py')
    instability_calculator = os.path.join(mantis_folder, 'calculate_instability.py')


    # First, use the k-mer repeat counter to count how many repeats of each
    # repeat unit (k-mer) each locus has, both in the normal and tumor file.
    # Results will be saved into a file for passing to the next phase.
    if '.' in output_filepath:
        kmer_count_output = output_filepath.split('.')
        kmer_count_output[-1] = '{0}.{1}'.format('kmer_counts', kmer_count_output[-1])
        kmer_count_output = '.'.join(kmer_count_output)
    else:
        kmer_count_output = output_filepath + '.kmer_counts'

    # First job - generate file with repeat counts for each n-multiple of kmers
    kmer_count_filepath = output_filepath + '.kmer_counts'
    command = [
        'python {0} '.format(kmer_repeat_counter),
        '-b {0} '.format(os.path.abspath(config['bedfile'])),
        '-n {0} '.format(os.path.abspath(config['normal_filepath'])),
        '-t {0} '.format(os.path.abspath(config['tumor_filepath'])),
        '-o {0} '.format(os.path.abspath(kmer_count_output)),
        '--min-read-quality {0} '.format(config['mrq']),
        '--min-locus-quality {0} '.format(config['mlq']),
        '--min-read-length {0} '.format(config['mrl']),
        '--genome {0} '.format(config['genome']),
        '--threads {0} '.format(config['threads']),
    ]

    print('\\\n'.join(command))

    print('Getting repeat counts for repeat units (k-mers) ...')
    sp = subprocess.Popen([' '.join(command)], stdout=subprocess.PIPE, shell=True)
    response = sp.communicate()[0]
    if sp.returncode != 0:
        print('Error with k-mer repeat count calculations; terminating program.')
        exit(1)
    print(response)        
    print('done.')
    

    # Run the outlier filter to get rid of outliers that are likely caused
    # due to sequencing artifacts or one-off biological events.

    filtered_kmer_counts = kmer_count_output.replace('kmer_counts', 'kmer_counts_filtered')

    command = [
        'python {0} '.format(kmer_count_filter),
        '-i {0} '.format(os.path.abspath(kmer_count_output)),
        '-o {0} '.format(os.path.abspath(filtered_kmer_counts)),
        '-mlc {0} '.format(config['mlc']),
        '-mrr {0} '.format(config['mrr']),
        '-sd {0}'.format(config['sd']),
    ]

    print('\\\n'.join(command))

    print('Filtering out outliers ... ')
    sp = subprocess.Popen([' '.join(command)], stdout=subprocess.PIPE, shell=True)
    response = sp.communicate()[0]
    if sp.returncode != 0:
        print('Error with k-mer repeat count filter; terminating program.')
        exit(1)
    print(response)        
    print('done.')


    # Calculate instability scores.
    # Analyze locus instability
    command = [
        'python {0} '.format(instability_calculator),
        '-i {0} '.format(filtered_kmer_counts),
        '-o {0} '.format(output_filepath),
        '--difference-threshold {0}'.format(config['dif_threshold']),
        '--distance-threshold {0}'.format(config['euc_threshold']),
        '--dissimilarity-threshold {0}'.format(config['cos_threshold']),
        ]

    print('\\\n'.join(command))

    print('Calculating instability scores ... ')
    sp = subprocess.Popen([' '.join(command)], stdout=subprocess.PIPE, shell=True)
    response = sp.communicate()[0]
    if sp.returncode != 0:
        print('Error with instability score calculations; terminating program.')
        exit(1)
    print(response)
    print('done.')

    print('\nMANTIS complete.')

