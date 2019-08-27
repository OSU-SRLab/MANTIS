Output Files
============

For each run, MANTIS generates four output files. The file ending in .status contains the final MANTIS result and MSI call. The other three files contain microsatellite locus-specific information used in the process. File formats are described below.

*_kmer_counts.txt
-----------------

This file contains, for each locus with sufficient coverage in normal and tumor, the number of reads in the tumor and normal sample which support different numbers of k-mer repeats. Columns are as follows:

|Column Number|Column Name|Description|
|---|---|---|
|1|Locus|Position of the microsatellite locus within the genome|
|2|Repeats|Number of k-mer repeats supported by the reads summarized in this line|
|3|Normal|Number of reads in the normal sample which support this number of repeats at this locus|
|4|Tumor|Number of reads in the tumor sample which support this number of repeats at this locus|

Note that minimum quality and coverage filters have been enforced, but outliners have not been filtered in this file.

*_kmer_counts_filtered.txt
--------------------------

Similarly to kmer\_counts.txt, this file summarizes the numbers of reads in the tumor and normal sample which support different numbers of k-mer repeats. However, the contents of this file reflect outlier filtering, and therefore consist of all reads that were used to compute the MANTIS score.  
File format is identical to that of kmer\_counts.txt.

*.txt
-----

This file contains the locus-specific scores computed by MANTIS with the difference, distance, and dissimilarity metrics used to compare them. The last line contains the average score for each metric. The average difference is used by MANTIS to call a sample MSS or MSI. Columns are as follows:

|Column Number|Column Name|Description|
|---|---|---|
|1|Locus|Position of the microsatellite locus within the genome|
|2|Normal_Reads|Number of reads in the normal sample passing all filters used to compute this locus-specific score|
|3|Tumor_Reads|Number of reads in the tumor sample passing all filters used to compute this locus-specific score|
|4|Difference|Locus-specific score computed by the stepwise difference metric (see the [MANTIS publication](https://www.ncbi.nlm.nih.gov/pubmed/?term=27980218) for algorithmic details)|
|5|Distance|Locus-specific score computed by the Euclidean distance metric|
|6|Dissimilarity|Locus-specific score computed by the cosine dissimilarity metric|

*.status
--------

This file contains the final MANTIS scores for the tumor-normal pair with each of the three distance metrics computed, along with the recommended thresholds for MSI calls. The MSI status per these recommendations is also included, along with a notice that we recommend usage of the stepwise difference call.