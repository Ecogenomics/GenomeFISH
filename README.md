# GenomeFISH

Scripts are contained in the scripts directory and are written in Python 2.x.



## Simulating <i>in silico</i> probes

Script: in_silico_probes.py

Requirements:
* [biolib](https://github.com/dparks1134/biolib) Python library
* [BLASTn](http://blast.ncbi.nlm.nih.gov) >= 2.6.0+ to be on your system path
* [melt.pl](http://unafold.rna.albany.edu/?q=unafold-man-pages/melt.pl) from the [UNAFold](http://unafold.rna.albany.edu/?q=unafold-man-pages) software library to be on your system path

Calculates the percent identity and free energy error between in silico probes from a reference genome and a target genome. It can be run using:
> in_silico_probes.py <genome_dir> <ani_matrix> <output_dir>

where <genome_dir> contain 2 genomic FASTA files for two or more genomes, <ani_matrix> is a file with pairwise average nucleotide identity between genomes, and <output_dir> is the desired output directory for results. In order to reduce computational requirements <i>in silico</i> probes are only simulated every 120 bp (i.e. same length as the probe). At this step size it takes ~2 hours to compare a pair of genomes. You can do multiple genome comparisons in parallel though so this should help a bit. The length of the probes can be changed with the "--probe_size" parameter and the spacing between probes changed with the "--probe_step_size" parameters. The number of CPUs to use can be specified with the "--threads" parameter. Using multiple CPUs will substantially reduce processing time when multiple pairs of genomes are being processed. For other optional parameters see the command line help (i.e. in_silico_probes.py -h).
