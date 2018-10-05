#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__prog_name__ = 'in_silico_probes.py'
__prog_desc__ = 'Identify in silico probes which hybridize between genomes.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2018'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import ntpath
import argparse
import logging
import string
import itertools
import shutil
import tempfile
import uuid
import multiprocessing as mp
from collections import defaultdict, namedtuple

from biolib.external.execute import check_dependencies
from biolib.seq_tk import rev_comp
from biolib.seq_io import read_fasta
from biolib.logger import logger_setup


class ProbeMatches(object):
    """Identify near identical matches of probes between genomes."""
    
    def __init__(self, 
                    temp,
                    na_plus,
                    ct,
                    free_energy_threshold,
                    output_dir):
        """Initialization."""
        
        check_dependencies(['melt.pl', 'blastn', 'makeblastdb'])
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        self.output_dir = output_dir
        
        logger_setup(output_dir, "in_silico_probes.log", "in_silico_probes", __version__, False)
        self.logger = logging.getLogger('timestamp')

        self.temp = temp
        self.na_plus = na_plus
        self.ct = ct
        self.free_energy_threshold = free_energy_threshold
        
        self.output_fmt = '6 qseqid qlen qseq sseqid slen sseq length mismatch gaps pident bitscore evalue'

        self.BlastHit = namedtuple('BlastHit', """query_id
                                                query_len
                                                query_aln_seq
                                                subject_id
                                                subject_len
                                                subject_aln_seq
                                                aln_len
                                                mismatch
                                                gaps
                                                perc_identity
                                                bitscore
                                                evalue""")
        
    def _blastn(self, query_seqs, nucl_db, output_file, evalue=1e-3, max_matches=500, task='megablast'):
        """Apply blastn to query file.

        Finds homologs to query sequences using blastn homology search
        against a nucleotide database. Hit can be reported using  either
        the 'standard' table 6 format or the following 'custom' format:
            qseqid qlen sseqid slen length pident evalue bitscore


        Parameters
        ----------
        query_seqs : str
            File containing query sequences.
        nucl_db : str
            File containing blastn formatted database.
        output_file : str
            Output file containing blastn results.
        evalue : float
            E-value threshold used to identify homologs.
        max_matches : int
            Maximum hits per query sequence.
        output_fmt : str
            Specified output format of blast table: standard or custom.
        """

        cmd = "blastn -num_threads %d" % self.threads_per_job
        cmd += " -query %s -db %s -out %s -evalue %g" % (query_seqs, nucl_db, output_file, evalue)
        cmd += " -max_target_seqs %d" % max_matches
        cmd += " -task %s" % task
        cmd += " -outfmt '%s'" % self.output_fmt
        cmd += " -gapopen 1000 -gapextend 1000 -penalty -2 -reward 3"
        cmd += " -dust no -soft_masking false"
        os.system(cmd)
        
    def _read_hit(self, table):
        """Generator function to read hits from a blast output table.

        Parameters
        ----------
        table : str
            Name of table to read.
        table_fmt : str
            Specified output format of blast table: standard or custom.

        Yields
        ------
        namedtuple
            Information about blast hit.
        """

        if table.endswith('.gz'):
            open_file = gzip.open
        else:
            open_file = open

        for line in open_file(table):
            line_split = line.split('\t')
            hit = self.BlastHit(query_id=line_split[0],
                            query_len=int(line_split[1]),
                            query_aln_seq=line_split[2],
                            subject_id=line_split[3],
                            subject_len=int(line_split[4]),
                            subject_aln_seq=line_split[5],
                            aln_len=int(line_split[6]),
                            mismatch=int(line_split[7]),
                            gaps=int(line_split[8]),
                            perc_identity=float(line_split[9]),
                            bitscore=float(line_split[10]),
                            evalue=float(line_split[11]))

            yield hit
        
    def _free_energy(self, seq1, seq2):
        """Calculate free energy between two sequences."""
        
        seq1_file = 'seq1_' + str(uuid.uuid4())
        fout = open(seq1_file, 'w')
        fout.write(seq1 + '\n')
        fout.close()
        
        seq2_file = 'seq2_' + str(uuid.uuid4())
        fout = open(seq2_file, 'w')
        fout.write(seq2 + '\n')
        fout.close()
        
        cmd = 'melt.pl -n DNA -t %g -N %g -C %g -p %s %s > /dev/null' % (self.temp, 
                                                                        self.na_plus, 
                                                                        self.ct,
                                                                        seq1_file, 
                                                                        seq2_file)
        os.system(cmd)
        
        lines = open('%s-%s.dG' % (seq1_file, seq2_file)).readlines()
        results = lines[-1].split()
        dG = float(results[1])

        os.remove(seq1_file)
        os.remove(seq2_file)
        cmd = 'rm %s-%s.*' % (seq1_file, seq2_file)
        os.system(cmd)
        
        return dG
        
    def _free_energy_of_formation(self, probe, target_seq): 
        """Determine free energy of formation between probe and target sequence."""
        
        assert(len(probe) == len(target_seq))
        
        # dG_perfect_complement; free energy for perfect complement
        pc = self._free_energy(probe, rev_comp(probe))
        
        # dG_mismatched; free energy for target sequence
        mm = self._free_energy(probe, target_seq)
        
        # ddG = -dG_perfect_complement + dG_mismatched
        ddG = -pc + mm
        
        return ddG

    def __workerThread(self, 
                        probe_size,
                        probe_step_size,
                        mismatch,
                        min_aln_len,
                        keep_fragments,
                        results_dir,
                        queueIn, 
                        queueOut):
        """Process each data item in parallel.
        
        The reference genome is the genome from which probes are being
        designed. The aim is to determine how many of these 
        reference probes will hybridize to the target genome. 

        To determine the number of reference probes which will hybridize
        to the target genome, the target genome is fragmented
        into probe sized windows to determine how many of these are 
        nearly identical to the reference genome.
        """

        while True:
            ref_genome, target_genome = queueIn.get(block=True, timeout=None)
            if ref_genome == None:
                break

            ref_name = ntpath.basename(ref_genome).replace('.fasta', '').replace('.fna', '')
            target_name = ntpath.basename(target_genome).replace('.fasta', '').replace('.fna', '')
            
            if keep_fragments:
                fragment_dir = os.path.join(results_dir, 'fragments')
            else:
                fragment_dir = tempfile.mkdtemp()

            # count total number of reference genome probes
            ref_seqs = read_fasta(ref_genome)
            ref_genome_size = 0
            num_ref_probes = 0
            for seq in ref_seqs.values():
                num_ref_probes += (len(seq)-probe_size)/probe_step_size + 1 #sum([1 for i in range(0, len(seq)-probe_size, probe_step_size)])
                ref_genome_size += len(seq)

            # fragment target genome into probe sized windows
            window_file = os.path.join(fragment_dir, ref_name + '~' + target_name + '.fna')
            fout = open(window_file, 'w')
            target_seqs = read_fasta(target_genome)
            num_target_probes = 0
            target_windows = {}
            target_genome_size = 0
            for seq in target_seqs.values():
                target_genome_size += len(seq)
                
                for i in range(0, len(seq)-probe_size, probe_step_size):
                    fout.write('>probe_%d\n' % num_target_probes)
                    fout.write(seq[i:i+probe_size] + '\n')
                    target_windows[str(num_target_probes)] = seq[i:i+probe_size]
                    num_target_probes += 1
            fout.close()
                  
            # BLAST target probes against reference genome
            output_table = os.path.join(results_dir, ref_name + '~' + target_name + '.blast_hits.tsv')
            self._blastn(window_file, 
                            ref_genome, 
                            output_table, 
                            evalue=1e-2, 
                            max_matches=1, 
                            task='dc-megablast')

            window_hits = set()
            failed_similarity_test = set()
            failed_free_energy_test = set()
            output_file = os.path.join(results_dir, ref_name + '~' + target_name + '.probe_hits.tsv')
            fout = open(output_file, 'w')
            fout.write('Probe ID\tSubject ID\tProbe percent alignment\tPercent identity\tAdjusted percent identity\tFree energy of formation\n')
            for hit in self._read_hit(output_table):
                adj_aln_len = hit.aln_len - hit.gaps
                query_aln_frac = adj_aln_len * 100.0 / hit.query_len
                adjusted_perc_identity = (adj_aln_len - hit.mismatch) * 100.0 / hit.query_len

                if (query_aln_frac >= (100*min_aln_len)
                    and adjusted_perc_identity >= (100*(1.0 - mismatch))):

                    if hit.query_id not in window_hits:
                        probe = hit.subject_aln_seq
                        target_seq = rev_comp(hit.query_aln_seq) # don't want genomic region on same strand, but likely hybridization on the other strand
                        ddG = self._free_energy_of_formation(probe, target_seq)
                        if ddG <= self.free_energy_threshold:
                            window_hits.add(hit.query_id)
                            fout.write('%s\t%s\t%.1f\t%.1f\t%.1f\t%.2f\n' % (hit.query_id, 
                                                                                hit.subject_id, 
                                                                                query_aln_frac,
                                                                                hit.perc_identity,
                                                                                adjusted_perc_identity,
                                                                                ddG))
                        else:
                            failed_free_energy_test.add(hit.query_id)
                else:
                    failed_similarity_test.add(hit.query_id)
            fout.close()
            
            num_failed_free_energy_test = len(failed_free_energy_test - window_hits)
            num_failed_similarity_test = len(target_windows) - num_failed_free_energy_test - len(window_hits) # need to account for probes with no blast hit
            
            output_file = os.path.join(results_dir, ref_name + '~' + target_name + '.summary.tsv')
            fout = open(output_file, 'w')
            fout.write('Reference ID\tReference genome size (bp)')
            fout.write('\tTarget ID\tTarget genome size (bp)')
            fout.write('\tNo. reference probes\tNo. target probes\tNo. hybridized probes\tNo. probes failing genomic similarity test\tNo. probes failing free energy test\n')
            
            fout.write('%s\t%d' % (ref_name, ref_genome_size))
            fout.write('\t%s\t%d' % (target_name, target_genome_size))
            fout.write('\t%d\t%d\t%d\t%d\t%d' % (num_ref_probes, len(target_windows), len(window_hits), num_failed_similarity_test, num_failed_free_energy_test))
            fout.write('\n')
            fout.close()
            
            if not keep_fragments:
                shutil.rmtree(fragment_dir)

            # allow results to be processed or written to file
            queueOut.put(ref_name)

    def __writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a == None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) genome pairs.' % (processedItems, numDataItems, float(processedItems)*100/numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')
        
    def read_ani_matrix(self, ani_matrix):
        """Read ANI matrix."""
        
        ani = defaultdict(lambda: defaultdict(float))
        of = defaultdict(lambda: defaultdict(float))
        with open(ani_matrix) as f:
            header = f.readline()
            gids = header.strip().split('\t')

            row_index = 0
            for line in f:
                line_split = line.strip().split('\t')
                
                row_gid = line_split[0]
                for c, v in enumerate(line_split[1:]):
                    if row_index < c:
                        of[row_gid][gids[c]] = float(v)
                        of[gids[c]][row_gid] = float(v)
                    elif row_index > c:
                        ani[row_gid][gids[c]] = float(v)
                        ani[gids[c]][row_gid] = float(v)
                    
                row_index += 1
                    
        return ani, of
        
    def build_table(self, results_dir, ani_matrix):
        """Build table with results."""
        
        ani, of = self.read_ani_matrix(ani_matrix)
        
        output_table = os.path.join(self.output_dir, 'probe_table.tsv')

        fout = open(output_table, 'w')
        fout.write('Reference ID\tReference genome size (bp)')
        fout.write('\tTarget ID\tTarget genome size')
        fout.write('\tANI (%)\tShared genes (%)')
        fout.write('\tNo. reference probes\tNo. target probes\tNo. hybridized probes\tNo. probes failing genomic similarity test\tNo. probes failing free energy test\n')
        for result_file in os.listdir(results_dir):
            if not result_file.endswith('.summary.tsv'):
                continue
               
            r = os.path.join(results_dir, result_file)
            if r == output_table:
                continue
                
            with open(r) as f:
                f.readline()
                
                for line in f:
                    line_split = line.strip().split('\t')
                    ref_name = line_split[0]
                    ref_genome_size = line_split[1]
                    target_name = line_split[2]
                    target_genome_size = line_split[3]
                    num_ref_probes = int(line_split[4])
                    num_target_probes = int(line_split[5])
                    hybridized_probes = int(line_split[6])
                    failed_similarity_test = int(line_split[7])
                    failed_free_energy_test = int(line_split[8])
   
                    fout.write('%s\t%s' % (ref_name, ref_genome_size))
                    fout.write('\t%s\t%s' % (target_name, target_genome_size))
                    fout.write('\t%.2f\t%.2f' % (ani[ref_name][target_name], of[ref_name][target_name]))
                    fout.write('\t%d\t%d\t%d\t%d\t%d' % (num_ref_probes, num_target_probes, hybridized_probes, failed_similarity_test, failed_free_energy_test))
                    fout.write('\n')
                 
        fout.close()

    def run(self, 
                genome_dir, 
                probe_size,
                probe_step_size,
                mismatch,
                min_aln_len,
                ani_matrix, 
                keep_fragments,
                threads):
                
        # create results and fragments directory
        results_dir = os.path.join(self.output_dir, 'results')
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
                
        if keep_fragments:
            fragment_dir = os.path.join(results_dir, 'fragments')
            if not os.path.exists(fragment_dir):
                os.makedirs(fragment_dir)

        # read gene files to process
        genome_files = []
        for genome_file in os.listdir(genome_dir):
            if genome_file == 'GOUP_E_MS2403.fasta':
                continue
                
            if genome_file.endswith('.fasta') or genome_file.endswith('.fna'): 
                genome_files.append(os.path.join(genome_dir, genome_file))
                
        # create BLAST databases
        self.logger.info('Creating BLAST databases.')
        for gf in genome_files:
            os.system('makeblastdb -dbtype nucl -in %s > /dev/null' % gf)
            
        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        num_pairs = 0
        for gf1, gf2 in itertools.permutations(genome_files, 2):
            workerQueue.put((gf1,gf2))
            num_pairs += 1
            
        self.threads_per_job = max(1, threads/num_pairs)
        self.logger.info('Processing %d pairs with %d CPUs and %d CPUs per job.' % (num_pairs, threads, self.threads_per_job))

        for _ in range(threads):
            workerQueue.put((None, None))

        try:
            workerProc = [mp.Process(target = self.__workerThread, args = (probe_size,
                                                                            probe_step_size,
                                                                            mismatch,
                                                                            min_aln_len,
                                                                            keep_fragments,
                                                                            results_dir,
                                                                            workerQueue, 
                                                                            writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target = self.__writerThread, args = (num_pairs, writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()
        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate()
            
        self.build_table(results_dir, ani_matrix)
        
        self.logger.info('Done.')

if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genome_dir', help='directory with genomes to process')
    parser.add_argument('ani_matrix', help='matrix with ANI results')
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('--probe_size', help='probe/window size', type=int, default=120)
    parser.add_argument('--probe_step_size', help='step size for generating in silico probes', type=int, default=120)
    parser.add_argument('--mismatch', help='maximum percent mismatch between probe and target for hybridization [0,1]', type=float, default=0.15)
    parser.add_argument('--min_aln_len', help='minimum percent alignment length of probe for hybridization [0,1]', type=float, default=0.85)
    parser.add_argument('--temp', help='temperature for calculating the free energy penalty', type=float, default=48)
    parser.add_argument('--na_plus', help='concentration of Na+ for calculating the free energy penalty', type=float, default=0.825)
    parser.add_argument('--ct', help='nucleic acid concentration in M for calculating the free energy penalty', type=float, default=0.000000005)
    parser.add_argument('--free_energy_threshold', help='threshold for free energy penalty', type=float, default=2)
    parser.add_argument('--keep_fragments', action='store_true', help='retain file with query fragments')
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=1)

    args = parser.parse_args()

    try:
        p = ProbeMatches(args.temp,
                            args.na_plus,
                            args.ct,
                            args.free_energy_threshold,
                            args.output_dir)
        p.run(args.genome_dir, 
                args.probe_size,
                args.probe_step_size,
                args.mismatch,
                args.min_aln_len,
                args.ani_matrix, 
                args.keep_fragments,
                args.threads)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
