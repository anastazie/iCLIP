#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
@author: Anastassiya Zidkova and Martin Zidek
"""
import csv
import numpy as np
import pandas as pd
import numpy.random as nprnd
from Kmers import Kmers
from RefSeq import RefSeq


def get_frequencies(values, margin):
    unique_values = np.unique(np.array(values))
    
    freq = np.c_[unique_values, np.zeros(unique_values.shape[0])]
    
    for val in values:
        indices = np.where(((freq[:,0] - margin) <= val) & ((freq[:,0] + margin) >= val))[0]
        if (len(indices) > 0):
            for index in indices:
                freq[index][1] += 1
            
    return freq
    

def get_reads(input_file, n, dump_to_file=False):
    data = pd.io.parsers.read_csv(input_file, sep='\t', skiprows=0, header=None)
    
    #get column with read position as numpy series
    reads = data.ix[:,3]
  
    freq = get_frequencies(reads, n)
    
    freq = freq[freq[:,0].argsort()]

    return reads
        

def get_intervals(input_file):
    df = pd.io.parsers.read_csv(input_file, sep='\t', skiprows=0, header=None)
    
    limits = df.ix[:,3:4]
    
    return limits
    
def randomize(limits, reads, seed=0):
    freq = np.c_[limits, np.zeros(limits.shape[0])]
    
    for read in reads:
        a = np.where((limits[3] < read) & (limits[4] > read))
        
        if (a[0].size > 0):
            freq[a[0][0]][2] += 1
        #else:
            #print read
    rand = []
        
    for row in freq:
        nprnd.seed(seed)
        rand += nprnd.randint(row[0], row[1], size=row[2]).tolist()
    
    return rand

def add_freq_tables(orig, new):
    orig_index = orig.shape[1]
    orig = np.c_[orig, np.zeros(orig.shape[0])]
    
    for row in orig:
        index = (np.where(new[:,0] == row[0]))[0]
        #print index
        if (len(index) > 0):
            row[orig_index] = new[index[0]][1]
        else:
            row[orig_index] = 0
            
    return orig
 
def dump_to_file(data, filename):
    with open(filename, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for line in data:
            writer.writerow(line)   


def get_positions_above_threshold(threshold, freq_table):
    pos_array = np.where(freq_table[:,1] >= threshold)[0]
    
    if len(pos_array) > 0:
        return freq_table[pos_array]
    else:
        print 'Error: no positions above threshold!'

def get_reference_subsequences(positions_table, reference, l):
    sequences = []    
    for row in positions_table:
        s = reference.get_position(int(row[0]), l)
        sequences.append(s)
    return sequences



def compute_z_score(input_file, ref_seq, gtf_file, threshold, range_number, rand_iter, k, l):
    rand_freq_list = []
    reads = get_reads(input_file,range_number)       
    p = Kmers(k)
    
    r = RefSeq(ref_seq)
    
    real_freq = get_frequencies(reads, range_number)
    positions = get_positions_above_threshold(threshold, real_freq)
    
    references = get_reference_subsequences(positions, r, l)
    output_file = input_file[:-4] + '_' + str(k) + '_seq.txt'
    dump_to_file(references, output_file)
    real_pentamere_freq = p.calculate_frequencies(references)
    limits = get_intervals(gtf_file)
    rand = randomize(limits, reads)
    n = []
    for i in range(0, rand_iter):
        rand = randomize(limits, reads, i+1)
        n.append(len(np.unique(rand)))
        rand_freq = get_frequencies(rand, range_number)
        
        rand_freq_list.append(rand_freq)

    #rand_freq = get_frequencies(rand, 15)
    
    rand_pent_table = np.array(np.empty((p.high+1, len(rand_freq_list))))
    #rand_pent_table[:,0] = p.__list__
    
    i = 0
    for rand in rand_freq_list:
        r_positions = get_positions_above_threshold(threshold, rand)
        references = get_reference_subsequences(r_positions, r, l)
        rand_p_freq = p.calculate_frequencies(references)
        ar = np.array(rand_p_freq)
        rand_pent_table[:,i] = ar[:,1]
        i += 1

    ps = []        
    for i in range(0, p.count):
        mean = np.mean(rand_pent_table[i])
        dev = np.std(rand_pent_table[i])
        row = [p.__list__[i], mean, dev]
        ps.append(row)
    
    final_p = []        
    for i in range(0, p.count):
        val = (real_pentamere_freq[i][1] - ps[i][1])/ps[i][2]
        row = [ps[i][0], real_pentamere_freq[i][1], val]
        final_p.append(row)
    output_file = input_file[:-4] + '_' + str(k) + '_zscores.csv'    
    dump_to_file(final_p, output_file)

if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description="""
        ------------> iCLIPit <------------
        Script z-scores for k-mers.
        Two output files are produced:
        <kmer>_seq.txt contains sequence for positions above threshold +- *l*
        <kmer>_zscores.csv contains three columns: <kmer><occurence in real data above threshold><zscore>
        Created by Anastassiya Zidkova and Martin Zidek, 2014
        """, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-f", "--file", dest="input_file",required=True,
                        help="input file created by extracting first 4 columns from sam file without header. \n Requred", metavar="file")
    parser.add_argument("-fa", "--ref", dest="ref_file",required=True,
                        help="Fasta file containing reference nucleotide sequence in fa.tab format with two columns: <organism name><reference sequence> \n Requred", metavar="file")                    
    parser.add_argument("-g", "--gtf", dest="gtf_file", required=True,
                        help="input file in gtf format. \n Requred", metavar="file")  
    parser.add_argument("-t", "--threshold", dest="t",required=True,
                        help="FDR threshold value \n Requred", 
                        default=None, type=int, metavar="integer")                    
    parser.add_argument("-r", "--range_number", dest="range_number",required=True,
                        help="range of nucleotides in one direction from the position to calculate number of reads starting in the range: position +- range_number", 
                        default=15, type=int, metavar="integer")
    parser.add_argument("-i", "--rand_iter", dest="rand_iter",required=True,
                        help="number of iterations to generate randomized datasets from input file", default=100, type=int, metavar="integer")
    parser.add_argument("-k", "--kmer", dest="kmer",required=True,
                        help="k-mer type: pentamer, hexamer, etc.", default=5, type=int, metavar="integer") 
    parser.add_argument("-l", "--seq_length", dest="seq_length",required=True,
                        help="Length of sequence to output into seq.txt file: position +- length", default=10, type=int, metavar="integer")                                           
    args = parser.parse_args()
    compute_z_score(args.input_file, args.ref_file, args.gtf_file, args.t, args.range_number, args.rand_iter, args.kmer, args.seq_length)