#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
@author: Anastassiya Zidkova and Martin Zidek
"""

import csv
import numpy as np
import pandas as pd
import numpy.random as nprnd
from scipy.signal import argrelextrema


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
  

def compute_fdr(input_file, gtf_file, range_number, rand_iter):
   limits = get_intervals(gtf_file)
   reads = get_reads(input_file, range_number)    
   real_freq = get_frequencies(reads, range_number)
   
   n = len(real_freq[:,0])
   real_freq_freq = get_frequencies(real_freq[:,1], 0)
   a = []
   for i in range (0, real_freq_freq.shape[0]):
        row = real_freq_freq[i].tolist()
        ph = np.sum(real_freq_freq[i:,1]) / n
        row.append(ph)
        a.append(row)
        
   real_freq_ph = np.array(a)
    
   rand = randomize(limits, reads)
    
   rand_freq = get_frequencies(rand, range_number)
   maxima_indices = argrelextrema(rand_freq[:,1], np.greater)
   maxima = (rand_freq[maxima_indices[0].tolist(),:])[:,1]
   maxima_freq_first = get_frequencies(maxima, 0)
   rand_freq_list = [] 
   
   n = []
   for i in range(0, rand_iter):
        rand = randomize(limits, reads, i+1)
        n.append(len(np.unique(rand)))
        rand_freq = get_frequencies(rand, range_number)
        
        rand_freq_list.append(rand_freq)
        
        maxima_indices = argrelextrema(rand_freq[:,1], np.greater)
        maxima = (rand_freq[maxima_indices[0].tolist(),:])[:,1]
        maxima_freq = get_frequencies(maxima, 0)
        maxima_freq_first = add_freq_tables(maxima_freq_first, maxima_freq)
    
   p_table = []
   for r in range(0, maxima_freq_first.shape[0]):    
        row = maxima_freq_first[r]
        p_row = [row[0]]
        for c in range(1, len(row)):
            ph = sum(maxima_freq_first[r:,c]) / n[i]
            p_row.append(ph)
        p_table.append(p_row)
    
   md_table = []
   for i in range(0, len(p_table)):
        row = [p_table[i][0]]
        mean = np.mean(p_table[i][1:])
        dev = np.std(p_table[i][1:])
        row.append(mean)
        row.append(dev)
        md_table.append(row)
    
   fdr_table = []
   for i in range(0, len(p_table)):
        row = [p_table[i][0]]
    
        index = np.where(real_freq_ph[:,0] == row[0])[0]
        if (len(index) > 0):
            fdr = (p_table[i][1] + p_table[i][2]) / real_freq_ph[index[0]][2]
            row.append(fdr)     
            fdr_table.append(row)
   output_file = input_file[:-4] + '_fdr.csv' 
   dump_to_file(fdr_table, output_file)
if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description="""
        ------------> iCLIPit <------------
        Script to calculate FDR (False Discovery Rate) threshold.
        Output format is in .csv format and has two columns: <reads number><FDR probability>
        Output has "_fdr.csv" extension.
        Created by Anastassiya Zidkova and Martin Zidek, 2014
        """, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-f", "--file", dest="input_file",required=True,
                        help="input file created by extracting first 4 columns from sam file without header. \n Requred", metavar="file")
    parser.add_argument("-g", "--gtf", dest="gtf_file", required=True,
                        help="input file in gtf format. \n Requred", metavar="file")    
    parser.add_argument("-r", "--range_number", dest="range_number",required=True,
                        help="range of nucleotides in one direction from the position to calculate number of reads starting in the range: position +- range_number", 
                        default=15, type=int,metavar="integer")
    parser.add_argument("-i", "--rand_iter", dest="rand_iter",required=True,
                        help="number of iterations to generate randomized datasets from input file", default=100, type=int,metavar="integer")
                        
    args = parser.parse_args()
    compute_fdr(args.input_file, args.gtf_file, args.range_number, args.rand_iter)
