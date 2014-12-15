#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
@author: Anastassiya Zidkova and Martin Zidek
"""

import csv
import numpy as np
import pandas as pd

def get_frequencies(values, margin):
    unique_values = np.unique(np.array(values))    
    freq = np.c_[unique_values, np.zeros(unique_values.shape[0])]
    
    for val in values:
        indices = np.where(((freq[:,0] - margin) <= val) & ((freq[:,0] + margin) >= val))[0]
        if (len(indices) > 0):
            for index in indices:
                freq[index][1] += 1            
    return freq
    

def get_reads(input_file, dump_to_file=False):
    data = pd.io.parsers.read_csv(input_file, sep='\t', skiprows=0, header=None)
    reads = data.ix[:,3]  
    freq = get_frequencies(reads, 0)    
    freq = freq[freq[:,0].argsort()]
    return reads
        
def dump_to_file(data, filename):
    with open(filename, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for line in data:
            writer.writerow(line)   


def get_real_freq(input_file):
    reads = get_reads(input_file, dump_to_file=False)
    real_freq = get_frequencies(reads, 0)
    output_file = input_file[:-4] + '_real_freq.csv'
    dump_to_file(real_freq, output_file)
    
if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description="""
        ------------> iCLIPit <------------
        Script to calculate number of reads starting at each position.
        Output format is in .csv format and has two columns: <position><number of reads starting in this position>
        Output has "_real_freq.csv" extension.
        Created by Anastassiya Zidkova and Martin Zidek, 2014
        """, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-f", "--file", dest="input_file",required=True,
                        help="input file created by extracting first 4 columns from sam file without header. \n Requred", metavar="file")
                        
    args = parser.parse_args()
    get_real_freq(args.input_file)