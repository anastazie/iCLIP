#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: Anastassiya Zidkova and Martin Zidek
"""


import time
import csv
import os
import re
import numpy as np
import pandas as pd
from pylab import *

#get start and stop positions for each gene and gene name as list of lists
def get_gene_interval(input_file):
    start=[]
    stop=[]
    name=[]
    with open(input_file) as fd:
        fd.readline()
        for line in fd:
            start.append(line.split()[3])
            stop.append(line.split()[4])
            name.append(line.replace(';', '\t').replace('"', '\t').split('\t')[9])
        start=map(float, start)
        #print start        
        stop=map(float, stop)
        return name, start, stop

#compute number of total mapped reads from *_freq_real.csv
def get_total_reads(input_file):
    with open(input_file) as fd:
        reader = csv.reader(fd,delimiter="\t")
        n=0        
        for line in reader:
            n = n + float(line[1])
        return n
        
#get fraction of reads mapped to each gene(interval) in percent
def get_read_counts(input_file, intervals):
    #list command is important to iterate over .csv file several times
        reader = list(csv.reader(open(input_file), delimiter="\t"))
        n=get_total_reads(input_file)
        e=[] 
              
        for i in range(len(intervals[0])):
            count=0            
            for line in reader:
                if float(line[0]) in range(int(intervals[1][i]), int(intervals[2][i]+1)):
                    #print intervals[1][i], intervals[2][i]
                    count=count + float(line[1])
            #print count
            e.append((count/n)*100)#count percentage of mapped reads per each gene
        #print e 
        return e

# calculate number of genes
def calc_gene_number(gtf_file):
    count=0
    with open(gtf_file) as f:
        for line in f:
            if "##" not in line:
                count+=1
        return count
      
#create barplots with y - Mapped reads percentage and x - gene name and save as .png file
def create_barplot(input_file, gtf_file, title, l, number_gene=None):
    figure()
    if number_gene == None:
        n=calc_gene_number(gtf_file)
    elif number_gene != None:
        n=number_gene
    ind = np.arange(n)
    bar(ind, input_file[:n], color="green")
    ylabel('Reads fraction, %')
    xticks(ind, l[0][:n], rotation = 40)
    suptitle(title)
    yticks(range(0,110,10))
    png_file = title + '_barplot.png' 
    savefig(png_file, bbox_inches='tight')
    
if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description="""
        ------------> iCLIPit <------------
        Script to create barplot for iCLIP analysis.
        Created by Anastassiya Zidkova and Martin Zidek, 2014
        """, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-f", "--file", dest="input_file",required=True,
                        help="input file in csv format with two columns, no header. \n Format: <crosslink position> <cDNA count> \n Requred", metavar="file")
    parser.add_argument("-g", "--gtf", dest="gtf_file",required=True,
                        help="input file in gtf format. \n Requred", metavar="file")
    parser.add_argument("-t", "--title", dest="title", default="Barplot",
                        help="barplot title", metavar="name")
    parser.add_argument("-n", "--nIntervals", default=None, type=int, dest="n",
                        help="number of intervals in genome", metavar="integer")
                        
    args = parser.parse_args()
    l=get_gene_interval(args.gtf_file)
    a=get_read_counts(args.input_file, l)
    create_barplot(a, args.gtf_file, args.title, l, args.n)

