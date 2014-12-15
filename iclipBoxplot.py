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


#extract start position of the gene - not needed to make script easier
#def get_start_pos(input_file):
#    start=[]
#    with open(input_file) as fd:
#        fd.readline()
#        for line in fd:
#            start.append(line.split()[3])
#        start=map(float, start)
#        return start

#get size of each gene by subtractiong start from stop
def get_gene_size(input_file):
    start=[]
    stop=[]
    name=[]
    with open(input_file) as fd:
        fd.readline()
        for line in fd:
            start.append(line.split()[3])#extract start use this part instead of extracting start positions using get_start_pos
            stop.append(line.split()[4])#extract stop
            name.append(line.replace(';', '\t').replace('"', '\t').split('\t')[9])#extract gene name
        start=map(float, start)
        #print start        
        stop=map(float, stop)
        col=[i - j for i, j in zip(stop, start)] #compute gene size   
        return name, col, start, stop#return list of lists

#divide gene ino n equal parts
def divide_gene(input_file, n):
    parts = [x/n for x in input_file]
    return parts

#compute gene intervals by adding k*parts to the start position
def determine_gene_intervals(input_file, table, n):
    pos=[]    
    for i in range (0, (len(l[0]))): 
        #print i
        pos0=[]
        for k in range(1,(n+1)):
            #print k
            interval=table[2][i]+k*parts[i]#compute each interval for each gene
            pos0.append(interval)
        #print pos0
        pos0=map(int, pos0)
        pos.append(pos0)
    #pos=zip(l, pos)
    #print pos 
    return l, pos

#count cDNA occurence in each computed interval
#create intervals distribution by multipliing interval number by cDNA occurence
#Example: 5 cDNA in 7th interval, e will be extended by following: 7,7,7,7,7
#creates one list per each gene 
def count_intervals(input_file, gtf_file, interval, p, number_gene):
    if number_gene == None:
        n=calc_gene_number(gtf_file)
    elif number_gene != None:
        n=number_gene
    with open(input_file) as fd:
        reader = csv.reader(fd,delimiter=",")
        c=[]        
        for line in reader:
            c.append(float(line[0]))
        d=[]
        for i in range(0,(len(l[0][:n]))):
            e=[]
            for k in range(0,p-1):
                e.extend([k]*len([x for x in c if interval[1][i][k] <= x < interval[1][i][k+1]]))
            #print e            
            #e=list("".join(e))
            e=map(int, e)
            d.append(e)        
        return d
        
# calculate number of genes (lines) in GTF file
def calc_gene_number(gtf_file):
    count=0
    with open(gtf_file) as f:
        for line in f:
            if "##" not in line:
                count+=1
        return count

#create boxplot, add x and y labels and save as png file
def create_boxplot(input_file, gtf_file, title, l, number_gene=None):
    figure()    
    if number_gene == None:
        n=calc_gene_number(gtf_file)
    elif number_gene != None:
        n=number_gene
    boxplot(input_file[:n],vert=False)
    yticks(range(1,n+1), l[0][:n])
    suptitle(title)
    xlabel('Distance from the start of the gene')
    png_file = title + '_boxplot.png' 
    savefig(png_file, bbox_inches='tight')      

if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description="""
        ------------> iCLIPit <------------
        Script to create boxplot for iCLIP analysis.
        Created by Anastassiya Zidkova and Martin Zidek, 2014
        """, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-f", "--file", dest="input_file", required=True,
                        help="input file in csv format with two columns, no header. \n Format: <crosslink position> <cDNA count> \n Requred", metavar="file")
    parser.add_argument("-g", "--gtf", dest="gtf_file", required=True,
                        help="input file in gtf format. \n Requred", metavar="file")
    parser.add_argument("-t", "--title", dest="title", default="Boxplot",
                        help="boxplot title", metavar="name")
    parser.add_argument("-n", "--nIntervals", default=None, type=int, dest="n",
                        help="number of intervals in genome", metavar="integer")
    parser.add_argument("-p", "--nParts", default=100, type=int, dest="p",
                        help="number of parts of gene to divide", metavar="integer")

    args = parser.parse_args()
    
    l=get_gene_size(args.gtf_file)
    parts=divide_gene(l[1], args.n)
    interval = determine_gene_intervals(args.gtf_file,l,args.p)   
    a=count_intervals(args.input_file, args.gtf_file, interval, args.p, args.n)
    create_boxplot(a, args.gtf_file, args.title, l, args.n)
