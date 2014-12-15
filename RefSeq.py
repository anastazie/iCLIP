# -*- coding: utf-8 -*-
"""
Created on Sun Jul 20 09:59:32 2014

@author: nastya
"""

class RefSeq:
    def __init__(self, filename):
        in_file = open(filename,'r')
        line = in_file.readline()
        data = line.strip().split('\t')
        
        self.__sequence__ = data[1].upper()
    
    def get_position(self, position, length):
        if ((position + length) > len(self.__sequence__)) | ((position - length) < 0):
            raise Exception('Requested position out of reference data set boundaries')
        else:
            subseq = self.__sequence__[position - length:position + length + 1]
            return subseq
            