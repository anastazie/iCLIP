# -*- coding: utf-8 -*-
"""
Created on Sun Jul 20 09:43:12 2014

@author: Anastassiya Zidkova and Martin Zidek
"""

import itertools

class Kmers:
    __list__ = []
    
    def __init__(self,n):
        self.current = 0
        x = ['A', 'C', 'T', 'G']
        y = [p for p in itertools.product(x, repeat=n)]
        for t in y:
            self.__list__.append(''.join(t))
        self.high = len(self.__list__) - 1
    
    def __iter__(self):
        return self
    
    def next(self):
        if self.current > self.high:
            raise StopIteration
        else:
            self.current += 1
            return self.__list__[self.current - 1]
            
    def calculate_frequencies(self, seq_list):
        freq_table = []
        
        for pentamere in self.__list__:
            pentamere_record = [pentamere, 0]
            for seq in seq_list:
                if pentamere in seq:
                    pentamere_record[1] += 1
            freq_table.append(pentamere_record)
        return freq_table
    
    @property
    def count(self):
        return len(self.__list__)