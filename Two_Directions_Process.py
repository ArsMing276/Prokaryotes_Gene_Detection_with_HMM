# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:29:58 2016

@author: stanley
"""
def drtprocess(genomes, annotations, reverse = False):
    #if reversed, change C to N and reverse all genomes and annotations sequences
    #otherwise, only change R to N
    if reverse:
        genomes = [x[::-1] for x in genomes]
        annotations = [x[::-1] for x in annotations]
        
        annotations = [map(lambda x: 'N' if x == 'C' else x, y) for y in annotations]
        annotations = [map(lambda x: 'C' if x == 'R' else x, y) for y in annotations]
    else:
        annotations = [map(lambda x: 'N' if x == 'R' else x, y) for y in annotations]
    
    return (genomes, annotations)