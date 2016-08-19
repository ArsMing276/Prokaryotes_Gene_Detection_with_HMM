# -*- coding: utf-8 -*-
"""
Each genome has two sides, they may contain different genes at differnet location. We
need to consider the two direction separately. This script helps reverse genome if necessary.
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