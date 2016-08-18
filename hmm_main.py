# -*- coding: utf-8 -*-
"""
Main function
"""
import numpy as np
import pandas as pd
import glob
import re
from Two_Directions_Process import drtprocess
from BW_init import init_sim
from BW_pars import BW_alpha, BW_beta, BW_r, BW_s, annotation_map, genome_map
from BW_update import update


def read_genomes(path):
    genomes = []
    
    files = glob.glob(path + '/genome*.txt')
    sortfun = lambda x: int(re.findall('genome(.*).txt', x).pop())
    files = sorted(files, key = sortfun)
     
    for file in files:
        genome = open(file).readlines()
        genome = [x.strip('\n') for x in genome]
        genomes.append(list(''.join(genome)))
    
    return genomes

def read_annotation(path):
    annotations = []
    
    files = glob.glob(path + '/annotation*.txt')
    sortfun = lambda x: int(re.findall('annotation(.*).txt', x).pop())
    files = sorted(files, key = sortfun)
     
    for file in files:
        annotation = open(file).readlines()
        annotation = [x.strip('\n') for x in annotation]
        annotations.append(list(''.join(annotation)))
    
    return annotations

def main():
    genomes = read_genomes('./data')
    annotations = read_annotation('./data')
    
    ##change False to True when finding reversed genes
    genomes, annotations = drtprocess(genomes, annotations, False)
    genomes = [pd.Series(x) for x in genomes]
    annotations = [pd.Series(x) for x in annotations]
    
    p0, A0, B0 = init_sim(len(np.unique(annotations[0])), len(np.unique(genomes[0])))
    while True:
        alphas = BW_alpha(genomes, A0, B0, p0)
        betas = BW_beta(genomes, A0, B0)
        rs = BW_r(alphas, betas, A0, B0)
        ss = BW_s(genomes, alphas, betas, A0, B0)
        p, A, B = update(genomes, rs, ss)
        if (abs(B - B0) < 0.01).all() and (abs(A - A0) < 0.01).all() and (abs(p - p0) < 0.01).all():
            break
        else:
            p0 = p
            A0 = A
            B0 = B
    
    return (p, A, B)
    
    
    
if __name__ == "__main__":
    main()