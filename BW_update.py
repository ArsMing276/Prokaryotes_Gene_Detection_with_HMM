# -*- coding: utf-8 -*-
"""
This is the following-up step of parameter matrices calculation. After we have updated 
alpha, beta, r and s, we could update p , A and B to finish a iteration.
"""
import numpy as np
import pandas as pd

annotation_map = {'C' : 0, 'N' : 1}
genome_map = {'A' : 0, 'T' : 1, 'C' : 2, 'G' : 3}

def update(genomes, rs, ss):
    genomes_train = [x.map(genome_map) for x in genomes]
    D = len(genomes_train)
    
    ##update p
    rs0 = [x[:,0] for x in rs]
    p = np.add.reduce(rs0) / D

    ##update A
    ssum = sum([x.sum(2) for x in ss])
    rs0 = [x[:, :-1] for x in rs]
    rsum = sum([x.sum(1) for x in rs0])
    A = ssum / rsum.reshape((rsum.shape[0], 1))
    
    ##update B
    obs_states = np.unique(genomes_train[0])
    n = len(obs_states)
    
    B = np.zeros((rs[0].shape[0],n))
    for j in range(n):
        idx = [x.index[x == obs_states[j]] for x in genomes_train]
        rcond = [x[:,y] for x,y in zip(rs, idx)]
        Bnom = sum([x.sum(1) for x in rcond])
        Bdenom = sum([x.sum(1) for x in rs])
        B[:,j] = Bnom / Bdenom
    
    return (p, A, B)