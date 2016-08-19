# -*- coding: utf-8 -*-
"""
After Baum Welch training, we use viterbi to find the most possible hidden sequences for 
the first five genome sequences, and evaluate the accruacy
"""
import numpy as np
import pandas as pd

##suppose we have already have the converged initial probbability, transition matrix and emission matrix
annotation_map = {'C' : 0, 'N' : 1}
genome_map = {'A' : 0, 'T' : 1, 'C' : 2, 'G' : 3}

genomes_test = genomes[0:len(annotations)]
genomes_test = [x.map(genome_map) for x in genomes_test]
annotations = [x.map(annotation_map) for x in annotations]

def viterbi(genomes, A, B, p):
    D = len(genomes)
    T = [len(x) for x in genomes]
    
    sequences = []
    for d in range(D):
        obs = genomes[d]
        
        seqs = np.zeros((A.shape[0], T[d]))
        traces = np.zeros((A.shape[0], T[d]), dtype = np.int)
        best_seq = np.zeros(T[d], dtype = np.int)
        
        ##to avoid underflow, we work in log space 
        seqs[:, 0] = np.log(p * B[:,obs[0]])
        traces[:, 0] = [0,1]
        for t in range(1, T[d]):
            for j in range(A.shape[0]):
                seqs[j,t] = (seqs[:,(t-1)] + np.log(A[:,j]) + np.log(B[j, obs[t]])).max()
                traces[j,t] = (seqs[:,(t-1)] + np.log(A[:,j]) + np.log(B[j, obs[t]])).argmax()
        
        ##back-tracking to find the best hidden sequence
        best_seq[T[d] - 1] = (seqs[:,T[d] - 1]).argmax()
        for i in reversed(range(T[d] - 1)):
            from_idx = traces[:,(i+1)][best_seq[i+1]]
            best_seq[i] = from_idx

        sequences.append(best_seq)
    
    return sequences


def main():
    annotations_test = viterbi(genomes_test, A, B, p)
    precision = np.mean([np.mean(x == y) for x,y in zip(annotations_test, annotations)])  
    
    print('Precision of Hidden States prediction for the first five genomes is %.2f' %precision)
    return annotations_test
    
if __name__ == "__main__":
    main()