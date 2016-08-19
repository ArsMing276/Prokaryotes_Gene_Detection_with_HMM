import numpy as np
import pandas as pd
# -*- coding: utf-8 -*-
"""
Baum Welch parameters calculation

Baum Welch is actually a special form of Expectation - Maximization where we embed 
forward backward into EM to avoid calculating probabilities that require knowledge of 
hidden states. In the while process, we need to calculate the following matrices:

1. Forward probabilities   ----  BW_alpha
2. Backward probabilities  ----  BW_beta
3. Hidden parameter matrix r and Hidden parameter matrix s, these two matrices are
just transformation of the original hidden parameter matrix. 

"""

annotation_map = {'C' : 0, 'N' : 1}
genome_map = {'A' : 0, 'T' : 1, 'C' : 2, 'G' : 3}

def BW_alpha(genomes, A, B, p):
    genomes_train = [x.map(genome_map) for x in genomes]
    D = len(genomes_train)
    T = [len(x) for x in genomes_train]
    
    alphas = []
    for d in range(D):
        obs = genomes_train[d]
        alpha = np.zeros((A.shape[0], T[d]))
        alpha[:, 0] = p * B[:,obs[0]]
        for t in range(1, T[d]):
            for j in range(A.shape[0]):
                alpha[j,t] = sum(alpha[:,(t-1)] * A[:,j] * B[j, obs[t]])
        
            ##scaling to avoid underflow
            alpha[:,t] = alpha[:,t] / sum(alpha[:,t])
        alphas.append(alpha)
    
    return alphas

def BW_beta(genomes, A, B):
    genomes_train = [x.map(genome_map) for x in genomes]
    D = len(genomes_train)
    T = [len(x) for x in genomes_train]
    betas = []
    for d in range(D):
        obs = genomes_train[d]
        beta = np.zeros((A.shape[0], T[d]))
        beta[:, (T[d] - 1)] = 1
        for t in reversed(range(T[d] - 1)):
            for i in range(A.shape[0]):
                beta[i,t] = sum(beta[:, (t+1)] * A[i,:] * B[:, obs[t+1]])
            beta[:,t] = beta[:,t] / sum(beta[:,t])
            
        betas.append(beta)
    
    return betas
    
def BW_r(alphas, betas, A, B):
    D = len(alphas)
    T = [x.shape[1] for x in alphas]      
    
    rs = []
    for d in range(D):
        r = np.zeros((A.shape[0], T[d]))
        alpha = alphas[d]
        beta = betas[d]
        
        for t in range(T[d]):
            r[:,t] = alpha[:,t] * beta[:,t] / sum(alpha[:,t] * beta[:,t])
        
    
        rs.append(r)
    
    return rs
    

def BW_s(genomes, alphas, betas, A, B):
    genomes_train = [x.map(genome_map) for x in genomes]
    D = len(alphas)
    T = [x.shape[1] for x in alphas]      
    
    ss = []
    for d in range(D):
        s = np.zeros((A.shape[0], A.shape[0], (T[d] - 1)))
        alpha = alphas[d]
        beta = betas[d]
        obs = genomes_train[d]
        
        for t in range(T[d] - 1):
            fb = A * B[:, obs[t+1]] * beta[:,(t+1)] * alpha[:,t].reshape((alpha.shape[0], 1))
            s[:,:,t] = fb / sum(sum(fb))
    
        ss.append(s)
    
    return ss
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        