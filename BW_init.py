# -*- coding: utf-8 -*-
import numpy as np

"""
Initalize Baum - Welch parameters p, A, B
"""

def init_sim(hid_num, obs_num):
    """
    hid_num is the number of hidden states, in our project only two states, coding(C)
               and non-coding (N). In the reversed sequence, there is also a reversed state(R)
    obs_num is the number of observation states, we have four states (ATCG) in our project.
    """    

## initialize initial probability of hidden states p
    p = np.random.random(hid_num)
    p = p / p.sum()

## initialize transition matrix A
    A = np.random.random((hid_num, hid_num))
    A = A / A.sum(1).reshape((hid_num, 1))

## initialize emission matrix B
    B = np.random.random((hid_num, obs_num))
    B = B / B.sum(1).reshape((hid_num, 1))
    
    return (p, A, B)