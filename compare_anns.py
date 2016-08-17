#
# compare_anns.py <true> <pred>
#
# compares a predicted gene structure against the true gene structure and computes
# various statistics summarizing the quality of the prediction. The argument <true> is 
# the true gene structure in faste format, and <pred> is the predicted gene structure 
# in fasta format, e.g.
#
# > python compare_anns.py ./annotation1.fa ./pred1.fa
# > Only Cs (tp=728238, fp=0, tn=505177, fn=249):
# > Sn = 0.9997, Sp = 1.0000, CC = 0.9996, AC = 0.9996
# > Only Rs (tp=618777, fp=0, tn=505426, fn=0):
# > Sn = 1.0000, Sp = 1.0000, CC = 1.0000, AC = 1.0000
# > Both (tp=1347015, fp=0, tn=505177, fn=249):
# > Sn = 0.9998, Sp = 1.0000, CC = 0.9997, AC = 0.9997
# >
#
# Christian Storm <cstorm@birc.au.dk>
#

import os
import sys
import string
import math

def read_ann(filename):
    lines = []
    for l in open(filename).readlines():
        if l[0] != ">":
            lines.append(string.strip(l))
    return string.join(lines, "")

def count_c(true, pred):
    total = tp = fp = tn = fn = 0
    for i in range(len(true)):
        if pred[i] == 'C' or pred[i] == 'c':
            total = total + 1
            if true[i] == 'C' or true[i] == 'c':
                tp = tp + 1
            else:
                fp = fp + 1
        if pred[i] == 'N' or pred[i] == 'n':
            if true[i] == 'N' or true[i] == 'n' or true[i] == 'R' or true[i] == 'r':
                tn = tn + 1
            else:
                fn = fn + 1
    return(total, tp, fp, tn, fn)

def count_r(true, pred):
    total = tp = fp = tn = fn = 0
    for i in range(len(true)):
        if pred[i] == 'R' or pred[i] == 'r':
            total = total + 1
            if true[i] == 'R' or true[i] == 'r':
                tp = tp + 1
            else:
                fp = fp + 1
        if pred[i] == 'N' or pred[i] == 'n':
            if true[i] == 'N' or true[i] == 'n' or true[i] == 'C' or true[i] == 'c':
                tn = tn + 1
            else:
                fn = fn + 1
    return(total, tp, fp, tn, fn)

def count_cr(true, pred):
    total = tp = fp = tn = fn = 0
    for i in range(len(true)):
        if pred[i] == 'C' or pred[i] == 'c' or pred[i] == 'R' or pred[i] == 'r':
            total = total + 1
            if (pred[i] == 'C' or pred[i] == 'c') and (true[i] == 'C' or true[i] == 'c'):
                tp = tp + 1
            elif (pred[i] == 'R' or pred[i] == 'r') and (true[i] == 'R' or true[i] == 'r'):
                tp = tp + 1                
            else:
                fp = fp + 1
        if pred[i] == 'N' or pred[i] == 'n':
            if true[i] == 'N' or true[i] == 'n':
                tn = tn + 1
            else:
                fn = fn + 1
    return(total, tp, fp, tn, fn)

def print_stats(tp, fp, tn, fn):
    sn = float(tp) / (tp + fn)
    sp = float(tp) / (tp + fp)
    cc = float((tp*tn - fp*fn)) / math.sqrt(float((tp+fn)*(tn+fp)*(tp+fp)*(tn+fn)))
    acp = 0.25 * (float(tp)/(tp+fn) + float(tp)/(tp+fp) + float(tn)/(tn+fp) + float(tn)/(tn+fn))
    ac = (acp - 0.5) * 2
    print("Sn = %.4f, Sp = %.4f, CC = %.4f, AC = %.4f" % (sn, sp, cc, ac))

def print_all(true, pred):
    (totalc, tp, fp, tn, fn) = count_c(true, pred)
    if totalc > 0:
        print "Only Cs (tp=%d, fp=%d, tn=%d, fn=%d):" % (tp, fp, tn, fn)
        print_stats(tp, fp, tn, fn)

    (totalr, tp, fp, tn, fn) = count_r(true, pred)
    if totalr > 0:
        print "Only Rs (tp=%d, fp=%d, tn=%d, fn=%d):" % (tp, fp, tn, fn)
        print_stats(tp, fp, tn, fn)

    (total, tp, fp, tn, fn) = count_cr(true, pred)
    if totalc > 0 and totalr > 0:
        print "Both (tp=%d, fp=%d, tn=%d, fn=%d):" % (tp, fp, tn, fn)
        print_stats(tp, fp, tn, fn)
    
# Read true annotation
true_ann = read_ann(sys.argv[1])

# Read predicted annotations
pred_ann = read_ann(sys.argv[2])

# Check annoation length
error = 0
if len(true_ann) != len(pred_ann):
    print("ERROR: The lengths of two predictions are different")
    print("Expected %d, but found %d" % (len(true_ann), len(pred_ann)))
    sys.exit(1)    

# Print stats
print_all(true_ann, pred_ann)


        
    





