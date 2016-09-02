# Prokaryotes_Gene_Detection_with_HMM

This project is about using a hidden Markov model for gene finding in prokaryotes.

We have a data set containing 11 Staphylococcus genomes, each containing several 
genes (i.e. substring) obeying the "gene syntax". The genomes are between 1.8 million 
and 2.8 million nucleotides. For 5 of the genomes, we also know the location of 
the genes. For the remaining 6 genomes, we only know that they contain genes according 
to the "gene syntax". The genomes and their annontations are given in FASTA format. 

In this project, we will train a hidden markov model on the 11 genomes with **Baum-Welch 
Algorithm**. After that, we will use **Viterbi Algorithm** to infer the most possible hidden states
(i.e, genes or not) of the fist five genomes. Finally we will compare the predictions with the true 
gene coding and calculate the precision of the prediction.

Something needs our attention...

1. Multiple start-codons are possible in each genome which means we shouldn't fix any probability when training HMM model. In general, if we know what hidden states are or we know some probabilities should be fixed, we could train HMM *By Counting*, otherwise Baum-Welch algorithm is necessary in training. 

2. To predict the full gene structure, we face the fact that the genomes have genes in both directions (i.e. a nucleotide can be C or R). We can make a HMM which only models genes in one direction and then use it twice to predict the genes in each direction or we can make a model which models genes in both directions.

3. Since genomes are long, i.e. running Viterbi and Baum-Welch takes a lot of memory and time. Thus this implementation is far from perfect and need further improvement.


# References
http://cs.au.dk/~cstorm/courses/MLiB_f14/project3.html, Department of Computer Science, AARHUS UNIVERSITY
