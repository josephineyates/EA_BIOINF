# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 14:53:43 2019

@author: josep
"""

################################# N sequence alignment ###################################

# Workflow:
# 1/ do NW to align pairs of sequences 
# 2/ create matrix with distances 
# 3/ create tree linked to distances
# 4/ do a NW for alignments, starting at the most profound leaf of tree

import TD2 
import clustering as cls
import Bio.SeqIO
import networkx as nx
import matplotlib as plt
import TD3

def ouvertureFASTA(fichier):
    handle = open(fichier+".fasta")
    seqs = []
    for seqrec in Bio.SeqIO.parse(handle, "fasta"):
        seqs.append(seqrec)
    handle.close()
    return seqs 
	
sequences = ouvertureFASTA("balibase/RV11.unaligned/BBS11002")
BLOSUM = TD2.extract_file_bla('blosum62.bla')

def print_tree(tree):
    G = nx.Graph()
    for node in  tree:
        G.add_node(node)
    for node in tree:
        if node != tree[node]:
            G.add_edge(tree[node], node)
    nx.draw_networkx(G, prog='dot')
    #plt.show()
    
def reverse_tree(tree):
    rev = {}
    for node in tree:
        parent = tree[node]
        if not parent in rev:
            rev[parent] = []
        rev[parent].append(node)
    return rev

def aux(nom1, nom2, g, e, sequences, tree):
    if not nom1 in tree.keys():
        seq1 = [sequences[int(nom1)]]
    else:
        g1, g2 = tree[nom1][0], tree[nom1][1] 
        seq1 = aux(g1, g2, g, e, sequences, tree)
    if not nom2 in tree.keys():
        seq2 = [sequences[int(nom2)]]
    else:
        g1, g2 = tree[nom2][0], tree[nom2][1] 
        seq2 = aux(g1, g2, g, e, sequences, tree)
    traceback = TD3.NW_affine_multi(seq1,seq2,g,e,TD3.cout_blosum)[1]
    l1,l2 = TD3.affiche_multi(seq1,seq2,traceback)
    alignments = l1+l2
    return alignments

def alignNseq(sequences):
    N = len(sequences)
    g,e = 30,5
    cost_matrix = [[0]*N]*N
    for i in range(N): 
        for j in range(i+1):
            if i!=j:
                cost_matrix[i][j] = TD2.NW_affine(sequences[i].seq,sequences[j].seq,g,e,TD2.cout_blosum)[0]
    table_dist = cls.to_table (cost_matrix)
    tree_NJ = cls.neighbor_joining(table_dist)
    edges = cls.edges_tree(tree_NJ)
    ed,res = cls.find_root(tree_NJ,table_dist)
    cls.place_root(tree_NJ,ed,edges,table_dist,res)
    print_tree(tree_NJ)
    good_tree = reverse_tree(tree_NJ)
    def f(a, b):
        return TD3.cout_blosum(a, b, g=g)
    def aux(nom1, nom2, g, e, sequences, tree):
        if not nom1 in tree.keys():
            seq1 = [sequences[int(nom1)]]
        else:
            g1, g2 = tree[nom1][0], tree[nom1][1] 
            seq1 = aux(g1, g2, g, e, sequences, tree)
        if not nom2 in tree.keys():
            seq2 = [sequences[int(nom2)]]
        else:
            g1, g2 = tree[nom2][0], tree[nom2][1] 
            seq2 = aux(g1, g2, g, e, sequences, tree)
        traceback = TD3.NW_affine_multi(seq1,seq2,g,e,f)[1]
        l1,l2 = TD3.affiche_multi(seq1,seq2,traceback)
        alignments = l1+l2
        return alignments
    #cls.print_gvtree(tree_NJ,table_dist,[seq.name for seq in sequences])
    g1, g2 = good_tree["root"][0], good_tree["root"][1] 
    alignments = aux(g1, g2, g, e, sequences, good_tree)
    
#    alignments = [str(sequences[0].seq)]
#    for i in range(1,N):
#        traceback = TD3.NW_affine_multi(alignments,[sequences[i].seq],g,e,TD2.cout_blosum)[1]
#        l1,l2 = TD3.affiche_multi(alignments,[sequences[i].seq],traceback)
#        alignments = l1+l2
    return alignments

#align = alignNseq(sequences)
#for a in align:
#    print(a)