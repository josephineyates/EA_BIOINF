# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:02:26 2019

@author: pepou
"""

from clustering import to_table, neighbor_joining, edges_tree, find_root, place_root
from TD2 import NW_affine, cout_blosum
from Bio.Seq import Seq
from Bio import SeqIO
import os
import time
import networkx as nx
from matplotlib import pyplot as plt

def import_fasta():
    file_list = os.listdir("balibase/RV11.unaligned/")
    s = []
    for file in file_list:
        seqs = []
        handle = open("balibase/RV11.unaligned/" + file)
        for seqrec in SeqIO.parse(handle, "fasta"):
            seqs.append(seqrec)
        s.append(seqs)
    return s

def NW_func(s1, s2):
    return NW_affine(s1, s2, 8, 1, cout_blosum)

def distance_matrix(seqs, fun=NW_func): # fun prend 2 arguments (2 sequences) et renvoie la distance entre les 2
    N = len(seqs)
    t = [[0 for i in range(N)] for j in range(N)]
    for i in range(N):
        for j in range(i):
            d = fun(seqs[i].seq, seqs[j].seq)[0]
            t[i][j] = d
            t[j][i] = d
    return t

def NJ_tree(t):
    table = to_table(t)
    # Construction de l'arbre en utilisant Neighbor-Joining
    start = time.clock()
    tree = neighbor_joining (table)
    # Enracinement de l'arbre
    edges = edges_tree(tree)
    e, res = find_root(tree, table)
    place_root(tree, e, edges, table, res)
    print("\nTIME Neighbor Joining : " + str(time.clock() - start))
    return tree

def print_tree(tree):
    G = nx.Graph()
    for node in  tree:
        G.add_node(node)
    for node in tree:
        if node != tree[node]:
            G.add_edge(tree[node], node)
    nx.draw_networkx(G, prog='dot')
    plt.show()


