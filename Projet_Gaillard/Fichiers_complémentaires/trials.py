# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 16:11:34 2019

@author: pepou
"""

from methods import *

BLOSUM = extract_file_bla("blosum62.bla")

def cout_blosum(seq1, i, seq2, j, g=lambda x : 11, e=lambda x : 1, mat=BLOSUM): #mat doit être la matrice BLOSUM
    a = seq1["seq"][i]
    b = seq2["seq"][j]
    if a == "-":
        if b == "-":
            return 0
        if j > 1 and seq2["seq"][j-1] == "-":
            return -e(b)
        return -g(b)
    if b == "-":
        if i > 1 and seq1["seq"][i-1] == "-":
            return -e(a)
        return -g(a)
    return mat[a][b]
    
def cout_bidon(seq1, i, seq2, j):
    a = seq1["seq"][i]
    b = seq2["seq"][j]
    if a==b:
        return 8
    return -2

def cout_enf(seq1, i, seq2, j):
    if "enf" in seq1 and "enf" in seq2:
        c = -abs(seq1["enf"][i] - seq2["enf"][j])
        return c
    return 0

def cout_hse(seq1, i, seq2, j):
    if "hse" in seq1 and "hse" in seq2 and seq1["hse"][i] != None and seq2["hse"][j] != None:
        c = -(abs(seq1["hse"][i][0] - seq2["hse"][j][0]) + abs(seq1["hse"][i][1] - seq2["hse"][j][1]))
        return c
    return 0

class NeedlemanWunsch:
    ''' Class for Needleman-Wunsch algorithm
        Can be initialized with different parameters :
            - cost_functions : list of cost functions that applies when comparing two residues
            - cost_coefs : list of coefficients to multiply the previous cost functions. Sum must be 1
            - gap_opening_function and gap_extending_function : how to penalize the alignment between a residue and a gap
            - clustering : method for clustering, ie "NJ" for Neighbor Joining or "UPGMA"
    '''
    
    def __init__(self, cost_functions=cout_blosum, gap_opening_function=11, gap_extending_function=1, clustering="NJ", cost_coefs=None):
        if isinstance(cost_functions, list):
            N = len(cost_functions)
            if cost_coefs == None:
                cost_coefs = [1./float(N)] * N
            else:
                assert len(cost_coefs) == N,  \
                "{} coefs for cost functions are given, whereas {} cost functions are given".format(len(cost_coefs), N)
            self.cost_coefs = cost_coefs
        else:
            self.cost_coefs = [1.]
            cost_functions = [cost_functions]
        self.costs = cost_functions
        if isinstance(gap_opening_function, int):
            gap_ = lambda x, y : gap_opening_function
            self.gap = gap_
        else:
            self.gap = gap_opening_function
        if isinstance(gap_extending_function, int):
            extend_ = lambda x, y : gap_extending_function
            self.extend = extend_
        else:
            self.extend = gap_extending_function
        self.clustering = clustering
        
    def align(self, sequences1, sequences2):
        ''' Aligns two groups of sequences already aligned among each group
            sequences are list of dictionaries containing the sequence and sometimes descriptors '''
        n = len(sequences1[0]["seq"])
        m = len(sequences2[0]["seq"])
        nb1 = len(sequences1)
        nb2 = len(sequences2)
        M = [[0 for j in range(m+1)] for i in range(n+1)]
        Ix = [[0 for j in range(m+1)] for i in range(n+1)]
        Iy = [[0 for j in range(m+1)] for i in range(n+1)]
        tM = [[0 for j in range(m+1)] for i in range(n+1)]
        tIx = [[0 for j in range(m+1)] for i in range(n+1)]
        tIy = [[0 for j in range(m+1)] for i in range(n+1)]

        M[1][0] = -sum([self.gap(sequences1[j], 0) for j in range(nb1)])
        Ix[1][0] = M[1][0]
        Iy[1][0] = M[1][0]
        tM[1][0] = 1 #up
        tIx[1][0] = 1 #up
        tIy[1][0] = 1 #up    
        for i in range(2, n+1):
            M[i][0] = M[i-1][0] - sum([self.extend(sequences1[j],i-1) for j in range(nb1)])
            Ix[i][0] = M[i][0]
            Iy[i][0] = M[i][0]
            tM[i][0] = 1 #up
            tIx[i][0] = 1 #up
            tIy[i][0] = 1 #up
            
        M[0][1] = -sum([self.gap(sequences2[j], 0) for j in range(nb2)])
        Ix[0][1] = M[0][1]
        Iy[0][1] = M[0][1]
        tM[0][1] = 2 #left
        tIx[0][1] = 2 #left
        tIy[0][1] = 2 #left
        for j in range(2, m+1):
            M[0][j] = M[0][j-1] - sum([self.extend(sequences2[i], j-1) for i in range(nb2)])
            Ix[0][j] = M[0][j]
            Iy[0][j] = M[0][j]
            tM[0][j] = 2 #left
            tIx[0][j] = 2 #left
            tIy[0][j] = 2 #left
        
        for i in range(1, n+1):
            for j in range(1, m+1):
                #M
                m_ = M[i-1][j-1]
                ix = Ix[i-1][j-1]
                iy = Iy[i-1][j-1]
            
                for index1 in range (nb1):
                    for index2 in range (nb2):
                        cout_ab = sum([cout(sequences1[index1], i-1, sequences2[index2], j-1) for (cout, coef) in zip(self.costs, self.cost_coefs)])
                        m_ += cout_ab
                        ix += cout_ab
                        iy += cout_ab
            
                if ix > iy:
                    if ix > m_:
                        M[i][j] = ix
                        tM[i][j] = 4 # dans Ix -1,-1
                    else:
                        M[i][j] = m_
                        tM[i][j] = 3 #diag
                else:
                    if iy > m_:
                        M[i][j] = iy
                        tM[i][j] = 5 # dans Iy -1,-1
                    else:
                        M[i][j] = m_
                        tM[i][j] = 3 #diag
                #Ix
                m_ = M[i-1][j] - sum([self.gap(sequences2[k], j-1) for k in range(nb2)])
                ix = Ix[i-1][j] - sum([self.extend(sequences2[k], j-1) for k in range(nb2)])
                if m_ > ix:
                    Ix[i][j] = m_
                    tIx[i][j] = 1
                else:
                    Ix[i][j] = ix
                    tIx[i][j] = 6 # up dans Ix -1, 0
                #Iy
                m_ = M[i][j-1] - sum([self.gap(sequences1[k], i-1) for k in range(nb1)])
                iy = Iy[i][j-1] - sum([self.extend(sequences1[k], i-1) for k in range(nb1)])
                if m_ > iy:
                    Iy[i][j] = m_
                    tIy[i][j] = 2
                else:
                    Iy[i][j] = iy
                    tIy[i][j] = 7 # left dans Iy 0, -1    
    
        # Now the traceback :
        i = n
        j = m
        out = ""
        pos = None
        if Ix[i][j] > Iy[i][j]:
            if Ix[i][j] > M[i][j]:
                s = Ix[i][j]
                pos = tIx
            else:
                s = M[i][j]
                pos = tM
        else:
            if Iy[i][j] > M[i][j]:
                s = Iy[i][j]
                pos = tIy
            else:
                s = M[i][j]
                pos = tM
        while (i + j > 0):
            if pos[i][j] == 1:
                out = "1" + out
                i -= 1
                pos = tM
            elif pos[i][j] == 2:
                out = "2" + out
                j -= 1
                pos = tM
            elif pos[i][j] == 3:
                out = "*" + out
                i -= 1
                j -= 1
                pos = tM
            elif pos[i][j] == 4:
                out = "*" + out
                i -= 1
                j -= 1
                pos = tIx
            elif pos[i][j] == 5:
                out = "*" + out
                i -= 1
                j -= 1
                pos = tIy
            elif pos[i][j] == 6:
                out = "1" + out
                i -= 1
                pos = tIx
            else:
                out = "2" + out
                j -= 1
                pos = tIy
        return s, out
    
    @staticmethod
    def affiche_multi(sequences1, sequences2, traceback):
        L = len(traceback)
        nb1 = len(sequences1)
        nb2 = len(sequences2)
        l1 = [""]*nb1
        l2 = [""]*nb2
        i, j = 0, 0
        for k in range(L):
            if traceback[k] == "1":
                for index1 in range(nb1):
                    l1[index1] = l1[index1] + sequences1[index1]["seq"][i]
                i += 1
                for index2 in range(nb2):
                    l2[index2] = l2[index2] + "-"
            elif traceback[k] == "2":
                for index1 in range(nb1):
                    l1[index1] = l1[index1] + "-"
                for index2 in range(nb2):
                    l2[index2] = l2[index2] + sequences2[index2]["seq"][j]
                j += 1
            else:
                for index1 in range(nb1):
                    l1[index1] = l1[index1] + sequences1[index1]["seq"][i]
                i += 1
                for index2 in range(nb2):
                    l2[index2] = l2[index2] + sequences2[index2]["seq"][j]
                j += 1
        for index1 in range(nb1):
            sequences1[index1]["seq"] = l1[index1]
        for index2 in range(nb2):
            sequences2[index2]["seq"] = l2[index2]
    
    def tree_align(self, nom1, nom2, sequences, tree):
            if not nom1 in tree.keys(): # si nom1 est une feuille
                seq1 = [sequences[int(nom1)]]
            else:
                g1, g2 = tree[nom1][0], tree[nom1][1] 
                seq1 = self.tree_align(g1, g2, sequences, tree)
            if not nom2 in tree.keys():
                seq2 = [sequences[int(nom2)]]
            else:
                g1, g2 = tree[nom2][0], tree[nom2][1] 
                seq2 = self.tree_align(g1, g2, sequences, tree)
            traceback = self.align(seq1,seq2)[1]
            self.affiche_multi(seq1,seq2,traceback)  
            return seq1 + seq2
    
    def run(self, sequences):
        
        N = len(sequences)
        
        # Computing the distance matrix between all sequences
        cost_matrix = [[0]*N]*N
        for i in range(N): 
            for j in range(i+1):
                if i!=j:
                    cost_matrix[i][j] = self.align([sequences[i]], [sequences[j]])[0]
        table_dist = to_table (cost_matrix)
        
        #Clustering
        if self.clustering == "NJ":
            tree_NJ = neighbor_joining(table_dist)
            edges = edges_tree(tree_NJ)
            ed,res = find_root(tree_NJ,table_dist)
            place_root(tree_NJ,ed,edges,table_dist,res)
            print_tree(tree_NJ)
            good_tree = reverse_tree(tree_NJ) # good_tree stocke la liste des fils de chaque noeud s'il n'est pas une feuille
        elif self.clustering == "UPGMA":
            tree_UPGMA = tree_upgma(table_dist)
            print_tree(tree_UPGMA)
            good_tree = reverse_tree(tree_UPGMA)
        
        g1, g2 = good_tree["root"][0], good_tree["root"][1] 
        return self.tree_align(g1, g2, sequences, good_tree)
    
def test():
    rec = readFASTA("balibase/RV11.unaligned/BBS11001.fasta")
    sequences = []
    for r in rec:
        seq = {"seq" : r.seq}
        desc = get_descriptors("PDB/"+r.name[:4]+".cif")
        seq["hse"] = desc["hse"]
        values = decode_enf(r.name[:4])
        seq["name"] = r.name[:4]
        seq["enf"] = values
        assert len(desc["hse"]) == len(values), "Pas bon r={} : {} vs {}".format(r.name[:4], len(desc), len(values))
        sequences.append(seq)
    NW = NeedlemanWunsch()        
    sequences = NW.run(sequences)
    aln = [seq["seq"] for seq in sequences]
    return aln
        