# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 14:17:52 2019

@author: pepou
"""

import os
import TD3
import clustering as cls
from Bio.SeqRecord import SeqRecord
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import ExposureCN, HSExposureCB
from TD3 import cout_blosum
from NW_align_N_seq import print_tree, reverse_tree
from pdb_analysis import saveFASTA, decode_enf

def get_atoms(file):
    parser = MMCIFParser()
    structure = parser.get_structure(file.split('.')[0],file)
    pos = []
    model = structure[0]
    for chain in model:
        pos_c = []
        for residue in chain:
            if residue.has_id('CA'):
                vca = residue['CA'].get_vector()
                pos_c.append((residue.get_resname(), vca))
        pos.append(pos_c)
    return pos

def get_descriptors(file):
    parser = MMCIFParser()
    structure = parser.get_structure(file.split('.')[0],file)
    pos = []
    model = structure[0]
    hse = HSExposureCB(model)
    for chain in model:
        pos_c = []
        for residue in chain:
            dic = {}
            dic["name"] = residue.get_resname()
            if residue.has_id('CA'):
                vca = residue['CA'].get_vector()
                dic["coord"] = vca
                hse_ = hse[(chain.id, residue.id)]
                dic["hse"] = (hse_[0], hse_[1])
            pos_c.append(dic)
        pos = pos + pos_c
    return pos

def cout_structural(a, ind_a, b, ind_b, desc_a, desc_b):
    c = cout_blosum(a, b)
    # Différence d'enfouissement
    if "enf" in desc_a and "enf" in desc_b:
        c -= abs(desc_a["enf"] - desc_b["enf"]) // 10
    # Pénalité HSExposure
    if "hse" in desc_a and "hse" in desc_b:
        c -= (abs(desc_a["hse"][0] - desc_b["hse"][0]) + abs(desc_a["hse"][1] - desc_b["hse"][1])) // 4
    return c

def NW_affine_structure(seq1, seq2, g, e, cout, desc): # coût est une fonction, g et e doit être positif
    n = len(seq1)
    m = len(seq2)
    M = [[0 for i in range(m+1)] for j in range(n+1)]
    Ix = [[0 for i in range(m+1)] for j in range(n+1)]
    Iy = [[0 for i in range(m+1)] for j in range(n+1)]
    tM = [[0 for i in range(m+1)] for j in range(n+1)]
    tIx = [[0 for i in range(m+1)] for j in range(n+1)]
    tIy = [[0 for i in range(m+1)] for j in range(n+1)]
    
    for i in range(1, n+1):
        M[i][0] = -g+(1-i)*e
        Ix[i][0] = -g+(1-i)*e
        Iy[i][0] = -g+(1-i)*e
        tM[i][0] = 1 #up
        tIx[i][0] = 1 #up
        tIy[i][0] = 1 #up
    for j in range(1, m+1):
        M[0][j] = -g+(1-j)*e
        Ix[0][j] = -g+(1-j)*e
        Iy[0][j] = -g+(1-j)*e
        tM[0][j] = 2 #left
        tIx[0][j] = 2 #left
        tIy[0][j] = 2 #left
        
    for i in range(1, n+1):
        for j in range(1, m+1):
            #M
            cout_ab = cout(seq1[i-1], i-1, seq2[j-1], j-1, desc[0], desc[1])
            m_ = M[i-1][j-1] + cout_ab
            ix = Ix[i-1][j-1] + cout_ab
            iy = Iy[i-1][j-1] + cout_ab
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
            m_ = M[i-1][j] - g
            ix = Ix[i-1][j] - e
            if m_ > ix:
                Ix[i][j] = m_
                tIx[i][j] = 1
            else:
                Ix[i][j] = ix
                tIx[i][j] = 6 # up dans Ix -1, 0
            #Iy
            m_ = M[i][j-1] - g
            iy = Iy[i][j-1] - e
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
            out = "-" + out
            i -= 1
            j -= 1
            pos = tM
        elif pos[i][j] == 4:
            out = "-" + out
            i -= 1
            j -= 1
            pos = tIx
        elif pos[i][j] == 5:
            out = "-" + out
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

def NW_affine_multi_structure(seq1, seq2, g, e, cout, desc): # coût est une fonction, g et e doit être positif
    n = len(seq1[0])
    m = len(seq2[0])
    nb1 = len(seq1)
    nb2 = len(seq2)
    M = [[0 for i in range(m+1)] for j in range(n+1)]
    Ix = [[0 for i in range(m+1)] for j in range(n+1)]
    Iy = [[0 for i in range(m+1)] for j in range(n+1)]
    tM = [[0 for i in range(m+1)] for j in range(n+1)]
    tIx = [[0 for i in range(m+1)] for j in range(n+1)]
    tIy = [[0 for i in range(m+1)] for j in range(n+1)]
    
    for i in range(1, n+1):
        M[i][0] = -g+(1-i)*e
        Ix[i][0] = -g+(1-i)*e
        Iy[i][0] = -g+(1-i)*e
        tM[i][0] = 1 #up
        tIx[i][0] = 1 #up
        tIy[i][0] = 1 #up
    for j in range(1, m+1):
        M[0][j] = -g+(1-j)*e
        Ix[0][j] = -g+(1-j)*e
        Iy[0][j] = -g+(1-j)*e
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
                    cout_ab = cout(seq1[index1][i-1], i-1, seq2[index2][j-1], j-1, desc[0][index1], desc[1][index2])
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
            m_ = M[i-1][j] - g
            ix = Ix[i-1][j] - e
            if m_ > ix:
                Ix[i][j] = m_
                tIx[i][j] = 1
            else:
                Ix[i][j] = ix
                tIx[i][j] = 6 # up dans Ix -1, 0
            #Iy
            m_ = M[i][j-1] - g
            iy = Iy[i][j-1] - e
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

def alignNseq_structure(sequences, desc):
    N = len(sequences)
    g,e = 50,30
    cost_matrix = [[0]*N]*N
    for i in range(N): 
        for j in range(i+1):
            if i!=j:
                cost_matrix[i][j] = NW_affine_structure(sequences[i].seq,sequences[j].seq,g,e,cout_structural, [desc[i], desc[j]])[0]
    table_dist = cls.to_table (cost_matrix)
    tree_NJ = cls.neighbor_joining(table_dist)
    edges = cls.edges_tree(tree_NJ)
    ed,res = cls.find_root(tree_NJ,table_dist)
    cls.place_root(tree_NJ,ed,edges,table_dist,res)
    print_tree(tree_NJ)
    good_tree = reverse_tree(tree_NJ) # good_tree stocke la liste des fils de chaque noeud s'il n'est pas une feuille
    def aux(nom1, nom2, g, e, sequences, tree):
        if not nom1 in tree.keys(): # si nom1 est une feuille
            seq1 = [sequences[int(nom1)]]
            desc1 = [desc[int(nom1)]]
        else:
            g1, g2 = tree[nom1][0], tree[nom1][1] 
            seq1, desc1 = aux(g1, g2, g, e, sequences, tree)
        if not nom2 in tree.keys():
            seq2 = [sequences[int(nom2)]]
            desc2 = [desc[int(nom2)]]
        else:
            g1, g2 = tree[nom2][0], tree[nom2][1] 
            seq2, desc2 = aux(g1, g2, g, e, sequences, tree)
        traceback = NW_affine_multi_structure(seq1,seq2,g,e,cout_structural, desc)[1]
        l1,l2 = TD3.affiche_multi(seq1,seq2,traceback)
        alignments = l1+l2
        new_desc = desc1 + desc2
        return alignments, new_desc
    #cls.print_gvtree(tree_NJ,table_dist,[seq.name for seq in sequences])
    g1, g2 = good_tree["root"][0], good_tree["root"][1] 
    alignments, _ = aux(g1, g2, g, e, sequences, good_tree)
    
    return alignments

def test():
    rec = saveFASTA("balibase/RV11.unaligned/BBS11002.fasta")
    descriptions = []
    for r in rec:
        desc = get_descriptors("PDB/"+r.name[:4]+".cif")
        values = decode_enf(r.name[:4])
        assert len(desc) == len(values), "Pas bon r={} : {} vs {}".format(r.name[:4], len(desc), len(values))
        for k in range(len(desc)):
            desc[k]["enf"] = values[k]
        descriptions.append(desc)        
    aln = alignNseq_structure(rec, descriptions)
    return aln
                   