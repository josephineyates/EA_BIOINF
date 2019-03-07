# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:04:50 2019

@author: josep
"""

#################################################### README ################################################################################
""" 
File structure :
    0/ Utils
        A/ BLOSUM functions
        B/ PDB analysis functions
        C/ Structure information retrieval
        D/ Clustering functions
    I/ 2 sequence alignment with linear Needleman-Wunsch
    II/ N sequence alignment with complex linear Needleman-Wunsch
    III/ Scoring of alignment quality
    IV/ N sequence alignment with complex linear Needleman-Wunsch using structural information

Method list : 
    I/ 
    
    II/ 
    
    III/
    
    IV/
    
To obtain alignments with corresponding scores, use main.py [README for more information]
"""
##################################################### 0/ Utils ###############################################################################

######################################################### A/ BLOSUM functions ###################################################################

def extract_file_bla(filename):
    ''' Extracts blosum from file 'blosum35.bla' with corresponding AA sequence
        Returns list shape(n,n) '''
    table = []
    file = open(filename, "r")
    for i in range(6):
        file.readline()
    aa = []
    l = file.readline()
    for x in l[1:-1].split("  "):
        if x != "\n":
            aa.append(x)
    n = len(aa)
    for i in range(n):
        l = file.readline()[1:-1].split("  ")
        L = []
        for k in range(len(l)):
            L = L + l[k].split(" ")
        for j in range(n):
            L[j] = int(L[j])
        table.append(L[:-1])
    t = {}
    for i in range(n):
        t[aa[i]] = {}
        for j in range(n):
            t[aa[i]][aa[j]] = table[i][j]
    return t

# Initialisation of BLOSUM matrix gloablly, so that it only needs to be charged once   
BLOSUM = {}
def blosum():
    ''' Initialisation of BLOSUM matrix from file blosum62.bla if it's not already done '''
    global BLOSUM
    if BLOSUM == {}:
        BLOSUM = extract_file_bla('blosum62.bla')
    
def cout_blosum(a, b, mat=None): 
    """ returns BLOSUM cost associated to AA a and b """
    if mat == None:
        blosum()
        mat = BLOSUM
    return mat[a][b]

############################################################ B/ PDB analysis functions #####################################################"
    
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import numpy as np

def readFASTA(doc):
    ''' Reads a fasta file '''
    records = []
    for seq_record in SeqIO.parse(doc, "fasta"):
        records.append(seq_record)
    return records

def ouvertureFASTA(fichier):
    ''' Same as above '''
    handle = open(fichier+".fasta")
    seqs = []
    for seqrec in Bio.SeqIO.parse(handle, "fasta"):
        seqs.append(seqrec)
    handle.close()
    return seqs

def download_PDB_file():
    ''' Downloads all the PDB files related to the balibase list of alignments'''
    fichiers = os.listdir("balibase/RV11.unaligned")
    for file in fichiers:
        records = readFASTA("balibase/RV11.unaligned/"+file)
        ids=[]
        for record in records:
            ids.append(record.id.split("_")[0])
        for i in ids:
            pdbl = PDBList()
            pdbl.retrieve_pdb_file(i,pdir="PDB")
              
def change_PDB_files():
    ''' Puts all the files in the same directory '''
    noms_dir = os.listdir("PDB")
    for n in noms_dir:
        fichiers = os.listdir(n)
        for fich in fichiers:
            os.rename("PDB/"+n+"/"+fich,"PDB/"+fich)

############################################################### C/ Structure information retrieval ####################################################
        
def get_info_mmcif(file):
    ''' Computes the list of barycenters for each chain in the structure stored in file which must be a .cif
        Then computes the relative depth of a residue as the distance to the barycenter '''
    parser = MMCIFParser()
    structure = parser.get_structure(file.split('.')[0],file)
    coord_ca={}
    bary={}
    for chain in structure[0]:
        coord_ca[chain]=[]
        bary[chain]=0
        for residue in chain:
            if residue.has_id('CA'):
                coord_ca[chain].append(residue['CA'].get_coord())
            else:
                coord_moy=[0,0,0]
                for atom in residue:
                    coord_at=atom.get_coord()
                    coord_moy=[coord_at[i]/len(residue) for i in range(3)]
                coord_ca[chain].append(coord_moy)
        coord_ca[chain]=np.asarray(coord_ca[chain])
        bary[chain] = np.array([np.mean(coord_ca[chain][i]) for i in range(3)])
    enf={}
    for chain in structure[0]:
        enf[chain]=[]
        for coord in coord_ca[chain]:
            enf[chain].append(np.linalg.norm(coord-bary[chain]))
    #ppb = PPBuilder()
    #seqpdb = ppb.build_peptides(chain)[0].get_sequence()
    return bary,enf

def get_enf():
    ''' Saves all the barycenters and depth of all the files stored in the PDB folder '''
    fichiers = os.listdir("PDB")
    for fich in fichiers:
        bary,enf=get_info_mmcif("PDB/"+fich)
        with open("enfouissements.txt","a") as f:
            f.write("new_file \n")
            f.write(fich+"\n")
            for chain in enf.keys():
                f.write("new_chain \n")
                f.write(str(chain)+"\n")
                for i in range(len(enf[chain])):
                    f.write(str(i)+","+str(enf[chain][i])+"\n")
        f.close()
        
def decode_enf(file_id):
    ''' Gets all the values of the residue's depth 
        Returns the list of "enf" values in the sequential order '''
    file = open("enfouissements.txt", "r")
    content = file.read()
    lines = content.split("\n")[:-1]
    N = len(lines)
    good_file = False
    next_file = False
    values = []
    for k in range(N):
        if next_file:
            if lines[k].split(".cif")[0] == file_id:
                good_file = True
            else:
                good_file = False
            next_file = False
        elif lines[k] == "new_file ":
            next_file = True
        elif not "hain" in lines[k] and good_file:
            tokens = lines[k].split(",")
            values.append(float(tokens[1]))
    return values

############################################################## D/ Clustering functions ######################################################"

import networkx as nx
import graphviz as gv
import matplotlib.pyplot as plt

def to_table(t):
    ''' Pour le rendre utilisable on transforme la matrice des distances 
    en une table à deux clés qui sont des String d'entiers (plus facilement manipulables par
    l'algorithme Neighbor Joining.'''
    table = {}
    for i in range(len(t)):
        l = {}
        for j in range(len(t)):
            l[str(j)] = t[i][j]
        table[str(i)] = l
    return table

def mini(table):
    ''' Calcule le minimum de la table delta dans Neighbor Joining et renvoie le couple des indices correspondants '''
    mini = None
    i_ = None
    j_ = None
    for i in table:
        if table[i][0]:
            for j in table[i][1]:
                if i != j:
                    if table[j][0]:
                        if mini == None:
                            mini = table[i][1][j]
                            i_ = i
                            j_ = j
                        if table[i][1][j]<mini :
                            mini =table[i][1][j]
                            i_=i
                            j_=j
    return (i_, j_)


def neighbor_joining(t):
    ''' Implémentation de l'algorithme Neighbor Joining à partir d'une table des distances.
    Renvoie un arbre non-enraciné tree où tree(i) = j indique que i pointe vers j. Nous verrons par la suite
    que cette relation devra parfois être inversée afin d'enraciner l'arbre.'''
    
    n = len(t)
    tree = {}
    D = {}
    delta = {}
    # Initialisation de l'arbre : au début tout le monde pointe vers lui-même
    for i in t:
        tree[i] = i
    for i in t:
        delta[i] = [True, {}]
        # True veut dire que i est encore à traiter. lorsque i et j fusionnent, on a delta[i][0] = False
    
    # On effectue la boucle n-2 afin de faire fusionner les taxons en deux arbres binaires
    for k in range(n, 2, -1):
        
        # Calcul de la pseudo-distance delta à minimiser entre tous les noeuds du graphe
        for i in t:
            # Calcul de D : représente intuitivement la distance moyenne d'un noeud aux autres noeuds
            if delta[i][0]: # si i est encore à traiter
                s = 0.
                for j in t[i]: 
                    if delta[j][0]: # de même pour j
                        s += t[i][j]
                D[i] = s / (k - 2)
        for i in t:
            if delta[i][0]:
                for j in t:
                    if i!= j and delta[j][0]:
                        delta[i][1][j] = t[i][j] - D[i] - D[j]
                        
        # on extrait les deux noeuds les plus proches au sens de la distance delta
        (i, j) = mini(delta)
        
        # on crée leur père
        ij = i +" "+ j
        tree[ij] = ij
        tree[i] = ij
        tree[j] = ij
        
        # on met à jour les distances des autres noeuds vers ce père
        l = {}
        for u in t:
            if u == i:
                d = (t[i][j] + D[i] - D[j]) / 2
                t[i][ij] = d
                l[i] = d
            elif u == j:
                d = (t[j][i] + D[j] - D[i]) / 2
                t[j][ij] = d
                l[j] = d
            else :
                d = (t[i][u] + t[j][u] - t[i][j]) / 2
                l[u] = d
                t[u][ij] = d 
        t[ij] = l
        
        # on a initialisé les distances du père vers les autres noeuds, ij devient un noeud à traiter (True)
        delta[ij] = [True, {}]
        
        # maintenant on fait disparaître i et j (en mettant False)
        delta[i][0] = False
        delta[j][0] = False
        
    # A la fin des n - 2 itérations, il reste deux noeuds encore à traiter :
    # on relie maintenant les deux derniers noeuds entre eux
    reste = []
    for node in delta:
        if delta[node][0]:
            reste.append(node)
    tree[reste[0]] = reste[1]
    return tree

def print_tree(tree):
    ''' Pour un affichage rapide par Networkx. Si l'utilisateur possède pygraphviz, il est conseillé de plutôt
    utiliser le code de la fonction print_gvtree'''
    G = nx.Graph()
    for node in  tree:
        G.add_node(node)
    for node in tree:
        if node != tree[node]:
#            G.add_edge(tree[node], node, weight = t[node][tree[node]])
            G.add_edge(tree[node], node)
    nx.draw_networkx(G, prog='dot')
    
# Les fonctions suivantes ont pour but d'enraciner l'arbre calculé par Neighbor-Joining
    
def edges_tree(tree):
    ''' A partir de la table tree représentant l'arbre renvoyé par Neighbor-Joining, calcule l'ensemble des
    arêtes du graphe'''
    edges = {}
    for i in tree:
        father = tree[i]
        if father != i:
            edges[i +"-"+ father] = (i, father)
    return edges

def outgoing(edges, tree):
    ''' A partir de l'ensemble des arêtes calculé par la fonction edges_tree ci-dessus, pour chaque noeud
    calcule la liste des arêtes ayant ce noeud pour extrémité dans l'arbre vu comme un graphe non-orienté'''
    out = {}
    for i in tree:
        out[i] = []
    # Pour chaque arête on ajoute chaque extrémité à la liste des voisins de l'autre extrémité
    for e in edges:
        out[edges[e][0]].append(e)
        out[edges[e][1]].append(e)
    return out

def find_leaves(tree):
    ''' A partir de l'arbre construit par Neighbor-Joining, calcule l'ensemble des feuilles, c'est-à-dire les
    taxons de départ'''
    sons = {}
    for i in tree:
        sons[i] = 0
    # Pour chaque noeud on ajoute +1 à son père
    for i in tree:
        father = tree[i]
        sons[father] += 1
    leaves = []
    # Les feuilles sont celles n'ayant pas de fils donc ceux dont sons[i] == 0
    for i in sons:
        if sons[i] == 0:
            leaves.append(i)
    return leaves

def compute_weight(e, node, edges, out, leaves, weight, res):
    ''' Attention : on utilise abusivement le terme de "poids moyen" par la suite. Il correspond en fait à la
    distance moyenne d'un noeud node à toutes les feuilles situées de son côté de l'arbre par rapport à
    l'arête e dont il est une extrémité.
    
    Effectue le calcul récursivement la distance moyenne du noeud node aux feuilles de l'arbre du côté  
    opposé à l'arête e dans l'arbre.
    Le tableau res permet d'effectuer la mémoïzation : res[edge] = [a, b] signifie que si l'on considère 
    edges[edge] = (extremite_gauche, extremite_droite), alors le "poids moyen" des arêtes du côté de 
    l'extrémité gauche de l'arête edge vaut a, et que le poids moyen du côté de l'extrémité droite vaut b.
    
    Arguments : -e est une arête "a-b" du graphe entre deux sommets a et b
                -node est l'une des deux extrémité de l'arête (soit a soit b)
                -edges, out, leaves ont été calculés par les fonctions précédentes et fournissent les 
                informations nécessaires sur la structure du graphe
                -weight correspond au tableau des distances calculé au cours de Neighbor Joining, il sert
                à initialiser les poids des arêtes
                -res est le tableau de mémoïzation pour stocker les résultats
            
    Return : le poids moyen des arêtes du côté de l'extrémité node de l'arête e'''
    
    # Si le résultat demandé est déjà dans le tableau de mémoïzation on le renvoie
    if edges[e][0] == node:
        if res[e][0] != None:
            return res[e][0]
    elif res[e][1] != None:
        return res[e][1]
    
    # Si le noeud node est une feuille, alors on sait que sa distance moyenne aux feuilles est nulle.
    if node in leaves:
        # on note le résultat dans le tableau res avant de renvoyer 0
        if node == edges[e][0]:
            res[e][0] = 0 
        else:
            res[e][1] = 0
        return 0
    
    # Dans le cas général on regarde récursivement les autres arêtes liées à l'extrémité node
    score = 0
    count = 0
    # Pour chaque arête liée au noeud node...
    for edge in out[node]:
        # ...qui n'est pas l'arête e que l'on considère...
        if edge != e:
            # ...on détermine l'extrémité opposée à node...
            if edges[edge][0] == node:
                n = edges[edge][1]
            else:
                n = edges[edge][0]
            count += 1
            # ... et on appelle récursivement le poids de l'arbre de l'autre côté de l'arête edge,
            # auquel on ajoute le poids de l'arête edge.
            w = compute_weight(edge, n, edges, out, leaves, weight, res) + weight[(edges[edge])[0]][(edges[edge])[1]]
            # Enfin on veut sommer pour toutes les arêtes liées au noeud n opposé à node.
            score += w
    # Enfin on divise par le nombre de branches pour obtenir le poids moyen.
    score = score / count
    # On enregistre le résultat avant de le renvoyer
    if node == edges[e][0]:
        res[e][0] = score
    else:
        res[e][1] = score
    return score
    
def find_root(tree, t):
    ''' On veut placer la racine de manière à avoir un arbre équilibré :
        pour chaque arête edge on calcule donc la distance moyenne aux feuilles de part et d'autre de l'arête,
        et si la différence entre les moyennes de chaque côté est inférieure au poids de l'arête edge,
        alors on peut placer la racine sur cette arête.
        
        Arguments : l'arbre tree et la table des distances t mise à jour par Neighbor Joining.
                    t joue le rôle de weight dans la fonction compute_weight.
        
        Return : une arête où l'on va placer la racine et la table res où l'on a stocké les résulats de la fonction
        compute_weight décrite ci-dessus.'''
        
    edges = edges_tree(tree)
    out = outgoing(edges, tree)
    leaves = find_leaves(tree)
    # Initialisation du tableau de mémoïzation
    res = {}
    for edge in edges:
        res[edge] = [None, None]
    for edge in edges:
        compute_weight(edge, edges[edge][0], edges, out, leaves, t, res)
        compute_weight(edge, edges[edge][1], edges, out, leaves, t, res)
    # Recherche d'une arête où l'on peut placer la racine
    for e in res:
        d = abs(res[e][0] - res[e][1])
        if d < t[edges[e][0]][edges[e][1]]:
            return e, res
        
def place_root(tree, edge, edges, weight, res):
    ''' Place logiquement la racine dans l'arbre sur l'arête edge calculée par l'algorithme précédent.
    Cela modifie donc en place la structure de l'arbre tree et le poids des arêtes.'''
    
    node1 = edges[edge][0]
    node2 = edges[edge][1]
    w1 = res[edge][0]
    w2 = res[edge][1]
    
    # On calcule la distance moyenne de la racine aux feuilles
    mid = (w1 + w2 + weight[node1][node2]) / 2.
    # On ajuste les poids des deux arêtes reliant la racine à ses fils
    weight[node1]["root"] = abs(mid - w1)
    weight[node2]["root"] = abs(mid - w2)
    weight["root"] = {}
    weight["root"][node1] = abs(mid - w1)
    weight["root"][node2] = abs(mid - w2)
    
    # On insère la racine dans l'arbre. On inverse ensuite les arêtes qui ne pointent pas dans le bon sens :
    # on veut tree[node] soit égal au père de node dans l'arbre enraciné.
    tree[node1] = "root"
    temp1 = tree[node2]
    tree[node2] = "root"
    temp2 = node2    
    while temp1 != temp2:
        temp = tree[temp1]
        tree[temp1] = temp2
        temp2 = temp1
        temp1 = temp
    
def print_gvtree (tree,t, names) :
    ''' Nécessite pygraphviz et Graphviz.
    Trace l'arbre tree, avec les feuilles portant les noms donnés par la liste names, et 
    avec les poids de la table t qui provient du calcul par l'algorithme de construction de l'arbre.'''
    
    # Crée une table de référence des noeuds et de leur nom pour construire le graphe
    n = {}
    n["root"]="root"
    for k in range(len(names)):
        n[str(k)] =  names[k] 
        
    # Initialisation du graphe
    g= gv.Graph ("Graph")
    g.attr (splines="lines", rankdir='LR')
    
    # Création des noeuds
    c = 0
    for node in tree :
        if node in n:
            name = n[node]
            g.node(name)
        else :
            name = str(c)
            c += 1
            n[node] = name
            g.node (name, shape="point")   # shape = point permet de ne pas le faire apparaître

    # Création des arêtes              
    for node in tree :
        # Si le noeud n'est pas la racine
        if node!= tree[node]:
            g.edge (n[tree[node]], n[node], label =str(int (100*(t[node][tree[node]])) ))
            
    # Affiche le graphe
    g.render( view = True)
    
def tree_upgma (t) :
    ''' Construction de l'arbre par l'algorithme UPGMA à partir de la table t des distances entre les 
    séquences prises deux à deux.
    Renvoie un arbre enraciné tree où tree(i) = j indique que i pointe vers j.'''
    
    tbool = {}
    n= len (t)
    tree = {}
    # Initialise l'arbre
    for i in t :
        tree[i] = i
        
    # Initialise un tableau qui indique s'il faut encore traiter le noeud i (si tbool[i][1]), 
    # tout en gardant en mémoire la table des distances
    # Un autre avantage de cette structure est de réutiliser le code de mini utilisé par Neighbor Joining,
    # qui permet de calculer un minimum sur une table de poids avec des "flag" indiquant si les noeuds sont
    # encore à traiter.
    for i in t : 
        tbool [i] = [True, t[i]]
        
    # On effectue n-2 fusions :
    for k in range (1,n-1):
        # Récupéraion des indices des deux noeuds réalisant la distance minimale
        i, j = mini (tbool)
        # Création du père
        ij = i +" "+ j
        tree[ij] = ij
        tree[i] = ij
        tree[j] = ij
        # Calcul des distances entre le père et les autres noeuds 
        l = {}
        for u in tbool : 
            x = (tbool [u][1][i] + tbool[u][1][j]) /2 
            tbool[u][1][ij] = x
            t[u][ij] = x
            l[u] = x
        l[ij]=0
        t[ij] = l
        # Le père devient un noeud du graphe à traiter
        tbool [ij] = [True, l]
        # Les deux fils ne sont plus à traiter
        tbool [i][0] = False
        tbool [j][0] = False        
    # A la fin on obtient deux arbres qu'il nous faut relier tout en enracinant l'arbre
    i, j = mini (tbool)
    tree[i] = "root"
    tree[j] = "root"
    l = {}
    # On ajuste les distances de la racine aux différents noeuds
    for u in tbool : 
        x = (tbool [u][1][i] + tbool[u][1][j]) /2 
        tbool [u][1]["root"] =x
        t [u]["root"]=x
        l [u] = x
    t["root"] = l
    tbool ["root"] = [True, l]
    tbool [i][0] = False
    tbool [j][0] = False     
    tree ["root"]=  "root"
    return tree

def reverse_tree(tree):
    ''' Inverse un arbre codé suivant tree[fils] = père en tree[père] = [fils1, fils2]'''
    rev = {}
    for node in tree:
        parent = tree[node]
        if not parent in rev:
            rev[parent] = []
        rev[parent].append(node)
    return rev

###################################################### I/ 2 Sequence alignment with linear Needleman-Wunsch ##################################
    
def NW_affine(seq1, seq2, g, e, cout): 
    """ Linear Needleman-Wunsch """
    # cost is a function, g (gap opening cost) and e (gap pursuit cost) must be positive
    """ returns (score,alignment)"""
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
            m_ = M[i-1][j-1] + cout(seq1[i-1], seq2[j-1])
            ix = Ix[i-1][j-1] + cout(seq1[i-1], seq2[j-1])
            iy = Iy[i-1][j-1] + cout(seq1[i-1], seq2[j-1])
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

def affiche(seq1, seq2, traceback):
    """ print aligned sequences """
    L = len(traceback)
    l1 = ""
    l2 = ""
    i, j = 0, 0
    for k in range(L):
        if traceback[k] == "1":
            l1 = l1 + seq1[i]
            i += 1
            l2 = l2 + "-"
        elif traceback[k] == "2":
            l1 = l1 + "-"
            l2 = l2 + seq2[j]
            j += 1
        else:
            l1 = l1 + seq1[i]
            i += 1
            l2 = l2 + seq2[j]
            j += 1
    print(l1)
    print(l2)

###################################################################II/ N sequence alignment with complex linear Needleman-Wunsch ############################
def NW_affine_multi(seq1, seq2, g, e, cout): # coût est une fonction, g et e doit être positif
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
                    m_ += cout(seq1[index1][i-1], seq2[index2][j-1])
                    ix += cout(seq1[index1][i-1], seq2[index2][j-1])
                    iy += cout(seq1[index1][i-1], seq2[index2][j-1])
            
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

def affiche_multi(seq1, seq2, traceback):
    L = len(traceback)
    nb1 = len(seq1)
    nb2 = len(seq2)
    l1 = [""]*nb1
    l2 = [""]*nb2
    i, j = 0, 0
    for k in range(L):
        if traceback[k] == "1":
            for index1 in range(nb1):
                l1[index1] = l1[index1] + seq1[index1][i]
            i += 1
            for index2 in range(nb2):
                l2[index2] = l2[index2] + "*"
        elif traceback[k] == "2":
            for index1 in range(nb1):
                l1[index1] = l1[index1] + "*"
            for index2 in range(nb2):
                l2[index2] = l2[index2] + seq2[index2][j]
            j += 1
        else:
            for index1 in range(nb1):
                l1[index1] = l1[index1] + seq1[index1][i]
            i += 1
            for index2 in range(nb2):
                l2[index2] = l2[index2] + seq2[index2][j]
            j += 1
    return l1,l2

################################################################### III/ Scoring of alignment quality ###################################################

def bali_score(reffile,testfile):

    # Read reference FASTA file
    refkeylist = []
    refseq = {}
    fh = open(reffile)
    # loop over lines
    for line in fh:
      line = line.rstrip()
      # read header
      if line[0] == '>':
        key = line[1:]
        refkeylist.append(key)
        refseq[key] = ''
      # read sequence
      else:
        refseq[key] += line 
    fh.close()
    
    # Read test FASTA file
    testkeylist = []
    testseq = {}
    fh = open(testfile)
    # loop over lines
    for line in fh:
      line = line.rstrip()
      # read header
      if line[0] == '>':
        key = line[1:]
        testkeylist.append(key)
        testseq[key] = ''
      # read sequence
      else:
        testseq[key] += line 
    fh.close()
    
    # Number of sequences
    nseqs = len(testkeylist)
    print('nseqs = %d' % nseqs)
    # verify that the number of sequences is identical between test and ref
    if len(testkeylist) != len(refkeylist):
      raise Exception('Different number of sequences in test (%d) and ref (%d) FASTA' % (len(testkeylist),len(refkeylist)))
    # verify that the keys are identical between test and ref
    if sorted(testkeylist) != sorted(refkeylist):
      raise Exception('Different keys in test and ref FASTA')
    
    # Number of columns in test alignment
    ntestcols = len(testseq[testkeylist[0]])
    print('ntestcols = %d' % ntestcols)
    # verify that all sequences have the same number of columns
    for key in testkeylist:
      if len(testseq[key]) != ntestcols:
        raise Exception('Different number of columns %d != %d for key %s in test FASTA' % (len(testseq[key]),ntestcols,key))
    
    # Number of columns in reference alignment
    nrefcols = len(refseq[refkeylist[0]])
    print('nrefcols = %d' % nrefcols)
    # verify that all sequences have the same number of columns
    for key in refkeylist:
      if len(refseq[key]) != nrefcols:
        raise Exception('Different number of columns %d != %d for key %s in ref FASTA' % (len(refseq[key]),nrefcols,key))
    
    # Calculate BALI scores
    
    # Find columns with gaps in the reference sequence
    # set refseqcol[i] = 0 if ngaps >= cutoff, nseqs - ngaps otherwise
    # here we use 20% of the number of sequences for cutoff
    cutoff = int(nseqs*20.0/100.0)
    if cutoff < 1: cutoff = 1
    refseqcol = [None] * nrefcols
    for i in range(nrefcols):
      ngaps = 0
      for key in refkeylist:
        if refseq[key][i] == '-': ngaps += 1
      if ngaps >= cutoff: refseqcol[i] = 0
      else: refseqcol[i] = nseqs - ngaps
    
    # Code the reference alignment
    # assign to each residue the number of the column it is in
    # gap positions are coded 0
    refcode = {}
    for key in refkeylist:
      refcode[key] = [None] * nrefcols
      for i in range(nrefcols):
        if refseq[key][i] == '-': refcode[key][i] = 0
        else: refcode[key][i] = i+1
    
    # Calculate the maximum score possible
    # ie the score for the reference alignment
    maxsp = 0
    for i in range(nrefcols):
      if refseqcol[i] > 1:
        maxsp += refseqcol[i]*(refseqcol[i]-1)/2.0
    if maxsp <= 0: raise Exception('Error in reference alignment')
    
    # Code the test alignment
    # look up each residue from the test alignment in the reference alignment
    # and assign the reference column number
    testcode = {}
    for key in testkeylist:
      testcode[key] = [0] * ntestcols
      # find the first residue in the reference sequence
      ix = 0
      for j in range(nrefcols):
        if refseq[key][j] != '-':
          ix = refcode[key][j]
          break
      for i in range(ntestcols):
        if testseq[key][i] == '-':
          testcode[key][i] = 0
        else:
          if refseqcol[ix-1] > 0 and refseq[key][ix-1] != '-':
            testcode[key][i] = ix
          for j in range(j+1,nrefcols):
            if refseq[key][j] != '-':
              ix = refcode[key][j]
              break
    
    # Calculate column scores
    ncols = 0
    tc = 0
    sp = 0
    for i in range(ntestcols):
      colscore1 = 0
      index = [None] * nseqs
      scores = [None] * nseqs
      for s in range(nseqs):
        scores[s] = 0
      n = 0
      for key in testkeylist:
        if testcode[key][i] != 0:
          found = False
          for s in range(n):
            if testcode[key][i] == index[s]:
              scores[s] += 1
              found = True
              break
          if found == False:
            scores[n] = 1
            index[n] = testcode[key][i]
            n += 1
      for s in range(nseqs):
        if scores[s] > 1: sp += scores[s]*(scores[s]-1)/2.0
      if testcode[testkeylist[0]][i] > 0 and scores[0] >= refseqcol[testcode[testkeylist[0]][i]-1]: colscore1 = 1
      if testcode[testkeylist[0]][i] != 0: ncols += 1
      tc += colscore1
    if ncols > 0: tc = int(100*float(tc)/float(ncols))/100.0
    sp /= maxsp
    
    # Print results
    print('SP = %.3f' % sp)
    print('TC = %.3f' % tc)
    
    return sp, tc
    
######################################### IV/ N sequence alignment with complex linear Needleman-Wunsch using structural information ###################

from Bio.PDB import HSExposureCB

def get_descriptors(file):
    ''' Gets descriptors for each residue in the protein's structure, 
        including the HSExposure and the coordinates of the alpha carbon '''
    parser = MMCIFParser()
    structure = parser.get_structure(file.split('.')[0],file)
    vca = []
    HSE = []
    names = []
    model = structure[0]
    hse = HSExposureCB(model)
    for chain in model:
        for residue in chain:
            names.append(residue.get_resname())
            if residue.has_id('CA'):
                vca.append(residue['CA'].get_vector())
                hse_ = hse[(chain.id, residue.id)]
                HSE.append((hse_[0], hse_[1]))
            else:
                vca.append(None)
                HSE.append(None)
    dic = {"name" : names, "hse" : HSE, "coord" : vca}
    return dic
                     