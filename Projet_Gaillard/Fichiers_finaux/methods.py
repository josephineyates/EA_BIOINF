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
        E/ Template class for Needleman-Wunsch algorithm
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
BLOSUM = extract_file_bla('blosum62.bla')

############################################################ B/ PDB analysis functions #####################################################"
    
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
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
    for seqrec in SeqIO.parse(handle, "fasta"):
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
            
from Bio.PDB import HSExposureCB
        
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
        
def decode_enf(file_id, chain=None):
    ''' Gets all the values of the residue's depth 
        Returns the list of "enf" values in the sequential order '''
    file = open("enfouissements.txt", "r")
    content = file.read()
    lines = content.split("\n")[:-1]
    N = len(lines)
    good_file = False
    if chain == None:
        chain="A"
    good_chain=False
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
        elif "Chain" in lines[k]:
            id_ = lines[k].split("id=")[1].split(">")[0]
            good_chain = (id_ == chain)          
        elif not "hain" in lines[k] and good_file and good_chain:
            tokens = lines[k].split(",")
            values.append(float(tokens[1]))
    return values

AA3to1 = {"ALA":"A", "ARG":"R", "ASP":"D", "ASN":"N", "CYS":"C", "GLU":"E", "GLN":"Q", "GLY":"G", "HIS":"H", "ILE":"I",
          "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V"}

def get_descriptors(file, seq, chain=None):
    ''' Gets descriptors for each residue in the protein's structure, 
        including the HSExposure and the coordinates of the alpha carbon '''
    print(file)
    parser = MMCIFParser()
    structure = parser.get_structure(file.split('.')[0],file)
    vca = []
    HSE = []
    names = []
    model = structure[0]
    good = False
    if chain == None:
        chain="A"
    i = 1
    while not good and i <= len(structure):
        try:
            hse = HSExposureCB(model)
            good = True
        except:
            if i == len(structure):
                raise Exception("No model worked for HSE")
            print("HSE model error")
            model = structure[i]
            i += 1
            good = False
    for chain_ in model:
        if chain==chain_.id:           
            for residue in chain_:
                if residue.get_resname() in AA3to1:
                    nem = AA3to1[residue.get_resname()]
                else:
                    nem = "X"
                names.append(nem)
                if residue.has_id('CA'):
                    vca.append(residue['CA'].get_vector())
                else:
                    vca.append(None)
                HSE.append(None)
    for residue in hse:
        hse_ = residue[1]
        if residue[0].get_full_id()[2] == chain and residue[0].get_full_id()[3][1] <= len(HSE):
            HSE[residue[0].get_full_id()[3][1] - 1] = (hse_[0], hse_[1])
    s = ""
    for c in names:
        s = s + c
    pos = s.find(str(seq))
    assert pos > -1, "No compatibility between PDB and fasta sequence for {}".format(file)
    N = len(seq)
    dic = {"name" : names[pos:pos+N], "hse" : HSE[pos:pos+N], "coord" : vca[pos:pos+N]}
    print("HSE errors : {} out of {}".format(HSE[pos:pos+N].count(None), N))
    return dic

############################################################## D/ Clustering functions ######################################################"

import networkx as nx
import graphviz as gv

def to_table(t):
    ''' Pour le rendre utilisable on transforme la matrice des distances 
    en une table à deux clés qui sont des String d'entiers (plus facilement manipulables par
    l'algorithme Neighbor Joining).'''
    table = {}
    for i in range(len(t)):
        l = {}
        for j in range(len(t)):
            l[str(j)] = t[i][j]
        table[str(i)] = l
    return table

def mini_t(table):
    ''' Calcule le minimum de la table t dans Neighbor Joining et renvoie sa valeur '''
    mini = None
    for i in table:
        for j in table[i]:
            if mini == None:
                mini = table[i][j]
            else:
                if mini > table[i][j]:
                    mini = table[i][j]
    return mini

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
    
    # On ramène toutes les valeurs en positif avant d'effectuer l'algorithme :
    offset = mini_t(t)
    if offset < 0:
        for i in t:
            for j in t[i]:
                t[i][j] = t[i][j] - 2 * offset
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
        if i != father:
            sons[father] += 1
    leaves = []
    # Les feuilles sont celles n'ayant pas de fils donc ceux dont sons[i] == 0, et dans le cas de deux noeuds seulement, il faut compter sons < 2 
    for i in sons:
        if sons[i] < 2:
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

###################################################### E/ Template class for Needleman-Wunsch algorithm ######################################
    
# Different cost functions that may be used in the NW algo
    
from time import time
from random import shuffle
    
def cout_blosum(seq1, i, seq2, j, mat=BLOSUM): #mat doit être la matrice BLOSUM
    a = seq1["seq"][i]
    b = seq2["seq"][j]
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

# Gap and extend functions taking the enf into account :

def gap_profiled(seq1, i, g=11, p=0.5):
    if "enf" in seq1:
        return g * p + g * (1. - p) * seq1["enf"][i] /seq1["enf_max"]
    else:
        return g
    
def extend_profiled(seq1, i, e=1, p=0.5):
    if "enf" in seq1:
        return e * p + e * (1. - p) * seq1["enf"][i] /seq1["enf_max"]
    else:
        return e 
    
class NeedlemanWunsch:
    ''' Class for Needleman-Wunsch algorithm
        Can be initialized with different parameters :
            - cost_functions : list of cost functions that applies when comparing two residues
            - cost_coefs : list of coefficients to multiply the previous cost functions. Sum must be 1
            - gap_opening_function and gap_extending_function : how to penalize the alignment between a residue and a gap
            - clustering : method for clustering, ie "NJ" for Neighbor Joining or "UPGMA"
    '''
    
    def __init__(self, cost_functions=cout_blosum, gap_opening_function=11, gap_extending_function=1, clustering="NJ", cost_coefs=None):
        if isinstance(gap_opening_function, int) or isinstance(gap_opening_function, float):
            gap_ = lambda x, y : gap_opening_function
            self.gap = gap_
        else:
            self.gap = gap_opening_function
        if isinstance(gap_extending_function, int) or isinstance(gap_extending_function, float):
            extend_ = lambda x, y : gap_extending_function
            self.extend = extend_
        else:
            self.extend = gap_extending_function
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
        self.hidden_costs = cost_functions
        self.costs = [self.wrap(cost_func) for cost_func in cost_functions]
        self.clustering = clustering
        
    def wrap(self, function):
        e = self.extend
        g = self.gap
        def f(seq1, i, seq2, j):
            a = seq1["seq"][i]
            b = seq2["seq"][j]
            if a == "-":
                if b == "-":
                    return 0
                if j > 1 and seq2["seq"][j-1] == "-":
                    return -e(seq2, j)
                return -g(seq2, j)
            if b == "-":
                if i > 1 and seq1["seq"][i-1] == "-":
                    return -e(seq1, i)
                return -g(seq1, i)
            return function(seq1, i, seq2, j)
        return f
        
    def set_coefs(self, new_coefs):
        self.cost_coefs = new_coefs
        
    def set_ge(self, g, e):
        if isinstance(g, int) or isinstance(g, float):
            gap_ = lambda x, y : g
            self.gap = gap_
        else:
            self.gap = g
        if isinstance(e, int) or isinstance(e, float):
            extend_ = lambda x, y : e
            self.extend = extend_
        else:
            self.extend = e
        self.costs = [self.wrap(cost_func) for cost_func in self.hidden_costs]
            
    def set_costs(self, new_costs):
        self.costs = [self.wrap(cost_func) for cost_func in new_costs]
        
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

        M[1][0] = -sum([self.gap(sequences1[j], 0) for j in range(nb1)]) * nb2
        Ix[1][0] = M[1][0]
        Iy[1][0] = M[1][0]
        tM[1][0] = 1 #up
        tIx[1][0] = 1 #up
        tIy[1][0] = 1 #up    
        for i in range(2, n+1):
            M[i][0] = M[i-1][0] - sum([self.extend(sequences1[j],i-1) for j in range(nb1)])*nb2
            Ix[i][0] = M[i][0]
            Iy[i][0] = M[i][0]
            tM[i][0] = 6 #up
            tIx[i][0] = 6 #up
            tIy[i][0] = 6 #up
            
        M[0][1] = -sum([self.gap(sequences2[j], 0) for j in range(nb2)]) * nb1
        Ix[0][1] = M[0][1]
        Iy[0][1] = M[0][1]
        tM[0][1] = 2 #left
        tIx[0][1] = 2 #left
        tIy[0][1] = 2 #left
        for j in range(2, m+1):
            M[0][j] = M[0][j-1] - sum([self.extend(sequences2[i], j-1) for i in range(nb2)])*nb1
            Ix[0][j] = M[0][j]
            Iy[0][j] = M[0][j]
            tM[0][j] = 7 #left
            tIx[0][j] = 7 #left
            tIy[0][j] = 7 #left
        
        for i in range(1, n+1):
            for j in range(1, m+1):
                #M
                m_ = M[i-1][j-1]
                ix = Ix[i-1][j-1]
                iy = Iy[i-1][j-1]
            
                for index1 in range (nb1):
                    for index2 in range (nb2):
                        cout_ab = sum([coef * cout(sequences1[index1], i-1, sequences2[index2], j-1) for (cout, coef) in zip(self.costs, self.cost_coefs)])
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
                m_ = M[i-1][j] - sum([self.gap(sequences2[k], j-1) for k in range(nb2)]) * nb1
                ix = Ix[i-1][j] - sum([self.extend(sequences2[k], j-1) for k in range(nb2)])*nb1
                if m_ > ix:
                    Ix[i][j] = m_
                    tIx[i][j] = 1
                else:
                    Ix[i][j] = ix
                    tIx[i][j] = 6 # up dans Ix -1, 0
                #Iy
                m_ = M[i][j-1] - sum([self.gap(sequences1[k], i-1) for k in range(nb1)]) * nb2
                iy = Iy[i][j-1] - sum([self.extend(sequences1[k], i-1) for k in range(nb1)])*nb2
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
        mod = list(sequences1[0].keys())
        mod.remove("seq")
        mod.remove("name")
        desc1 = {}
        desc2 = {}
        for key in mod:
            desc1[key] = [[]]*nb1
            desc2[key] = [[]]*nb2
        i, j = 0, 0
        for k in range(L):
            if traceback[k] == "1":
                for index1 in range(nb1):
                    l1[index1] = l1[index1] + sequences1[index1]["seq"][i]
                    for key in desc1:
                        desc1[key][index1].append(sequences1[index1][key][i])
                i += 1
                for index2 in range(nb2):
                    l2[index2] = l2[index2] + "-"
                    for key in desc2:
                        desc2[key][index2].append(None)
            elif traceback[k] == "2":
                for index1 in range(nb1):
                    l1[index1] = l1[index1] + "-"
                    for key in desc1:
                        desc1[key][index1].append(None)
                for index2 in range(nb2):
                    l2[index2] = l2[index2] + sequences2[index2]["seq"][j]
                    for key in desc2:
                        desc2[key][index2].append(sequences2[index2][key][j])
                j += 1
            else:
                for index1 in range(nb1):
                    l1[index1] = l1[index1] + sequences1[index1]["seq"][i]
                    for key in desc1:
                        desc1[key][index1].append(sequences1[index1][key][i])
                i += 1
                for index2 in range(nb2):
                    l2[index2] = l2[index2] + sequences2[index2]["seq"][j]
                    for key in desc2:
                        desc2[key][index2].append(sequences2[index2][key][j])
                j += 1
        for index1 in range(nb1):
            sequences1[index1]["seq"] = l1[index1]
            for key in desc2:
                sequences1[index1][key] = desc1[key][index1]
#            print(l1[index1])
        for index2 in range(nb2):
            sequences2[index2]["seq"] = l2[index2]
            for key in desc2:
                sequences2[index2][key] = desc2[key][index2]
#            print(l2[index2])
#        print(" ")
    
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
        cost_matrix = [[0 for i in range(N)] for j in range(N)]
        for i in range(N): 
            for j in range(i):
                cost_matrix[i][j] = self.align([sequences[i]], [sequences[j]])[0]
                cost_matrix[j][i] = cost_matrix[i][j]
        table_dist = to_table (cost_matrix)
        
        #Clustering
        if self.clustering == "NJ":
            tree_NJ = neighbor_joining(table_dist)
            edges = edges_tree(tree_NJ)
            ed,res = find_root(tree_NJ,table_dist)
            place_root(tree_NJ,ed,edges,table_dist,res)
#            print_tree(tree_NJ)
            good_tree = reverse_tree(tree_NJ) # good_tree stocke la liste des fils de chaque noeud s'il n'est pas une feuille
        elif self.clustering == "UPGMA":
            tree_UPGMA = tree_upgma(table_dist)
#            print_tree(tree_UPGMA)
            good_tree = reverse_tree(tree_UPGMA)
        
        g1, g2 = good_tree["root"][0], good_tree["root"][1] 
        return self.tree_align(g1, g2, sequences, good_tree)
    
    def compute_score(self):
        t = time()
        files = os.listdir("balibase/RV11.unaligned/")
        SP, TC = [], []
        for file in files:        
            print(" ")
            rec = readFASTA("balibase/RV11.unaligned/{}".format(file))
            sequences = []
            for r in rec:
                seq = {"seq" : r.seq}
                seq["name"] = r.name
                chain = r.name.split("_")[1]
                if chain=="":
                    chain=None
                if cout_hse in self.costs:
                    desc = get_descriptors("PDB/"+r.name[:4]+".cif", chain=chain)
                    seq["hse"] = desc["hse"]
                if cout_enf in self.costs:
                    values = decode_enf(r.name[:4], chain=chain)
                    seq["enf"] = values
                    seq["enf_max"] = max(values)
                sequences.append(seq)
            sequences = self.run(sequences)
            save_to_fasta(sequences)
            sp, tc = bali_score("balibase/RV11.aligned/{}".format(file), "save_alignment.fasta")
            SP.append(sp)
            TC.append(tc)
        print("TIME : {}".format(time()-t))
        m_tc = np.mean(TC)
        m_sp = np.mean(SP)
        print("TC = {}\nSP = {}".format(m_tc, m_sp))
        return m_sp, m_tc
    
    def quick_score(self):
        t = time()
        files = os.listdir("balibase/RV11.unaligned/")
        SP, TC = [], []
        shuffle(files)
        for file in files:
            if time() - t > 300:
                break
            print(" ")
            rec = readFASTA("balibase/RV11.unaligned/{}".format(file))
            sequences = []
            for r in rec:
                seq = {"seq" : r.seq}
                seq["name"] = r.name
                chain=None
                if "_" in r.name:
                    chain = r.name.split("_")
                    if chain=="":
                        chain=None
                if cout_hse in self.costs:
                    desc = get_descriptors("PDB/"+r.name[:4]+".cif", r.seq, chain=chain)
                    seq["hse"] = desc["hse"]
                if cout_enf in self.costs:
                    values = decode_enf(r.name[:4], chain=chain)
                    seq["enf"] = values
                    seq["enf_max"] = max(values)
                sequences.append(seq)
            sequences = self.run(sequences)
            save_to_fasta(sequences)
            sp, tc = bali_score("balibase/RV11.aligned/{}".format(file), "save_alignment.fasta")
            SP.append(sp)
            TC.append(tc)
        m_tc = np.mean(TC)
        m_sp = np.mean(SP)
        print("TC = {}\nSP = {}".format(m_tc, m_sp))
        return m_sp, m_tc

###################################################### I/ 2 Sequence alignment with linear Needleman-Wunsch ##################################
    
def NW_affine(seq1, seq2, g, e, cout="blosum"): 
    """ Linear Needleman-Wunsch """
    # cost is a function, g (gap opening cost) and e (gap pursuit cost) must be positive
    """ returns (score,alignment)"""
    if cout == "blosum":
        cost_func = lambda sequ1, i, sequ2, j : cout_blosum(sequ1, i, sequ2, j, g=lambda x, y : g, e=lambda x, y : e)
    else:
        cost_func = cout
    NW = NeedlemanWunsch(cost_functions=cost_func, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[1])
    sequence1 = [{"seq":seq1}]
    sequence2 = [{"seq":seq2}]
    sequences = NW.run(sequence1, sequence2)
    print(sequences[0]["seq"])
    print(sequences[1]["seq"])
    return [sequences[0]["seq"], sequences[0]["seq"]]

###################################################################II/ N sequence alignment with complex linear Needleman-Wunsch ############################
    
def NW_affine_multi(seq1, seq2, g, e, cout="blosum"): # coût est une fonction, g et e doit être positif
    
    if cout == "blosum":
        cost_func = lambda sequ1, i, sequ2, j : cout_blosum(sequ1, i, sequ2, j, g=lambda x, y : g, e=lambda x, y : e)
    else:
        cost_func = cout
    NW = NeedlemanWunsch(cost_functions=cost_func, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[1])
    sequence1 = [{"seq":seq} for seq in seq1]
    sequence2 = [{"seq":seq} for seq in seq2]
    sequences = NW.run(sequence1, sequence2)
    for seq in sequences:
        print(seq["seq"])
    return [seq["seq"] for seq in sequences]

################################################################### III/ Scoring of alignment quality ###################################################
from Bio.Align.Applications import ClustalwCommandline   
from Bio import AlignIO
import Bio.SeqIO
import subprocess
import sys

def conversionAlignFASTA(fichier):
	alignment = AlignIO.read(fichier,"clustal")
	file = open(fichier+"_aligned.fasta","w")
	Bio.SeqIO.write(alignment,file,"fasta")
	file.close()

def score_align(clustalfile,reffile):
    #clustalw_exe = r"C:\Users\josep\Anaconda3\Lib\site-packages\Bio\Align\Applications\_Clustalw.py"
    #assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
    #cline=ClustalwCommandline(clustalw_exe, infile=reffile+".fasta", type="PROTEIN", output="FASTA", outfile=reffile+"_aligned.fasta", quiet=True)
    #child = subprocess.call(str(clustalw_exe)+" -align -infile="+reffile_name+".fasta -seqnos ON -output fasta_aln -type protein", shell=True)

    #print(cline())
    """
    #subprocess.check_call(args, *, stdin=None, stdout=None, stderr=None, shell=False)
    #clustalw_cline() 
    align = AlignIO.read(reffile_name+".aln", "clustal")
    print(align)
    conversionAlignFASTA(reffile_name+".aln")
    #Alignseq = align_fasta("BBS11001.fasta")
    bali_score("balibase/RV11.aligned/BBS11001.fasta", reffile_name+".aln_aligned.fasta")"""
        

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
        print(testkeylist, refkeylist)
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

from Bio.Seq import Seq

def align_fasta(file, algo="NW_blosum",silent=True):
    rec = readFASTA("balibase/RV11.unaligned/{}".format(file))
    sequences = []
    g, e = 11, 2
    if algo == "NW_blosum":
        NW = NeedlemanWunsch(cost_functions=cout_blosum, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[1]) 
    elif algo == "NW_structural":
        
        ############# TODO ##################
        NW = NeedlemanWunsch(cost_functions=cout_blosum, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[1]) 
    for r in rec:
        seq = {"seq" : r.seq}
        seq["name"] = r.name
        chain=None
        if "_" in r.name:
            chain = r.name.split("_")
            if chain=="":
                chain=None
        if cout_hse in NW.costs:
            desc = get_descriptors("PDB/"+r.name[:4]+".cif", r.seq, chain=chain)
            seq["hse"] = desc["hse"]
        if cout_enf in NW.costs:
            values = decode_enf(r.name[:4], chain=chain)
            seq["enf"] = values
            seq["enf_max"] = max(values)
        sequences.append(seq)
    sequences = NW.run(sequences)
    if not(silent):
        for seq in sequences:
            print(seq["seq"])
    return sequences

def save_to_fasta(sequences,outputfile):
    seqRec = []
    for seq in sequences:
        seqr = SeqIO.SeqRecord(Seq(seq["seq"]), id=seq["name"], name=seq["name"], description="", dbxrefs=[])
        seqRec.append(seqr)
    with open(outputfile+".fasta", "w") as output_handle:
        SeqIO.write(seqRec, output_handle, "fasta")

def align_score(file):
    seq = align_fasta(file)
    save_to_fasta(seq,file+"_aligned.fasta")
    return bali_score("balibase/RV11.aligned/{}".format(file), file+"_aligned.fasta")
    
##################################### V/ Optimization of cost coefficients ##########################################
    
from bayes_opt import BayesianOptimization

def optimize_unstructural(opt="quick"):
    ''' Optimization of gap and extend parameters in unstructural version '''
    pbounds = {'g': (0, 25), 'e': (0, 1)}
    g, e = 11, 1
    NW = NeedlemanWunsch(cost_functions=cout_blosum, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[1])
    if opt=="all":
        def compute(g, e):
            NW.set_ge(g, g*e)
            sp, tc = NW.compute_score()
            return sp
    elif opt=="quick":
        def compute(g, e):
            NW.set_ge(g, g*e)
            sp, tc = NW.quick_score()
            return sp
    else:
        raise Exception("Not valid option")
    optimizer = BayesianOptimization(
            f=compute,
            pbounds=pbounds,
            random_state=42,
            )
    optimizer.probe(params={"g":11, "e":0.1}, lazy=True)
    optimizer.probe(params={"g":15, "e":6/15}, lazy=True)
    optimizer.maximize(init_points=4, n_iter=30)
    return optimizer
    
def optimize_hse(opt="quick"):
    ''' Optimization of balance between cost functions in the structural case '''
    pbounds = {'x': (0, 1)}
    g, e = 11, 2
    costs = [cout_blosum, cout_hse]
    NW = NeedlemanWunsch(cost_functions=costs, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[0.4, 0.6])
#   NW = NeedlemanWunsch(cost_functions=costs, gap_opening_function=g, gap_extending_function=e, clustering="UPGMA", cost_coefs=[0.4, 0.4, 0.2])
    def compute(x):
        NW.set_coefs([x, (1-x)])
        if opt=="all":
            sp, tc = NW.compute_score()
        if opt=="quick":
            sp, tc = NW.quick_score()
        else:
            raise Exception("Not valid option")
        return sp
    optimizer = BayesianOptimization(
            f=compute,
            pbounds=pbounds,
            random_state=1,
            )
    optimizer.probe(params={"x":0.8}, lazy=True)
    optimizer.maximize(init_points=2, n_iter=10)
    return optimizer

def optimize_enf(opt="quick"):
    ''' Optimization of balance between cost functions in the structural case '''
    pbounds = {'x': (0, 1)}
    g, e = 11, 2
    alpha = 0.5
    costs = [cout_blosum, cout_hse, cout_enf]
    NW = NeedlemanWunsch(cost_functions=costs, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[0.4, 0.4, 0.2])
#   NW = NeedlemanWunsch(cost_functions=costs, gap_opening_function=g, gap_extending_function=e, clustering="UPGMA", cost_coefs=[0.4, 0.4, 0.2])
    def compute(x):
        NW.set_coefs([alpha * x, (1-alpha) * x, (1-x)])
        if opt=="all":
            sp, tc = NW.compute_score()
        if opt=="quick":
            sp, tc = NW.quick_score()
        else:
            raise Exception("Not valid option")
        return sp
    optimizer = BayesianOptimization(
            f=compute,
            pbounds=pbounds,
            random_state=1,
            )
    optimizer.probe(params={"x":0.8}, lazy=True)
    optimizer.maximize(init_points=2, n_iter=10)
    return optimizer

def optimize_gap(opt="quick"):
    ''' Optimization of the enf penalty for gaps '''
    pbounds = {'p' : (0, 1)}
    p = 0.5
    g_opt=11
    e_opt=1
    g = lambda x, y : gap_profiled(x, y, g=g_opt, p=p)
    e = lambda x, y : extend_profiled(x, y, e=e_opt, p=p)
    costs = [cout_blosum]
    NW = NeedlemanWunsch(cost_functions=costs, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[1])
#   NW = NeedlemanWunsch(cost_functions=costs, gap_opening_function=g, gap_extending_function=e, clustering="UPGMA", cost_coefs=[0.4, 0.4, 0.2])
    def compute(p):
        g = lambda x, y : gap_profiled(x, y, g=g_opt, p=p)
        e = lambda x, y : extend_profiled(x, y, e=e_opt, p=p)
        NW.set_ge(g, e)
        if opt=="all":
            sp, tc = NW.compute_score()
        if opt=="quick":
            sp, tc = NW.quick_score()
        else:
            raise Exception("Not valid option")
        return sp
    optimizer = BayesianOptimization(
            f=compute,
            pbounds=pbounds,
            random_state=1,
            )
    optimizer.maximize(init_points=2, n_iter=20)
    return optimizer

def grid_search(opt="quick"):
    g, e = 11, 1
    NW = NeedlemanWunsch(cost_functions=cout_blosum, gap_opening_function=g, gap_extending_function=e, clustering="NJ", cost_coefs=[1])
    if opt=="all":
        def compute(g, e):
            NW.set_ge(g, e)
            sp, tc = NW.compute_score()
            return sp, tc
    elif opt=="quick":
        def compute(g, e):
            NW.set_ge(g, e)
            sp, tc = NW.quick_score()
            return sp, tc
    else:
        raise Exception("Not valid option")
    res={}
    for (g,e) in [(8, 1), (10, 1), (10, 2), (10, 4), (11, 3), (12, 2), (13, 2), (13, 6), (14, 4), (15, 6)]:
        res[(g, e)] = compute(g, e)
    return res

def eval_clustalw():
    sp,tc = {"clustal":[],"align":[]},{"clustal":[],"align":[]}
    files = os.listdir("balibase/RV11.unaligned/")
    for file in files:
        sp_, tc_ = align_score(file)
        sp_clus, tc_clus = bali_score("balibase/RV11.aligned/{}".format(file), "balibase/ClustalW/{}".format(file))
        sp["clustal"].append(sp_clus)
        sp["align"].append(sp_)
        tc["clustal"].append(tc_clus)
        tc["align"].append(tc_)
    mean_sp_clus = np.sum(sp["clustal"])/len(files)
    mean_tc_clus = np.sum(tc["clustal"])/len(files)
    mean_sp = np.sum(sp["align"])/len(files)
    mean_tc = np.sum(tc["align"])/len(files)
    print("\nMean score clustal:")
    print("TC = {}".format(mean_tc_clus))
    print("SP = {}".format(mean_sp_clus))
    
    print("\nMean score alignment:")
    print("TC = {}".format(mean_tc))
    print("SP = {}".format(mean_sp))

eval_clustalw()
#res = {(8, 1): (0.39723254408655229, 0.21473684210526311),
# (10, 1): (0.41608800440831978, 0.25394736842105264),
# (10, 2): (0.50433622944747858, 0.27894736842105256),
# (10, 4): (0.48388987506593545, 0.25868421052631579),
# (11, 3): (0.49784464522354061, 0.26973684210526311),
# (12, 2): (0.50407271589054847, 0.28026315789473683),
# (13, 2): (0.49753416550527296, 0.29026315789473683),
# (13, 6): (0.46949131392905497, 0.255),
# (14, 4): (0.50545745015904853, 0.28078947368421053),
# (15, 6): (0.47925907268134837, 0.26236842105263153)}         
                     