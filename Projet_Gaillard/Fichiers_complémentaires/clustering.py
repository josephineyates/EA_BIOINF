# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:28:40 2018

@author: pepou
"""

import random
import time
#import distance
import networkx as nx
import graphviz as gv
import matplotlib.pyplot as plt

#def distancerapide_matrix(seqs, prot):
#    ''' Calcule la matrice des distances des séquences prises deux à deux, renvoie également la liste des 
#    noms correspondants 
#    Arguments : la liste des séquences à traiter et prot est soit "Hemo", soit "Myo", qui sert pour calculer
#    la distance entre les paires de séquences (l'algorithme est légèrement différent en cas de "Hémo".'''
#    n = len(seqs)
#    names = [x[0] for x in seqs]
#    t = [[0] * n for k in range(n)]
#    
#    # Calcul des distances des paires de séquences
#    for i in range(n):
#        for j in range(i+1, n):
#            d = distance.distancerapide(seqs[i][1], seqs[j][1], prot)
#            t[i][j] = d
#            t[j][i] = d
#    for x in t :
#        print (x)
#    return names, t
#
#def distanceBlosum_matrix(seqs, prot, name):
#    ''' Calcule la matrice des distances des séquences prises deux à deux, renvoie également la liste des 
#    noms correspondants 
#    Arguments : la liste des séquences à traiter et prot est soit "Hemo", soit "Myo", qui sert pour calculer
#    la distance entre les paires de séquences (l'algorithme est légèrement différent en cas de "Hémo".
#                Name est le nom de la matrice Blosum à utiliser : 'blosum62.bla' par exemple. '''
#    n = len(seqs)
#    aa, mat = distance.extract_file_bla(name)
#    names = [x[0] for x in seqs]
#    t = [[0] * n for k in range(n)]
#    
#    # Calcul des distances des paires de séquences
#    for i in range(n):
#        for j in range(i+1, n):
#            d = distance.distanceBlosum(seqs[i][1], seqs[j][1], prot, aa, mat)
#            t[i][j] = d
#            t[j][i] = d
#    for x in t :
#        print (x)
#    return names, t

def to_table(t):
    ''' Pour le rendre utilisable par une première version, on transforme la matrice des distances 
    en une table à deux clés qui sont des String d'entiers (plus facilement manipulables par
    l'algorithme Neighbor Joining).'''
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

def print_tree(tree, t):
    ''' Pour un affichage rapide par Networkx. Si l'utilisateur possède pygraphviz, il est conseillé de plutôt
    utiliser le code de la fonction main'''
    G = nx.Graph()
    for node in  tree:
        G.add_node(node)
    for node in tree:
        if node != tree[node]:
#            G.add_edge(tree[node], node, weight = t[node][tree[node]])
            G.add_edge(tree[node], node)
    nx.draw_networkx(G, prog='dot')
    plt.show()
    
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
    Trace l'arbre tree, avec les taxons (ou feuilles) portant les noms donnés par la liste names, et 
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
        # S'il fait partie de la liste des taxons onj crée un noeud portant le nom de l'espèce correspondante
        if node in n:
            name = n[node]
            g.node(name)
        # Sinon on lui attribue un nom (0 puis 1 puis 2 ...) mais il n'apparaît pas comme un noeud,
        # seulement comme une intersection d'arêtes (les noeuds internes ne nous intéressent pas dans un
        # arbre phylogénétique).
        else :
            name = str(c)
            c += 1
            n[node] = name
            g.node (name, shape="point")   # shape = point permet de ne pas le faire apparaître

    # Création des arêtes              
    for node in tree :
        # Si le noeud n'est pas la racine
        if node!= tree[node]:
            # On choisit une couleur au hasard et on colorie l'arête
            color =colorrandom()
            g.edge (n[tree[node]], n[node], color=color, label =str(int (100*(t[node][tree[node]])) ))
            
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


def main():
    ''' Ici il est possible de modifier le code pour changer des paramètres :
        - choisir la base de séquences protéiques (hémoglobine [ou Mammifères qui est plus complète] 
        ou myoglobine)
        - choisir l'algorithme de construction de l'arbre (Neighbor Joining ou UPGMA) '''
    
    start = time.clock()
    
    # Choisir l'un des trois fichiers pour en extraire les séquences des protéines
    seqs= distance.extract_prot('Mammifères.txt')
#    seqs = distance.extract_prot('Hemoglobin_database.txt')
#    seqs = distance.extract_prot('Myoglobin_database.txt')
    print("\nTIME Data extraction : " + str(time.clock() - start))
    
    # Construit la matrice des distances
    # A ajuster selon que l'on ait choisit de l'hémoglobine ou la myoglobine (Mammifères : hémoglobine)
    # ou que l'on veuille la distance construite à partir des scores BLOSUM
    names, t = distancerapide_matrix(seqs,  "Hemo") 
#    names, t = distancerapide_matrix(seqs,  "Myo") 
#    names, t = distanceBlosum_matrix(seqs,  "Hemo", 'blosum62.bla') 
#    names, t = distanceBlosum_matrix(seqs,  "Myo", 'blosum62.bla')
    
#    # Décommenter cette section pour avoir l'arbre combinant les distances calculées sur les deux
#    # types de protéines
#    seqs = distance.extract_prot('Myoglobin_database.txt')
#    names, t = distancerapide_matrix(seqs,  "Myo")
#    seqs = distance.extract_prot('Hemoglobin_database.txt')
#    t1=distancerapide_matrix(seqs,  "Hemo")[1]
#    for i in range (len(t)):
#        for j in range (len (t)):
#            t[i][j]=t[i][j]+t1[i][j]
    
    # crée la table utilisée par les algorithmes de construction des arbres
    table = to_table(t)

    # Construction de l'arbre en utilisant Neighbor-Joining
    start = time.clock()
    tree =neighbor_joining (table)
    # Enracinement de l'arbre
    edges = edges_tree(tree)
    e, res = find_root(tree, table)
    place_root(tree, e, edges, table, res)
    print("\nTIME Neighbor Joining : " + str(time.clock() - start))
    # Affichage de l'arbre
    print_gvtree( tree, table, names )
    
#    # Construction de l'arbre en utilisant UPGMA
#    start = time.clock()
#    tree1=tree_upgma(table)
#    print("\nTIME UPGMA : " + str(time.clock() - start))
#    print_gvtree( tree1, table, names )
    
def colorrandom ():
    ''' Retourne une couleur choisie au hasard. Cette fonction sert à colorier les arêtes du graphe'''
    return (COLORS[random.randint(0,388)]) 
    
COLORS = ['snow', 'ghost white', 'white smoke', 'gainsboro', 'floral white', 'old lace',
    'linen', 'antique white', 'papaya whip', 'blanched almond', 'bisque', 'peach puff',
    'navajo white', 'lemon chiffon', 'mint cream', 'azure', 'alice blue', 'lavender',
    'lavender blush', 'misty rose', 'dark slate gray', 'dim gray', 'slate gray',
    'light slate gray', 'gray', 'light grey', 'midnight blue', 'navy', 'cornflower blue', 'dark slate blue',
    'slate blue', 'medium slate blue', 'light slate blue', 'medium blue', 'royal blue',  'blue',
    'dodger blue', 'deep sky blue', 'sky blue', 'light sky blue', 'steel blue', 'light steel blue',
    'light blue', 'powder blue', 'pale turquoise', 'dark turquoise', 'medium turquoise', 'turquoise',
    'cyan', 'light cyan', 'cadet blue', 'medium aquamarine', 'aquamarine', 'dark green', 'dark olive green',
    'dark sea green', 'sea green', 'medium sea green', 'light sea green', 'pale green', 'spring green',
    'lawn green', 'medium spring green', 'green yellow', 'lime green', 'yellow green',
    'forest green', 'olive drab', 'dark khaki', 'khaki', 'pale goldenrod', 'light goldenrod yellow',
    'light yellow', 'yellow', 'gold', 'light goldenrod', 'goldenrod', 'dark goldenrod', 'rosy brown',
    'indian red', 'saddle brown', 'sandy brown',
    'dark salmon', 'salmon', 'light salmon', 'orange', 'dark orange',
    'coral', 'light coral', 'tomato', 'orange red', 'red', 'hot pink', 'deep pink', 'pink', 'light pink',
    'pale violet red', 'maroon', 'medium violet red', 'violet red',
    'medium orchid', 'dark orchid', 'dark violet', 'blue violet', 'purple', 'medium purple',
    'thistle', 'snow2', 'snow3',
    'snow4', 'seashell2', 'seashell3', 'seashell4', 'AntiqueWhite1', 'AntiqueWhite2',
    'AntiqueWhite3', 'AntiqueWhite4', 'bisque2', 'bisque3', 'bisque4', 'PeachPuff2',
    'PeachPuff3', 'PeachPuff4', 'NavajoWhite2', 'NavajoWhite3', 'NavajoWhite4',
    'LemonChiffon2', 'LemonChiffon3', 'LemonChiffon4', 'cornsilk2', 'cornsilk3',
    'cornsilk4', 'ivory2', 'ivory3', 'ivory4', 'honeydew2', 'honeydew3', 'honeydew4',
    'LavenderBlush2', 'LavenderBlush3', 'LavenderBlush4', 'MistyRose2', 'MistyRose3',
    'MistyRose4', 'azure2', 'azure3', 'azure4', 'SlateBlue1', 'SlateBlue2', 'SlateBlue3',
    'SlateBlue4', 'RoyalBlue1', 'RoyalBlue2', 'RoyalBlue3', 'RoyalBlue4', 'blue2', 'blue4',
    'DodgerBlue2', 'DodgerBlue3', 'DodgerBlue4', 'SteelBlue1', 'SteelBlue2',
    'SteelBlue3', 'SteelBlue4', 'DeepSkyBlue2', 'DeepSkyBlue3', 'DeepSkyBlue4',
    'SkyBlue1', 'SkyBlue2', 'SkyBlue3', 'SkyBlue4', 'LightSkyBlue1', 'LightSkyBlue2',
    'LightSkyBlue3', 'LightSkyBlue4', 'SlateGray1', 'SlateGray2', 'SlateGray3',
    'SlateGray4', 'LightSteelBlue1', 'LightSteelBlue2', 'LightSteelBlue3',
    'LightSteelBlue4', 'LightBlue1', 'LightBlue2', 'LightBlue3', 'LightBlue4',
    'LightCyan2', 'LightCyan3', 'LightCyan4', 'PaleTurquoise1', 'PaleTurquoise2',
    'PaleTurquoise3', 'PaleTurquoise4', 'CadetBlue1', 'CadetBlue2', 'CadetBlue3',
    'CadetBlue4', 'turquoise1', 'turquoise2', 'turquoise3', 'turquoise4', 'cyan2', 'cyan3',
    'cyan4', 'DarkSlateGray1', 'DarkSlateGray2', 'DarkSlateGray3', 'DarkSlateGray4',
    'aquamarine2', 'aquamarine4', 'DarkSeaGreen1', 'DarkSeaGreen2', 'DarkSeaGreen3',
    'DarkSeaGreen4', 'SeaGreen1', 'SeaGreen2', 'SeaGreen3', 'PaleGreen1', 'PaleGreen2',
    'PaleGreen3', 'PaleGreen4', 'SpringGreen2', 'SpringGreen3', 'SpringGreen4',
    'green2', 'green3', 'green4', 'chartreuse2', 'chartreuse3', 'chartreuse4',
    'OliveDrab1', 'OliveDrab2', 'OliveDrab4', 'DarkOliveGreen1', 'DarkOliveGreen2',
    'DarkOliveGreen3', 'DarkOliveGreen4', 'khaki1', 'khaki2', 'khaki3', 'khaki4',
    'LightGoldenrod1', 'LightGoldenrod2', 'LightGoldenrod3', 'LightGoldenrod4',
    'LightYellow2', 'LightYellow3', 'LightYellow4', 'yellow2', 'yellow3', 'yellow4',
    'gold2', 'gold3', 'gold4', 'goldenrod1', 'goldenrod2', 'goldenrod3', 'goldenrod4',
    'DarkGoldenrod1', 'DarkGoldenrod2', 'DarkGoldenrod3', 'DarkGoldenrod4',
    'RosyBrown1', 'RosyBrown2', 'RosyBrown3', 'RosyBrown4', 'IndianRed1', 'IndianRed2',
    'IndianRed3', 'IndianRed4', 'sienna1', 'sienna2', 'sienna3', 'sienna4', 'burlywood1',
    'burlywood2', 'burlywood3', 'burlywood4', 'wheat1', 'wheat2', 'wheat3', 'wheat4', 'tan1',
    'tan2', 'tan4', 'chocolate1', 'chocolate2', 'chocolate3', 'firebrick1', 'firebrick2',
    'firebrick3', 'firebrick4', 'brown1', 'brown2', 'brown3', 'brown4', 'salmon1', 'salmon2',
    'salmon3', 'salmon4', 'LightSalmon2', 'LightSalmon3', 'LightSalmon4', 'orange2',
    'orange3', 'orange4', 'DarkOrange1', 'DarkOrange2', 'DarkOrange3', 'DarkOrange4',
    'coral1', 'coral2', 'coral3', 'coral4', 'tomato2', 'tomato3', 'tomato4', 'OrangeRed2',
    'OrangeRed3', 'OrangeRed4', 'red2', 'red3', 'red4', 'DeepPink2', 'DeepPink3', 'DeepPink4',
    'HotPink1', 'HotPink2', 'HotPink3', 'HotPink4', 'pink1', 'pink2', 'pink3', 'pink4',
    'LightPink1', 'LightPink2', 'LightPink3', 'LightPink4', 'PaleVioletRed1',
    'PaleVioletRed2', 'PaleVioletRed3', 'PaleVioletRed4', 'maroon1', 'maroon2',
    'maroon3', 'maroon4', 'VioletRed1', 'VioletRed2', 'VioletRed3', 'VioletRed4',
    'magenta2', 'magenta3', 'magenta4', 'orchid1', 'orchid2', 'orchid3', 'orchid4', 'plum1',
    'plum2', 'plum3', 'plum4', 'MediumOrchid1', 'MediumOrchid2', 'MediumOrchid3',
    'MediumOrchid4', 'DarkOrchid1', 'DarkOrchid2', 'DarkOrchid3', 'DarkOrchid4',
    'purple1', 'purple2', 'purple3', 'purple4', 'MediumPurple1', 'MediumPurple2',
    'MediumPurple3', 'MediumPurple4', 'thistle1', 'thistle2', 'thistle3', 'thistle4',
    'gray1',  'gray21',  'gray40', 'gray66',  'gray77',  'gray88',  'gray99']  

if __name__ == '__main__':
    main()