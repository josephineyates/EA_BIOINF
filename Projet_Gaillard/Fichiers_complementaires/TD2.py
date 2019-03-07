# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:04:55 2019

@author: pepou
"""

def cout1(a, b): # naif
    if a == b:
        return 1
    else:
        return -1
    
def extract_file_bla(filename):
    '''extrait la matrice blosum d'un fichier type 'blosum35.bla' avec la séquence des aa correspondant'''
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
        line = file.readline()
        if line[0] == " ":
            line = line[1:]
        l = line.split("  ")
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
    
BLOSUM = {}
def blosum():
    global BLOSUM
    if BLOSUM == {}:
        BLOSUM = extract_file_bla('blosum62.bla')
    
def cout_blosum(a, b, mat=None): #mat doit être la matrice BLOSUM
    if mat == None:
        blosum()
        mat = BLOSUM
    return mat[a][b]

def Needleman_Wunsch(seq1, seq2, pen, cout): # coût est une fonction, pen doit être positif
    n = len(seq1)
    m = len(seq2)
    score = [[0 for i in range(m+1)] for j in range(n+1)]
    traceback = [[0 for i in range(m+1)] for j in range(n+1)]
    
    for i in range(1, n+1):
        score[i][0] = -i*pen
        traceback[i][0] = 1 #up
    for j in range(1, m+1):
        score[0][j] = -j*pen
        traceback[0][j] = 2 #left
        
    for i in range(1, n+1):
        for j in range(1, m+1):
            left = score[i][j-1] - pen
            up = score[i-1][j] - pen
            diag = score[i-1][j-1] + cout(seq1[i-1], seq2[j-1])
            if left > up:
                if left > diag:
                    score[i][j] = left
                    traceback[i][j] = 2
                else:
                    score[i][j] = diag
                    traceback[i][j] = 3 #diag
            else:
                if up > diag:
                    score[i][j] = up
                    traceback[i][j] = 1
                else:
                    score[i][j] = diag
                    traceback[i][j] = 3 #diag
    
    # Now the traceback :
    i = n
    j = m
    out = ""
    while (i + j > 0):
        if traceback[i][j] == 1:
            out = "1" + out
            i -= 1
        elif traceback[i][j] == 2:
            out = "2" + out
            j -= 1
        else:
            out = "-" + out
            i -= 1
            j -= 1
    return score[n][m], out

def NW_affine(seq1, seq2, g, e, cout): # coût est une fonction, g et e doit être positif
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
            