# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 14:41:53 2019

@author: josep
"""

############################################ Open PDB files and use them ######################################
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import numpy as np
######################################## example mmcif #########################################################
"""mmcif_dict = MMCIF2Dict('fa/1BF5_A.cif')
mmcif_info = mmcif_dict['_atom_site.Cartn_y']"""

def saveFASTA(doc):
    records = []
    for seq_record in SeqIO.parse(doc, "fasta"):
        records.append(seq_record)
    return records

def download_PDB_file():
    fichiers = os.listdir("balibase/RV11.unaligned")
    for file in fichiers:
        records = saveFASTA("balibase/RV11.unaligned/"+file)
        ids=[]
        for record in records:
            ids.append(record.id.split("_")[0])
        for i in ids:
            pdbl = PDBList()
            pdbl.retrieve_pdb_file(i,pdir="PDB")
        
#download_PDB_file()
            
def change_PDB_files():
    noms_dir = os.listdir("PDB")
    for n in noms_dir:
        fichiers = os.listdir(n)
        for fich in fichiers:
            os.rename("C:/Users/josep/Documents/3A/BIOINF588/TD4/PDB/"+n+"/"+fich,"C:/Users/josep/Documents/3A/BIOINF588/TD4/PDB/"+fich)

#change_PDB_files()
        
def get_info_mmcif(file):
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


def check_corr(fastafile):
    rec = saveFASTA("balibase/RV11.unaligned/"+fastafile+".fasta")
    for r in rec:
        print(get_sequence_mmcif("PDB/"+r+".cif"))
"""
_,_,enf=get_info_mmcif("PDB/1ajs.cif")
for chain in enf.keys():
    print(len(enf[chain]))"""

################### Enregistrement enfouissement #################################
def get_enf():
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
        
            
