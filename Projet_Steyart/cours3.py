from Bio.ExPASy import Prosite
import requests
import gzip
import shutil
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def donneesProt():
    handle = open("prosite.dat")
    records = Prosite.parse(handle)
    save_file = open("prosite_entries.dat", "w") #Sauvegarde
    for record in records:
        save_file.write(record.accession+',')
        save_file.write(record.name+',')
        save_file.write(record.pattern+',')
        save_file.write(record.pdoc+'\r\n')
    save_file.close()
    records.close()

def prositeGetPloop():
    handle = open("prosite.dat")
    handle = open("prosite.dat")
    records = Prosite.parse(handle)
    for record in records:
        if (record.accession=="PS00017"):
            print(record.pattern)
	
#donneesProt()
#prositeGetPloop()

################################################"

def getFASTA(url):

    filename = url.split("/")[-1]
    with open(filename, "wb") as f:
        r = requests.get(url)
        f.write(r.content)

    with gzip.open(filename, 'rb') as f_in:
        with open('orf_coding_all.fasta', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    statinfo = os.stat('orf_coding_all.fasta')
    print("Size of file :",statinfo.st_size)
    
def saveFASTA(doc):
    records = []
    for seq_record in SeqIO.parse(doc, "fasta"):
        records.append(seq_record)
    return records

def make_RNA_record(records):
    ARN = []
    for record in records:
        ARN.append(SeqRecord(seq = record.seq.transcribe(),id = record.id, description = record.description))
    return ARN

def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq = nuc_record.seq.translate(table="Standard"), id = nuc_record.id, description = nuc_record.description)

def translateFASTA(records):
    #rec1=[]
    #for rec in records:
        #if (len(rec.seq)%3==0):
           #rec1.append(rec)
    proteins = (make_protein_record(nuc_rec) for nuc_rec in records)
    with open("orf-aa-all.fasta","w") as output_handle:
        SeqIO.write(proteins, output_handle, "fasta")

def compareTrans(doc1,doc2):
    records1=saveFASTA(doc1)
    records2=saveFASTA(doc2)
    falserecs = []
    for r1 in records1:
        for r2 in records2:
            if(r1.id==r2.id):
                if(r1.seq!=r2.seq):
                    falserecs.append([r1.id,r2.id,r1.seq,r2.seq])
    return falserecs
#getFASTA("http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz")
#getFASTA("http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz")
#recs = saveFASTA("orf_coding_all.fasta")
#ARN = make_RNA_record(recs)
#translateFASTA(ARN)
#print(len(compareTrans("orf_trans_all.fasta","orf-aa-all.fasta")))

#####################################################################################
def findPloop(seq):
    carac = []
    for car in seq:
        carac.append(car)
    for i in range(len(carac)-7):
        if(carac[i]=="A" or carac[i]=="G"):
            if(carac[i+5]=="G"):
                if(carac[i+6]=="K"):
                    if(carac[i+7]=="S" or carac[i+7]=="T"):
                        return(True)
    if(carac[-8]=="A" or carac[-8]=="G"):
            if(carac[-3]=="G"):
                if(carac[-2]=="K"):
                    if(carac[-1]=="S" or carac[-1]=="T"):
                        return(True)
    return(False)

def PloopProteins(doc):
    records = saveFASTA(doc)
    count=0
    with open("Ploop_Proteins.txt","w") as f:
        for record in records:
            if(findPloop(record.seq)):
                f.write(record.id+'\n')
                count+=1
    #print(count)
#PloopProteins("orf-aa-all.fasta")

def histo(doc):
    with open("Ploop_Proteins.txt","r") as f:
        prot=f.readlines()
    for i in range(len(prot)):
        prot[i]=prot[i][:-1]
    records = saveFASTA(doc)
    dct = {}
    dctploop = {}
    for record in records:
        if(record.id in prot):
            for car in record.seq:
                if not(car in dctploop.keys()):
                    dctploop[car]=100/(len(record.seq)*len(prot))
                else:
                    dctploop[car]+=100/(len(record.seq)*len(prot))
        else:
            for car in record.seq:
                if not(car in dct.keys()):
                    dct[car]=100/(len(record.seq)*(len(records)-len(prot)))
                else:
                    dct[car]+=100/(len(record.seq)*(len(records)-len(prot)))

    for key in dct.keys():
        if key in dctploop.keys():
            print(key+" ploops | "+"p"*int(dctploop[key]*5)+" "+str(dctploop[key])+" %")
            print(key+" genome | "+"*"*int(dct[key]*5)+" "+str(dct[key])+" %")
        else:
            print(key+" genome | "+"*"*int(dct[key]*5)+" "+str(dct[key])+" %")

#histo("orf-aa-all.fasta")

############################################################################################"
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
amb = IUPACAmbiguousDNA()

def findRestr(doc):
    with open("Ploop_Proteins.txt","r") as f:
        prot=f.readlines()
    for i in range(len(prot)):
        prot[i]=prot[i][:-1]
    records = saveFASTA(doc)
    prots = []
    for i in range(len(records)):
        if records[i].id in prot:
            prots.append([records[i].id,records[i].seq])
    enzsites = []
    nonrestrprot = []
    for prt in prots:
        restr = [EcoRI.search(prt[1]), XhoI.search(prt[1]),TaqI.search(prt[1])]
        enzsites.append([prt[0],restr])
        for r in restr:
            if len(r)==0:
                nonrestrprot.append(prt[0])
                continue
    return nonrestrprot,enzsites
"""
prots,enzsites = findRestr("orf_coding_all.fasta")
print("Non restrictive proteins : ",prots)
for enz in enzsites:
    print("ID ",enz[0])
    print("EcoRI ",enz[1][0])
    print("XhoI ",enz[1][1])
    print("TaqI ",enz[1][2])"""
    
    
#######################################################################################################""
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

pdbl = PDBList()
pdbl.retrieve_pdb_file("2GAA")

def readPDBFile(filename):
    mmcif_dict = MMCIF2Dict(filename)
    nbchains,nbres,nbatoms,res = mmcif_dict['_struct_sheet.number_strands'],mmcif_dict['_struct_site.pdbx_num_residues'],mmcif_dict['_refine_hist.number_atoms_total'],mmcif_dict['_exptl.method']
    return sum([int(nbchains[i]) for i in range(len(nbchains))]), nbres, nbatoms, res
print(readPDBFile("ga/2gaa.cif"))
                    
            
            
            



    