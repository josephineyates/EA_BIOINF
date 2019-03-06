
# coding: utf-8

# # Bioinfo JM SSteyaert
# 
# ## Détection de signatures sur un génome complet

# In[32]:


import requests

url1 = "ftp://ftp.expasy.org/databases/prosite/prosite.dat"
url2 = "http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz"

def load(url, name):
    r = requests.get(url)
    with open(name, 'wb') as handle:
        handle.write(r.content)

load(url2, "orf_coding_all.fasta.gz")


# In[36]:


from Bio.ExPASy import Prosite
handle = open("myprosite.dat")
records = Prosite.parse(handle)
for record in records:
    if "p-loop" in record.description.lower():
        print(record.name)
        print(record.accession)
        print(record.pattern)


# In[47]:


import shutil
import gzip

with gzip.open("orf_coding_all.fasta.gz", 'rb') as f_in, open("orf_coding_all.fasta", 'wb') as f_out:
    shutil.copyfileobj(f_in, f_out)


# In[40]:


import os
print("Taille du fichier : {} ko".format(os.stat("orf_coding_all.fasta").st_size/1000))


# In[54]:


from Bio import SeqIO
from Bio.Seq import Seq

handle = open("orf_coding_all.fasta")
seqRec = []
for seq in SeqIO.parse(handle, "fasta"):
    seqRec.append(seq)
seqRec


# In[61]:


seqTr = []
for seq in seqRec:   
    new_s = seq.seq.translate()
    new_seq = SeqRecord(new_s, id=seq.id, name=seq.name, description=seq.description, dbxrefs=seq.dbxrefs)
    seqTr.append(new_seq)

print(len(seqTr))

with open("orf_aa_all.fasta", "w") as output_handle:
        SeqIO.write(seqTr, output_handle, "fasta")


# In[75]:


url3 = "http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz"
load(url3, "orf_trans_all.fasta.gz")

with gzip.open("orf_trans_all.fasta.gz", 'rb') as f_in, open("orf_trans_all.fasta", 'wb') as f_out:
    shutil.copyfileobj(f_in, f_out)


# In[76]:


handle = open("orf_trans_all.fasta")
seqAA = []
for seq in SeqIO.parse(handle, "fasta"):
    seqAA.append(seq)
len(seqAA)


# In[79]:


for k in range(6713):
    if seqTr[k].seq != seqAA[k].seq:
        print(k)
        print(seqTr[k].seq, seqAA[k].seq)


# On note étrangement des différences pour les séquences 6681 à 6708 (inclues) entre la traduction depuis la séquence ADN et la séquence protéique correspondante. Le code génétique de la mitochondrie n'est-il pas différent du code génétique universel ?

# In[82]:


def transition(state, letter):
    if state == 0 and letter in "AG":
        return 1
    if state in [1, 2, 3, 4]:
        return state+1
    if state == 5 and letter == "G":
        return 6
    if state == 6 and letter == "K":
        return 7
    if state == 7 and letter in ["S", "T"]:
        return 8

def automate_search(seq):
    states = [0]
    for letter in seq:
        new_states = [0]
        for state in states:
            r = transition(state, letter)
            if r == 8:
                return True
            if r != None and !(r in new_states):
                new_states.append(r)
        states = new_states
    return False

output = open("P-loop_prot_list.txt", "w")
for seqaa in seqAA:
    if automate_search(seqaa.seq):
        output.write("{}\n".format(seqaa.id))
output.close()


# In[87]:


count_p = {}
count = {}
total = 0
total_p = 0

def add(seq, dic):
    for l in seq:
        if not l in dic:
            dic[l] = 0
        dic[l] += 1

for seqaa in seqAA:
    total += len(seqaa)
    add(seqaa.seq, count)
    if automate_search(seqaa.seq):
        total_p += len(seqaa)
        add(seqaa.seq, count_p)
        
for aa in count:
    count[aa] = 100. * float(count[aa]) / total
    
for aa in count_p:
    count_p[aa] = 100. * float(count_p[aa]) / total_p
    
print(len(count_p))
    
for aa in count_p:
    print("{} ploops | ".format(aa) + "p"*int(count_p[aa]*5) + " {:.3}%".format(count_p[aa]))
    print("{} genome | ".format(aa) + "*"*int(count[aa]*5) + " {:.3}%".format(count[aa]))


# In[111]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')

fig = plt.figure(figsize=(15, 10))
ax = fig.add_subplot(111)

## the data
N = 21
aa = list(count_p.keys())
p = [count_p[a] for a in aa]
tot = [count[a] for a in aa]

## necessary variables
ind = np.arange(N)                # the x locations for the groups
width = 0.45                      # the width of the bars

## the bars
rects1 = ax.bar(ind, p, width,
                color='darkblue')

rects2 = ax.bar(ind+width, tot, width,
                    color='darkorange')

# axes and labels
ax.set_xlim(-width,len(ind)+width)
ax.set_ylim(0,10)
ax.set_ylabel('Percentages')
ax.set_title('Percentage of genome')
xTickMarks = aa
ax.set_xticks(ind+width/2)
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=45, fontsize=10)

## add a legend
ax.legend( (rects1[0], rects2[0]), ('P-Loops', 'Total genome') )

plt.show()

