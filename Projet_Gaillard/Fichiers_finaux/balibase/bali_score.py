#!/usr/bin/python
# Calculate BALI scores from a test and a reference multiple sequence alignment in FASTA format
# reproduces the results of the bali_score C program
# Thomas Gaillard, March 2018
# ./bali_score.py ref.fasta test.fasta

from __future__ import print_function

import sys

reffile = sys.argv[1]
testfile = sys.argv[2]

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
