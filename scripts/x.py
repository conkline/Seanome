#!/usr/bin/env python2.7
# TODO
# takes the clustering results from stage 1 and from stage 2 and 
# by transitivity, finds the alignments between the initial reads the contigs in stage 2.
# ex in stage 1; seq_1 contig_1 100M + + 
#    in stage 2; contig_1 contig_1_2 100M + + 
#    in stage 3; seq_1 contig_1_2 100M + + 



# Filter the clusters that are rare or too abundant.
# Use a flag here as well to know whether we are working with an individual or with pool of individuals
# If mode is individual, filter parse each alignemnt and drop alignment where column do no spport diploidy (if diploid)



import re
import sys
from collections import defaultdict
from Bio import SeqIO
import argparse

from utils.cigar import expandCigar, compressCigar

# def getLociInfo(consFastaFile, lociInfo):
#     # will contain a dict with key is consensus and and list of lines
#     # from out files of reads contained in that sequence
#     accumulator= {}     
#     for seq in SeqIO.parse(consFastaFile, "fasta"):
#         accumulator[seq.id] = lociInfo[seq.id]
#         yield accumulator
#         accumulator = {}




def revComp(dnaSeq):
    revCompBases ={"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join([revCompBases[x] for x in dnaSeq[::-1]])

def getSeqFromCigar(dnaSeq, cigarSeq, strand):
    pos = 0
    newSeq = ""
    expandedCigar = expandCigar(cigarSeq)

    if strand == "-":
        dnaSeq = revComp(dnaSeq)
        
    for char in expandedCigar:
        if char == "I":
            newSeq+= " "
        elif char== "M":
            newSeq+=dnaSeq[pos]
            pos+=1
        else:
            sys.exit("UNEXPECTED CHAR %s" % char)
    return newSeq



def updateCigar2(it1_cigar, it2_cigar, it1_seq_strand, it2_seq_strand):


    it2_seq_strand_reverse = False
    if it2_seq_strand =="-":
        it2_seq_strand_reverse = True


    myCigar = it1_cigar
    orignewcig = it2_cigar

 
    if it2_seq_strand_reverse:
        myCigar = it1_cigar[::-1]

    posCigar = 0
    posRef = 0
    posRead = 0
    newCig = ""
    while posRef < len(it2_cigar) and posRead < len(myCigar):
        if it2_cigar[posRef] == 'I':
            newCig += 'I'
            posRef += 1
        elif it2_cigar[posRef] == 'D':
            # if we come across a 'D' in the new cigar, we need to keep it
            newCig += 'D'
            posRef += 1
        elif myCigar[posRead] == 'I':
            newCig += 'I'
            posRead += 1
            posRef += 1
        else: 
            newCig += 'M'
            posRef += 1
            posRead += 1
    while posRef < len(it2_cigar):
            newCig += "I"
            posRef += 1
    
    return newCig

def getStrand(or1, or2):
    # if we have two orientation (+,-) or (-,+) we return "-"
    # if we have one orientation (+,+) or (-,-) we return "+"
    # (-,-) means read is + by transitivity

    if len(set([or1, or2])) > 1:
        return "-"
    else:
        return "+"


if __name__ == "__main__":
    lvl1_clusts = defaultdict(list)
    print >> sys.stderr, "starting"
    # sys.arg[1]: samples.mapping_to_cons (concatenated : POP1S1.mapping_to_cons, POP1S1.mapping_to_cons, ... )
    for line in open(sys.argv[1]):
        data = line.split()
        lvl1_clusts[data[1]].append([ data[0], data[2], getStrand(data[3], data[4]) ])
    print >> sys.stderr, "finished loading the mapping"

    # sys.arg[2]: ALL.mapping_to_cons
    for line in open(sys.argv[2]):
        data = line.split()
        if data[0] in lvl1_clusts:
            lvl2_cigar = data[2]
            for entry in lvl1_clusts[data[0]]:
                lvl1_cigar = entry[1]
                lvl2_strand = getStrand(data[3], data[4])
                newCigar = compressCigar(updateCigar2(expandCigar(lvl1_cigar), expandCigar(lvl2_cigar), entry[2], lvl2_strand))
                print "%s\t%s\t%s\t%s" % (entry[0], data[1], newCigar, getStrand(entry[2], lvl2_strand))


                


