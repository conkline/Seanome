#!/usr/bin/env python2.7 

import argparse
import time
import logging
import multiprocessing
import itertools
from collections import defaultdict
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from utils.cigar import compressCigar
import sys

logging.basicConfig(level=logging.DEBUG,
                    format='(%(threadName)-9s) %(message)s',)

# Add this as a param
# minimum number of non ambiguous bases (not N ) that a seed needs to have
# or it is trashed
minNonNs = 70
maxFilteredBases = 70

def processWithCap3(msa):
    seqs = []
    # write the ungapped msa sequences to a file that will be assembled using cap3
    for s in msa:
        if s.id != "consensus":
            s.seq = s.seq.ungap("-")
            s.id = s.id.replace("*", "")
            seqs.append(s)
    SeqIO.write(seqs, open(os.path.join("cap3",seqs[0].id), 'w'), "fasta")

    # run Cap3 on sequeneces
    runCap3 =  shlex.split("/home/mahdi/programs/CAP3/cap3 %s -z 3 -j  31 -o  16 -p  77  -r  1 -s  251 -h 20" % os.path.join("cap3", seqs[0].id))
    with open(os.devnull, 'w') as FNULL:
        subprocess.call(runCap3, stdout=FNULL)

    # parse the cap3 ali
    # We are  only reading one contig (largest) and ignoring rest
    ace_gen = Ace.parse(open(os.path.join("cap3","%s.cap.ace" % seqs[0].id), 'r'))
    contig = ace_gen.next()

    # last sequence starts at THIS CAN BE NEGATIVE
    firstStart = min([contig.af[i].padded_start for i in range(len(contig.reads))])
    # last sequence ends at
    lastEnd = max([contig.af[i].padded_start + len(contig.reads[i].rd.sequence) for i in range(len(contig.reads))])
    aliReads = []
    for i in range(len(contig.reads)):
        seEffectiveLen = contig.af[i].padded_start + len(contig.reads[i].rd.sequence)
        seqStr = "-"* (contig.af[i].padded_start - firstStart) + contig.reads[i].rd.sequence.replace("*", "-")+"-"* (lastEnd - seEffectiveLen)
        aliReads.append(SeqRecord(Seq(seqStr, generic_dna), id=contig.af[i].name))
    for s in aliReads:
        pass
    # print s.id
        # print str(s.seq)

    for fileName in glob.glob(os.path.join("cap3", "%s*" % seqs[0].id)):
        # print "deleting file: %s " % fileName
        os.remove(fileName)

    return  MultipleSeqAlignment(aliReads)



def isProblematicAli(ali, maxErrNonAli=0):
    
    refSeqStr = str(ali[0].seq)

    startPos = min(set([refSeqStr.find("A"), refSeqStr.find("C"), \
                            refSeqStr.find("G") , refSeqStr.find("T"), \
                            refSeqStr.find("N")]).difference(set([-1])))

    refSeqStr = refSeqStr[::-1]
    endPos = len(refSeqStr) -  min(set([refSeqStr.find("A"), refSeqStr.find("C"), \
                                            refSeqStr.find("G") , refSeqStr.find("T"), \
                                            refSeqStr.find("N")]).difference(set([-1])))

    nbSeqsWithUpsGap = 0 # seqs with upstream gaps
    nbSeqsWithDownsGap = 0 # seqs with downstream gaps


    for s in ali:
        if s.id == "consensus":
            continue

        # Checking in the end for a gap starting righ before the start of the ref (upstream)
        beforeStart = str(s.seq[:startPos])[::-1]
        possibleStarts = set([beforeStart.find("A"), \
                                  beforeStart.find("C"), beforeStart.find("G") , \
                                  beforeStart.find("T"), beforeStart.find("N")]).difference(set([-1]))

        if len(possibleStarts) == 0:
            continue
        gapSize = min(possibleStarts)
        if gapSize > 4:
            nbSeqsWithUpsGap +=1

        # Checking in the end for a gap starting righ after the end of the ref (downstream)
        afterEnd = str(s.seq[endPos:])
        possibleStarts = set([afterEnd.find("A"), \
                                  afterEnd.find("C"), afterEnd.find("G") , \
                                  afterEnd.find("T"), afterEnd.find("N")]).difference(set([-1]))

        if len(possibleStarts) == 0:
            # all gaps from here on... strange but possible
            continue
        gapSize = min(possibleStarts)
        if gapSize > 4:
            nbSeqsWithDownsGap +=1

    if (nbSeqsWithUpsGap > maxErrNonAli) or (nbSeqsWithDownsGap > maxErrNonAli):
        return True

    return False

def getClustersIter(fastaFile):
    """
    Creates a generator from the MSA produced by usearch/vsearch.  Each
    MSA record ends with a consensus.  Use this consensus sequence to
    identify when a MSA record ends.
    """
    accumulator= []
    for seq in SeqIO.parse(fastaFile, "fasta"):
        accumulator.append(seq)
        if seq.id.lower().startswith("consensus"):
            yield MultipleSeqAlignment(accumulator)
            accumulator = []
    if accumulator:
        yield MultipleSeqAlignment(accumulator)

def generateConsSequence(msa, tag, threshold=0):
    # build the new cons sequence

    si = AlignInfo.SummaryInfo(msa)

    conSeq = si.dumb_consensus(threshold = threshold, ambiguous = 'N')
    consId =  msa[0].id.replace("*","")  + "_" + tag
    newCons = SeqRecord(conSeq, id = consId, description="")
    return newCons

def generateCigarFromSequence(sequence):
    """
    takes a sequence string (as aligned against a ref)
    and return a compressed cigar string
    """
    cigar = ""
    for c in sequence:
        if c == "-" or c == '+':
            cigar += 'I'
        else:
            cigar += 'M'
    return compressCigar(cigar)




def processCluster(data):
    # Takes one cluster msa info (data[1])
    # geneates the new cons either from vseach or from
    # Cap3 and updates the reads information
    # a cluster info in vsearch start with ">*" and end with a "consensus" sequence and an
    # empty line
    try:
        args = data[0]
        msa = data[1]
        oldSeed_id = msa[0].id.replace("*", "")
        
        newCigs = {}
        newCons = ""

        # the Seed is mostly Ns, so we can discard it.        
        if len(str(msa[0].seq).upper().replace("N", "")) < minNonNs:
            return (oldSeed_id, None)

        # if the sequence is mostly lower cases (filtered bases), then we also discard it
        if len(filter(str.islower, str(msa[0].seq))) >maxFilteredBases or \
                len(filter(str.islower, str(msa[0].seq))) > 0.8 * len( str(msa[0].seq) ):
            return (oldSeed_id, None)



        if isProblematicAli(msa):
            # pass
            # take the msa generated by CAP3 instead
            # msa = processWithCap3(msa)
            return (oldSeed_id, None)
        else:
            # print "%s: Not problematic" % msa[0].id
            pass
         
        # seqRecord consensus sequence generate either form vsearch output or from CAP3
        newCons =  generateConsSequence(msa, args.tag)
        
        newCigs =  dict([ \
                (s.id.replace("*", ""),[s.id.replace("*", ""), newCons.id, generateCigarFromSequence(str(s.seq))]) \
                    for s in msa \
                    ])
        if "consensus" in newCigs:
            del(newCigs["consensus"])
        
        return (newCons, newCigs)

    except KeyboardInterrupt:
        print "Interrupt called: not cleaning up!"
        



def parseOutFile(fileName):
    clustInfoLines =defaultdict(list)
    for line in open(fileName):
        data = line.rstrip().split()
        clustInfoLines[data[1]].append([data[0], data[-2], data[-1]])
    return clustInfoLines




if __name__ == "__main__":

    startTime  = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threads', type = int, default = 2, help = "Number of processing threads. (default: 2)" )
    parser.add_argument('-i1', '--input1',  required = True, help = "MSA output from vsearch/usearch")
    parser.add_argument('-o1', '--output1', required = True, help = " MSA consensus output")
    parser.add_argument('-i2', '--input2',  required = True, help = "userout from vsearch/usearch")
    parser.add_argument('-o2', '--output2', required = True, help = "modified userout")
    parser.add_argument('-g', '--tag', required = True, help = "what suffix to add the new consensus sequence's id. Ex. a clustering iteration number")
    parser.add_argument('-r', '--problemCtgs', required = True, help = "trcks problematic clusters that were assembled with CAP3")
    # parser.add_argument('-n', '--largecluster', required = False, type = int, default = 0,
    #                     help = "The size of a cluster that we consider too large, anything above this will be dropped")
    args = parser.parse_args()


    outLines = parseOutFile(args.input2)
    pool = multiprocessing.Pool(args.threads)
    logging.basicConfig(level=logging.DEBUG,
                        format='(%(threadName)-9s) %(message)s',)

    consensuses =  open(args.output1, "w")
    newOut =  open(args.output2, "w")
    probClusters =  open(args.problemCtgs, "w")

    # try:

    # FOR DEBUG, you can SWAP the following two lines with the next two
    # for msa in getClustersIter(args.input1):
    #     newCons, newCigs = processCluster([args, msa])
    for output in pool.imap(processCluster,  itertools.izip(itertools.repeat(args), getClustersIter(args.input1))):

        # print output
        #sys.exit(0)

        newCons, newCigs = output
        
        if newCigs == None: # for now, problematic CAP3 are not processed ... make it as a param
            print >> probClusters, newCons
            continue
        oldConsId = "_".join(newCons.id.split("_") [:-1])
 
        # write the new consensus sequence
        SeqIO.write(newCons, consensuses, "fasta")    

        consOutLines = outLines[oldConsId]
        # write the new out file (args.input2)
        # Add the strand as is to what was alrady generate in newCigs
        for seqInfo in consOutLines:

            if seqInfo[0] in newCigs:
                newCigs[seqInfo[0]].append(seqInfo[-2])
                newCigs[seqInfo[0]].append(seqInfo[-1])
                print >> newOut, "\t".join(map(str,  newCigs[seqInfo[0]]))
                del(newCigs[seqInfo[0]])
        # write the seed against itself


        print >> newOut, "\t".join(map(str,  newCigs[oldConsId]+["+", "+"]))
        del(newCigs[oldConsId])
        
        if len(newCigs) > 0:
            # whatever is left is either a seed or a singlet that does not fit
            print "WARNING: newCigs is not empty. It contains: %s" % newCigs.keys()
        


    # except:
    # print "Error in the pool"
    # finally:
    pool.close()
    consensuses.close()
    newOut.close()
    probClusters.close()

    endTime = time.time()
    print "Total run time is: %s" % (endTime - startTime)
