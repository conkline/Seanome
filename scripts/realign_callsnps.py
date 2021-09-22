import sys
import os
import math
import re
import multiprocessing
import subprocess
import operator
import pandas as pd
import numpy as np
import pysam
import glob
from collections import defaultdict
from Bio.Sequencing import Ace
from Bio import SeqIO
from htb import htb

#get reverse complement of sequence
def revComp(dnaSeq):
    revCompBases ={"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join([revCompBases[x] for x in dnaSeq[::-1]])

#make cigar machine readable
def expandCigar(cigar):
    cigstr = ""
    cigar_re = re.compile(r'([0-9]*)([M=XID])')
    for m in cigar_re.finditer(cigar):
        if m.group(1):
            cigstr +=  m.group(2) * int(m.group(1))
        else:
            cigstr +=  m.group(2)
    return cigstr

#produce sequence from cigar
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
        break
    return newSeq

def compressCigar(cigar):
    """ 
    We build the cigar as a string of I, M, and D.  We need to do a RLE on this 
    to compact it. 
    """
    cigar = cigar.rstrip("D")
    c = None
    cnt = 0
    ciglst = []
    for b in cigar:
        if b != c:
            if c!= None:
                ciglst.append("%s%s"%(cnt, c))
            cnt = 1
            c = b
        else:
            cnt += 1
    if c != None:
        ciglst.append("%s%s"%(cnt, c))
    return "".join(ciglst)

def generateCigarFromSequence(sequence):
    """
    takes a sequence string (as aligned against a ref)
    and return a compressed cigar string
    """
    cigar = ""
    for c in sequence:
        if c == " ":
            cigar += 'I'
        else:
            cigar += 'M'
    return compressCigar(cigar)

def filterAlignment(locusin, alignment, contig, cigars, strands, paramList, popnames, readnames):
  
  #pull out params
  minDP = paramList[4]
  maxDP = paramList[5]
  minDPSample = paramList[6]
  new_alignment = ""
  new_contig = ""
  new_cigars = ""
  new_strands = ""
  new_readnames = ""
  new_popnames = ""
  
  #first, make sure we're not below min or above max depth
  if minDP <= alignment.shape[0] <= maxDP:
    #next, make sure each sample is above the cutoff
    uniq_pops = list(set(popnames))
    pop_counts = {x:popnames.count(x) for x in uniq_pops}
    below_cutoff = [x for x in pop_counts.keys() if pop_counts[x] < minDPSample]
    
    #if we remove some, may need to adjust alignment and contig accordingly
    if below_cutoff != []:
      reads_to_keep = [i for i, e in enumerate(popnames) if e not in below_cutoff] 
      new_alignment = alignment.take(reads_to_keep) #new alignment with dropped reads
      new_strands = [strands[x] for x in reads_to_keep]
      new_readnames = [readnames[x] for x in reads_to_keep]
      new_popnames = [popnames[x] for x in reads_to_keep]
      
      #remove all columns with all NAs and edit consensus
      na_counts = new_alignment.isnull().sum().tolist()
      cols_to_keep = [i for i, e in enumerate(na_counts) if e != new_alignment.shape[0]] 
      new_alignment = new_alignment.dropna(1, how = "all")
      contig_list = [char for char in contig]
      new_contig = "".join([contig_list[i] for i in cols_to_keep])
      
      #and edit cigar sequences
      new_cigars = []; test_reads = []
      for x in range(len(new_alignment)):
        read = new_alignment.iloc[x]
        read[pd.isna(read)] = " "
        tmpread = "".join(read.tolist()) 
        new_cigars.append(generateCigarFromSequence(tmpread))
        
    else:
      #otherwise, good to go
      new_alignment = alignment
      new_contig = contig
      new_cigars = cigars
      new_strands = strands
      new_readnames = readnames
      new_popnames = popnames
    
  if type(new_alignment) != str:
    new_alignment[new_alignment == ' '] = pd.np.nan
  return new_alignment, new_contig, new_cigars, new_strands, new_readnames, new_popnames

def callCCandSNP(locus, alignment, read_names, contig, pops, allpops, paramList):
    
    #parameters to consider
    #minor allele count, minimum depth, minimum number of samples, max # SNPs / locus
    
    bp_sim = []
    SNPinfo = []
    SNPdict = []
    genotypes = {}
    contig = list(contig)
    minCount = paramList[0]
    minDepth = paramList[1]
    minNSamples = paramList[2]
    maxSNPs = paramList[3]*len(contig) #given as a percentage of contig length
    cc_threshold = paramList[7]
    ploidy = paramList[8]
    ploidy_test_snp = ""
    ploidy_test_genotype = ""
    cc_test = ""
    mm_test = ""
    pos = []
    
    for c in range(len(alignment.columns)):
    
        nucs = alignment.iloc[:,c].value_counts()
        
        try:
            #add cc for this column
            col_bp = {x:nucs[x] for x in dict(nucs.dropna()) if x in ['A','T','G','C']}
            tmp_cc = float(max(col_bp.values())) / float(sum(col_bp.values()))
            bp_sim.append(tmp_cc)
    
            #make sure we have enough support for contig at this location (>1 read)
            #otherwise, convert to NNNNs
            if sum(list(col_bp.values())) == 1:
                contig[c] = "N"
    
            #call SNPs in this column - only worry about SNPs for now
            #TODO: indels
            if tmp_cc < 1: #if we have a polymorphism
		
                #see if we consider this a real SNP
                posCount = min(col_bp.values())
                posDepth = sum(col_bp.values())
                tmppops = list(np.array(pops)[np.array([not x for x in alignment.iloc[:,c].isnull()])]) #see which reads are here
                posSamples = len(set(tmppops))
    
                #check for minor allele count
                if posCount >= minCount:
                    #check for min depth
                    if posDepth >= minDepth:
                        #check for min # samples 
                        if posSamples >= minNSamples:
                            #TODO: move this into its own function
                            #if we pass all checks, record SNP
                            #LOCUS POS TYPE DP NS MAJOR MINOR1 MINOR2 MINOR3 A:C:G:T A:C:G:T etc...
                            #ex:
                            #ACU_34 72 18 2 A C - - 0.3 ACU 8:2:0:0 DAM 0:8:0:0
                            
                            if len(col_bp.keys()) == 2:
                              type_snp = "bisnp" #biallelic
                            else:
                              type_snp = "msnp" #multiallelic
                                                     
                            #set up major and minor alleles
                            sortedAlleles =  sorted(col_bp.items(), key=operator.itemgetter(1), reverse=True)
                            major = sortedAlleles[0][0]
                            minorOne  = sortedAlleles[1][0]
                            if len(sortedAlleles)==2:
                                minorTwo = "-"
                                minorThree = "-"
                            elif len(sortedAlleles)==3:
                                minorTwo = sortedAlleles[2][0]
                                minorThree = "-"
                            elif len(sortedAlleles)==4:
                                minorTwo = sortedAlleles[2][0]
                                minorThree = sortedAlleles[3][0]
    
                            #format allele counts by population
                            popdict = {x:[0,0,0,0] for x in allpops}
    
                            #for pop in list(set(pops)): #iterate through populations
                            for n in range(len(alignment.iloc[:,c])): #iterate through row
                                r = alignment.iloc[n,c]
                                if not pd.isnull(r):
                                    thispop = pops[n]
                                    if r=="A":
                                        popdict[thispop][0] += 1
                                    elif r=="C":
                                        popdict[thispop][1] += 1
                                    elif r=="G":
                                        popdict[thispop][2] += 1
                                    elif r=="T":
                                        popdict[thispop][3] += 1	
				
                            allelecounts = [":".join([str(y) for y in popdict[x]]) for x in sorted(popdict.keys())]
    
                            #reset ploidy test
                            ploidy_test_snp = "passed"
    
                            #filter by ploidy (per SNP)
                            for l in popdict.values():
                              if len([x for x in l if x > 0]) > ploidy:
                                ploidy_test_snp = "failed"
                            
                            if ploidy_test_snp == "passed":
                              #save position
                              pos.append(c)
                              
                              #generate SNP string
                              SNPinfo.append("\t".join([locus, str(c), type_snp, str(posDepth), str(posSamples), major, minorOne, minorTwo, minorThree] + allelecounts))
    
                              #and SNP dict
                              SNPdict.append({"chrom": locus, 
                                              "pos": c + 1,
                                              "depth": posDepth,
                                              "ns": posSamples,
                                              "ref": major,
                                              "alt": [x for x in [minorOne, minorTwo, minorThree] if x != "-"],
                                              "type": type_snp,
                                              "counts": [x[1] for x in sortedAlleles[1::]],
                                              "pops": list(set(tmppops)),
                                              "pop_dp": {x: tmppops.count(x) for x in tmppops},
                                              "pop_counts": popdict})
        
        except KeyError: #if we have all N's at this location, skip it
            continue
        except ValueError: #same
            continue
        except IndexError: #ignore 1-read overhang
            continue
    
    if np.mean(bp_sim) < cc_threshold:
      cc_test = "failed"
    else:
      cc_test = "passed"
    
    #if we have more SNPs than we expect in this contig, discard (probable sequence or alignment errors)
    if len(SNPinfo) >= maxSNPs:
        mm_test = "failed"
    else:
        mm_test = "passed"
        
    #if passes other checks, call haplotypes to check for undersplitting
    if mm_test == "passed" and cc_test == "passed":
      haplo_counts = callHaplos(locus, alignment, read_names, pos, minCount)
      if any(i > ploidy for i in haplo_counts.values()):
        ploidy_test_genotype = "failed"
      else:
        ploidy_test_genotype = "passed"
    
    return cc_test, ploidy_test_genotype, mm_test, SNPinfo, SNPdict, "".join(contig)

def locusRealign(locusin, paramList, outPath):
    #send subprocess to run cap3 alignment 
    command = "/10tb_abyss/emily/seanome_project/programs/CAP3/cap3 " + outPath + locusin + ".fasta -g 1 -i 21 -j 31 -o 16 -p 87 -s 251 -k 0 -z 1 -y 6 > " + outPath + locusin + ".cons"
    process = subprocess.Popen(command, shell=True)
    process.wait() #wait for alignment to finish
    
    #min and max depth for alignment to be considered
    minDP = paramList[4]
    maxDP = paramList[5]
    
    alignment = ""
    contig = ""
    strands = ""
    cigars = ""
    popnames = ""
    readnames = ""
    result = "passed"
    
    try:
        #read in Ace format file
        acetmp = Ace.read(open(outPath + locusin + ".fasta.cap.ace"))
        c = acetmp.contigs[0]
        contig = c.sequence
        reads = []
        readnames = []
        strands = []
        cigars = []
        flag = True
        
        #initial check of alignment
        if len(acetmp.contigs) > 1: #split into multiple contigs
          alignment = ""
          result = "failed"
          flag = False
        if contig.count("*") > len(contig)/2: #failure to align
          alignment = ""
          result = "failed"
          flag = False
        if c.nreads < minDP or c.nreads > maxDP: #too shallow or too deep
          alignment = ""
          result = "depth"
          flag = False          

        if flag:
          #extract info from ace
          for readnum in range(c.nreads):
              tmpread = " "*(c.af[readnum].padded_start-1) + c.reads[readnum].rd.sequence.replace("*"," ")
              tmpname = c.af[readnum].name
              tmpstrand = c.af[readnum].coru
              if tmpstrand=="C":
                  tmpstrand = "-"
              elif tmpstrand=="U":
                  tmpstrand = "+"
              strands.append(tmpstrand)
              reads.append(tmpread)
              readnames.append(tmpname)
              cigars.append(generateCigarFromSequence(tmpread)) #convert seq to cigar

          #convert to pd dataframe 
          alis = [pd.Series(list(x)) for x in reads]
          alignment = pd.DataFrame(alis)
          alignment[alignment == ' '] = pd.np.nan

          #just keep populations from readnames
          popnames = [x.split("_F_")[0] for x in readnames]
    
    except (IndexError, ValueError, AttributeError) as e: 
        #failure to align
        alignment = ""
        contig = ""
        strands = ""
        cigars = ""
        popnames = ""
        readnames = ""
        result = "failed"
    
    #delete cons files
    toRemove = [".cons", ".fasta.cap.ace", ".fasta.cap.contigs", ".fasta.cap.contigs.links", ".fasta.cap.contigs.qual", ".fasta.cap.info", ".fasta.cap.singlets"]
    for i in toRemove:
        os.remove(outPath + locusin + i)
    
    #output alignment, contig, strands, cigars, reads
    return alignment, contig, strands, cigars, popnames, readnames, result

def getAlignment(locusin, tmp_read_dict, tmp_mapping, paramList):
    
    #min and max depth for alignment to be considered
    minDP = paramList[4]
    maxDP = paramList[5]
    
    alignment = ""
    read_strands = ""
    read_cigars = ""
    read_popnames = ""
    read_names = ""
    result = "passed"
    
    if minDP <= len(tmp_read_dict.keys()) <= maxDP:
      
      #get read info
      read_names = [x[0] for x in tmp_mapping]
      read_popnames = [x.split("_F_")[0] for x in read_names]
      read_cigars = [x[1] for x in tmp_mapping]
      read_strands = [x[2] for x in tmp_mapping]
      read_seqs = []
    
      for x in range(0, len(read_names)):
        tmp_cig = read_cigars[x]
        tmp_seq = str(tmp_read_dict[read_names[x]].seq)
        tmp_strand = read_strands[x]
        new_seq = getSeqFromCigar(tmp_seq, tmp_cig, tmp_strand)
        read_seqs.append(new_seq)
    
      #convert to pd dataframe
      alis = [pd.Series(list(x)) for x in read_seqs]
      alignment = pd.DataFrame(alis)
      alignment[alignment == ' '] = pd.np.nan
      
    else:
      result = "depth"
    
    #output alignment, contig, strands, cigars, reads
    return alignment, read_strands, read_cigars, read_popnames, read_names, result
    
#parallelize filtering
def worker(input_dir, output_dir, chunkedContigs, chunkedMapping, paramList, allpops, m_handle, c_handle, s_handle, b_handle, bam_header, v_handle, vcf_header):
    
    mappingFile = open(m_handle, "a")
    contigFile = open(c_handle, "a")
    snpFile = open(s_handle, "a")
    bamFile = pysam.AlignmentFile(b_handle, "wb", header=bam_header)
    outPath = output_dir + "/ALL/"
    
    bamOut = []
    vcfOut = []
    
    for locus in chunkedContigs.keys():
        print(locus)
        cap_result = ""
        sea_result = ""
        
        #pull reads from fastq
        command = "fgrep " + locus + " " + outPath + "ALL_FINAL.mapping_to_cons | cut -f1 | sed 's/$/\$/g' | grep --no-group-separator -f - -A3 " + outPath + "ALL.fastq > " + outPath + locus + ".fastq"
        command2 = "/10tb_abyss/emily/seanome_project/seqtk/seqtk seq -A " + outPath + locus + ".fastq > " + outPath + locus + ".fasta"
    
        process = subprocess.Popen(command, shell=True)
        process.wait() 
        process = subprocess.Popen(command2, shell=True)
        process.wait()
        
        #import mapping and reads
        read_dict = SeqIO.to_dict(SeqIO.parse(output_dir + "/ALL/" + locus + ".fastq", "fastq"))
      
        #first, realign reads and return alignment, filter based on depth
        alignment, contig, strands, cigars, popnames, readnames, cap_result = locusRealign(locus, paramList, outPath)
        
        #if cap3 can't realign, use original alignment, filter based on depth
        if cap_result == "failed":
          alignment, strands, cigars, popnames, readnames, sea_result = getAlignment(locus, read_dict, chunkedMapping[locus], paramList)
          contig = chunkedContigs[locus][1].upper()
        
        #take valid cap3 or seanome alignment and do further QC, call SNPs
        if cap_result == "passed" or sea_result == "passed":
          final_alignment, final_contig, final_cigars, final_strands, final_readnames, final_popnames = filterAlignment(locus, alignment, contig, cigars, strands, paramList, popnames, readnames)
          
          #if passed final checks, call SNPs and check for mismatch %, CC, and ploidy
          if type(final_alignment) != str:
            cc_test, ploidy_test, mm_test, SNPlist, SNPdict, strcontig = callCCandSNP(locus, final_alignment, final_readnames, final_contig, final_popnames, allpops, paramList) 
            
            #if we retain this locus, build output files
            if cc_test == "passed" and mm_test == "passed" and ploidy_test == "passed":
                mapping = ""
                for n in range(len(final_readnames)):
                    tmpname = final_readnames[n]
                    tmpcigar = final_cigars[n]
                    tmpstrand = final_strands[n]
                    mapping += tmpname + "\t" + locus + "\t" + tmpcigar + "\t" + tmpstrand + "\n"
                  
                #write fasta + mapping tsv
                #strcontig = "".join(n_contig)
                contigFile.write(">" + locus + "\n" + strcontig + "\n")
                mappingFile.write(mapping)
                #write BAM
                bamOut.append([chunkedContigs[locus][0], read_dict, final_cigars, final_readnames, strcontig])
                
                if SNPlist != []:
                  vcf_header.add_meta('contig', items=[('ID', locus)])
                  snpFile.write("\n".join(SNPlist)+"\n")
                  vcfOut.append([locus, SNPdict])
        
        #clean up files
        os.remove(outPath + locus + ".fastq")
        os.remove(outPath + locus + ".fasta")
    
    #TEST - write to BAM and VCF
    if bamOut != []:
      writeBAM(bamFile, bamOut)
      
    vcfFile = pysam.VariantFile(v_handle, "wu", header=vcf_header)
    if vcfOut != []:
      for i in range(len(vcfOut)):
        locus = vcfOut[i][0]
        snps = vcfOut[i][1]
        
        for snp in snps:
          #TODO - calculate quality (QUAL and GQ)
          #create new record
          r = vcfFile.new_record(contig=snp["chrom"], start=snp["pos"]-1, stop=snp["pos"], alleles=[snp["ref"]] + snp["alt"], id="_".join([snp["chrom"], str(snp["pos"])]))
          
          #add info            
          r.info["NS"] = snp["ns"]
          r.info["DP"] = int(snp["depth"])
          r.info["AC"] = ",".join([str(x) for x in snp["counts"]])
          r.info["TYPE"] = snp["type"]
          
          #add sample info
          for pop in snp["pops"]:
            popcounts, haplos = dictToCounts(snp["pop_counts"], pop, snp["ref"], snp["alt"])
            bytes_pop = pop.encode('utf-8')
            r.samples[bytes_pop]["DP"] = snp["pop_dp"][pop]
            r.samples[bytes_pop]["AD"] = ",".join([str(x) for x in popcounts])
            r.samples[bytes_pop]["GT"] = "/".join([str(x) for x in haplos])
            r.samples[bytes_pop]["GQ"] = "." 
          
          vcfFile.write(r)
          
    #writeVCF(vcfFile, vcf_header, vcfOut)
    mappingFile.close()
    contigFile.close()
    snpFile.close()
    bamFile.close()
    vcfFile.close()
    
    return

def main():

    #TODO: 
    #write to VCF
    #print out summary stats
    
    #get parameters
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    threads = sys.argv[3]
    minDP = sys.argv[4]
    maxDP = sys.argv[5]
    sampleMin = sys.argv[6]
    minCC = sys.argv[7]
    minCount = sys.argv[8]
    minDepth = sys.argv[9]
    minNSamples = sys.argv[10]
    maxMismatch = sys.argv[11]
    ploidy = sys.argv[12]
    
    #input_dir = "/10tb_abyss/emily/snail_seanome/2021_paolo_seanome"
    #output_dir = "/10tb_abyss/emily/snail_seanome/2021_paolo_seanome_indexed"
    os.mkdir(output_dir + "/FINAL_OUTPUT")
    #os.chdir(output_dir)
    
    #set filtering parameters
    #minDP = 3
    #sampleMin = 3
    #maxDP = 600 
    #minCC = .95
    #threads = 30
    
    #set SNP parameters
    #minCount = 2
    #minDepth = 6
    #minNSamples = 1
    #maxMismatch = .15
    #ploidy = 2
    
    #TODO: popmap support
    allpops = list(set([x.split(".")[0].split("/")[-1] for x in glob.glob(input_dir + "/*") if x.endswith("fastq") or x.endswith("fq")]))
    
    #read in input (fastq)
    #all_reads = SeqIO.to_dict(SeqIO.parse(output_dir + "/ALL/ALL.fastq", "fastq"))
    
    #read in input (mapping and contigs)
    seanome_mapping = open(output_dir + "/ALL/ALL_FINAL.mapping_to_cons", "r") #input
    contigsTmp = SeqIO.to_dict(SeqIO.parse(output_dir + "/ALL/FINAL.final.contigs", "fasta"))
    contigsList = sorted(list(contigsTmp.keys()))
    contigsDict = defaultdict(list) #for making BAM file
    counter = 0
    for x in contigsList:
        contigsDict[x] = [counter, str(contigsTmp[x].seq)] #[ID sequence]
        counter += 1
    
    #build depth data for clustering if no max dp specified
    depths = defaultdict(int)
    
    #set up mapping dict
    mappingDict = defaultdict(list)
    for line in seanome_mapping.readlines():
      myline = line.rstrip().split("\t") 
      depths[myline[1]] += 1
      mappingDict[myline[1]].append([myline[0], myline[2], myline[3]]) #[read cigar strand]
        
    if maxDP == "auto":  
      #define max depth using head-tail breaks, keeping ht-index of 4
      maxDP = htb(depths.values())[3]
      
    #consolidate parameters
    paramList = [minCount, minDepth, minNSamples, maxMismatch, minDP, maxDP, sampleMin, minCC, ploidy]
  
    #set up BAM header
    bam_header = {'HD': {'VN': '1.0', 
                  'SO': 'unsorted'},
                  'PG': [{'ID': 'Seanome'}],
                  'SQ': [{'LN': len(contigsDict[x][1]), 'SN': x} for x in contigsDict.keys()],
                  'RG': [{'ID':x} for x in allpops]}
                  
    #set up VCF header
    vcf_header = pysam.VariantHeader()
    vcf_header.add_meta('INFO', items=[('ID',"NS"), ('Number',1), ('Type','Integer'), ('Description','Number of Samples With Data')])
    vcf_header.add_meta('INFO', items=[('ID',"DP"), ('Number',1), ('Type','Integer'), ('Description','Total Depth')])
    vcf_header.add_meta('INFO', items=[('ID',"AC"), ('Number',1), ('Type','String'), ('Description','Counts of Alternate Alleles')])
    vcf_header.add_meta('INFO', items=[('ID',"TYPE"), ('Number',1), ('Type','String'), ('Description','Variant Type (bisnp, msnp)')])
    
    vcf_header.add_meta('FORMAT', items=[('ID',"DP"), ('Number',1), ('Type','Integer'), ('Description','Read Depth Per Sample')])
    vcf_header.add_meta('FORMAT', items=[('ID',"AD"), ('Number',1), ('Type','String'), ('Description','Allele Counts Per Sample')])
    vcf_header.add_meta('FORMAT', items=[('ID',"GT"), ('Number',1), ('Type','String'), ('Description','Genotype')])
    vcf_header.add_meta('FORMAT', items=[('ID',"GQ"), ('Number',1), ('Type','String'), ('Description','Genotype Quality')])
    
    for x in allpops:
      vcf_header.add_sample(x)
    
    #run realignment and filtering in parallel
    numKeys = len(contigsList)
    chunkSize = int(math.ceil(numKeys / threads))

    if __name__ == '__main__':
        currentLine = 0 #normally set to 0 except for troubleshooting
        jobs = []
        for i in range(threads):
            m_handle = output_dir + "/FINAL_OUTPUT/mapping" + str(i) + ".tsv"
            c_handle = output_dir + "/FINAL_OUTPUT/contigs" + str(i) + ".fasta"
            s_handle = output_dir + "/FINAL_OUTPUT/snp" + str(i) + ".tsv"
            b_handle = output_dir + "/FINAL_OUTPUT/raw_alignment" + str(i) + ".bam"
            v_handle = output_dir + "/FINAL_OUTPUT/raw_snps" + str(i) + ".vcf"
            
            try:
                dictChunk = contigsList[currentLine:(currentLine+chunkSize)]
            except IndexError:
                dictChunk = contigsList[currentLine:(numKeys+1)] #normally set to numKeys+1
            currentLine += chunkSize
            
            #break contigs and mapping dictionaries into chunks
            chunkedContigs = {x:contigsDict[x] for x in dictChunk} #[ID, sequence]
            chunkedMapping = {x:mappingDict[x] for x in dictChunk} #[read, cigar, strand]
            
            #send jobs to multiprocessing
            p = multiprocessing.Process(target=worker, args=(input_dir, output_dir, chunkedContigs, chunkedMapping, paramList, allpops, m_handle, c_handle, s_handle, b_handle, bam_header, v_handle, vcf_header,))
            jobs.append(p)
            p.start()
            
main()
