from Bio import SeqIO
from collections import defaultdict
import scipy.sparse as sp
import sys
from difflib import SequenceMatcher

#TODO: support for ddRAD

#get all parameters
dir_in = sys.argv[1]; dir_out = sys.argv[2] #input and output directories
pop_name = sys.argv[3] #population name
paired = sys.argv[4] #whether samples have paired reads
mixed = sys.argv[5] #if samples have different res. enzymes or pair status, csv file describing them
nb_output_files = int(sys.argv[6]) #number output files per sample (default: 4)
hash_len = int(sys.argv[7]) #hash length in bp (default: 6)
cutsite_seq = sys.argv[8] #expected cutsite associated with enzyme (default: GATC)
cutsite_seq_two = sys.argv[9] #second expected cutsite associated with enzyme (default: off)
keep_strand = sys.argv[10] #whether to keep a read that failed in one direction but not the other (default: off)
rescue = sys.argv[11] #whether to rescue a read that's one-off expected cutsite (default: off)
threads = sys.argv[12] #how many threads to use

print("Indexing reads for " + pop_name)

#flags for if we're paired-end or single-end; double enzyme or single enzyme
if paired == "off":
    p_end = False
else:
    p_end = True

if cutsite_seq_two == "off":
    p_enzyme = False
else:
    p_enzyme = True
    
#if we have mixed read types, read in file and get parameters
if mixed != "off":
    try:
        mixed_file = open(mixed, "r")
        for line in mixed_file.readlines():
            parts = line.split(",")
            if parts[0] == pop_name:
                cutsite_seq = parts[1]
                if parts[2] == "-":
                    p_enzyme = False
                else:
                    p_enzyme = True
                    cutsite_seq_two = parts[2]
                if parts[3] == "yes":
                    p_end = True
                else:
                    p_end = False
                break
        mixed_file.close()
        
    except:
        print("Sample info csv [-m] is missing or wrong format")
        exit()

#get all read files and make sure they are valid
#forward
try:
    forward_file_name = dir_in + "/" + pop_name + ".F.fq"
    test_file = open(forward_file_name, "r")
    fastq_forward = SeqIO.parse(forward_file_name, 'fastq')
except IOError:
    try:
        forward_file_name = dir_in + "/" + pop_name + ".F.fastq"
        test_file = open(forward_file_name, "r")
        fastq_forward = SeqIO.parse(forward_file_name, 'fastq')
    except IOError:
        print("No input files found, make sure they are in the format SAMPLENAME.F.fq [and SAMPLENAME.R.fq, if paired-end]")
        exit()

#reverse (if paired reads)
if p_end:
    try:
        reverse_file_name = dir_in + "/" + pop_name + ".R.fq"
        test_file = open(reverse_file_name, "r")
        fastq_reverse = SeqIO.parse(reverse_file_name, 'fastq')
    except IOError:
        try:
            reverse_file_name = dir_in + "/" + pop_name + ".R.fastq"
            test_file = open(reverse_file_name, "r")
            fastq_reverse = SeqIO.parse(reverse_file_name, 'fastq')
        except IOError:
            print("No reverse read file found for " + "pop_name")

#set up data structures for forward reads
valid_sequences = set()
seqKmers_F = defaultdict(list)
nb_non_gatc_f  = 0
nb_non_gatc_f_rescued = 0
nb_total_f = 0

output_files = defaultdict(list)
hash_overlap = defaultdict(list)
combo_count = defaultdict(int)

#begin processing
print("Processing forward reads for " + pop_name)
for s in fastq_forward:
    nb_total_f += 1
    validsite = False
    cutsite_length = len(cutsite_seq)
    cutsite = s.seq[:cutsite_length]

    if rescue == True: #allow for 1bp difference
        if p_enzyme == False: #if single-enzyme
            cutsite_match_ratio = SequenceMatcher(None, cutsite_seq, cutsite).ratio()
        else: #if double-enzyme
            cutsite_match_ratio_one = SequenceMatcher(None, cutsite_seq, cutsite).ratio()
            cutsite_two_length = len(cutsite_seq_two)
            cutsite_match_ratio_two = SequenceMatcher(None, cutsite_seq_two, s.seq[:cutsite_two_length]).ratio()
            cutsite_list = [cutsite_match_ratio_one, cutsite_match_ratio_two]
            cutsite_match_ratio = max(cutsite_list)
            probable_cutsite = [cutsite_seq, cutsite_seq_two][cutsite_list.index(max(cutsite_list))]
            
        if cutsite_match_ratio >=0.75:
            validsite = True
            if cutsite_match_ratio < 1:
                nb_non_gatc_f_rescued += 1
                
    else: #not allowing for 1bp difference
        if p_enzyme == False: #if single-enzyme
            if cutsite == cutsite_seq:
                validsite = True
        else: #if double-enzyme
            if cutsite == cutsite_seq:
                validsite = True
                probable_cutsite = cutsite_seq
            elif cutsite == cutsite_seq_two:
                validsite = True
                probable_cutsite = cutsite_seq_two

    if validsite: #if cutsite is as expected
        if p_enzyme == True: #if double-enzyme
            cutsite_length = len(probable_cutsite)
        
        seq_indx = str(s.seq)[:hash_len + cutsite_length]
        seqKmers_F[seq_indx].append(s.id)
        valid_sequences.add(s.id)

        #won't worry about multi GATC here since a lot are GATCGATC or GATCAGATC
        if s.seq.count(cutsite_seq) == 2:
            gatc_idx = s.seq.find(cutsite_seq, start = cutsite_length)
            gatc_overlap = str(s.seq)[gatc_idx:(gatc_idx+hash_len)]
            if len(gatc_overlap) == hash_len and gatc_overlap != seq_indx and (len(str(s.seq)) - gatc_idx >= 20):
                combo_str = "".join([gatc_overlap, seq_indx])
                combo_count[combo_str] += 1

                #for now only interested in common combinations
                if combo_count[combo_str] >= 2:
                    hash_overlap[gatc_overlap].append(seq_indx)
                    hash_overlap[seq_indx].append(seq_indx)

                else:
                    hash_overlap[seq_indx].extend([])

        else:
            hash_overlap[seq_indx].extend([])

    else:
        nb_non_gatc_f += 1

hash_overlap.update((k, list(set(v))) for k, v in hash_overlap.items())

print "%s forward reads out of %s dropped due to ambiguous RAD-tag" % (nb_non_gatc_f, nb_total_f)
if rescue == True:
    print "%s forward reads out of %s rescued" % (nb_non_gatc_f_rescued, nb_non_gatc_f)

if p_end == True:
    #set up data structures for reverse reads
    seqKmers_R = defaultdict(list)
    nb_non_gatc_r  = 0
    nb_non_gatc_r_rescued = 0
    nb_total_r = 0

    print("Processing Reverse " + pop_name)

    for s in fastq_reverse:
        
        nb_total_r += 1
        validsite = False
        cutsite_length = len(cutsite_seq)
        cutsite = s.seq[:cutsite_length]

        if rescue==True: #allow for 1bp difference
            if p_enzyme == False: #if single-enzyme
                cutsite_match_ratio = SequenceMatcher(None, cutsite_seq, cutsite).ratio()
            else: #if double-enzyme
                cutsite_match_ratio_one = SequenceMatcher(None, cutsite_seq, cutsite).ratio()
                cutsite_two_length = len(cutsite_seq_two)
                cutsite_match_ratio_two = SequenceMatcher(None, cutsite_seq_two, s.seq[:cutsite_two_length]).ratio()
                cutsite_list = [cutsite_match_ratio_one, cutsite_match_ratio_two]
                cutsite_match_ratio = max(cutsite_list)
                probable_cutsite = [cutsite_seq, cutsite_seq_two][cutsite_list.index(max(cutsite_list))]
                
            if cutsite_match_ratio >=0.75:
                validsite = True
                if cutsite_match_ratio < 1:
                    nb_non_gatc_r_rescued += 1
        else:
            if p_enzyme == False: #if single-enzyme
                if cutsite == cutsite_seq:
                    validsite = True
            else: #if double-enzyme
                if cutsite == cutsite_seq:
                    validsite = True
                    probable_cutsite = cutsite_seq
                elif cutsite == cutsite_seq_two:
                    validsite = True
                    probable_cutsite = cutsite_seq_two

        if validsite and (s.id in valid_sequences):
            if p_enzyme == True: #if double-enzyme
                cutsite_length = len(probable_cutsite)
            seq_indx = str(s.seq)[:hash_len + cutsite_length]
            seqKmers_R[seq_indx].append(s.id)

        else:
            if keep_strand == False:
                if s.id in valid_sequences:
                    valid_sequences.remove(s.id)
                    nb_non_gatc_f += 1
            nb_non_gatc_r +=1

    print "%s reverse reads out of %s dropped due to ambiguous RAD-tag" % (nb_non_gatc_r, nb_total_r)
    if rescue == True:
        print "%s reverse reads out of %s rescued" % (nb_non_gatc_r_rescued, nb_non_gatc_r)

for kmer in seqKmers_F:
    seqKmers_F[kmer] = list(set(seqKmers_F[kmer]).intersection(valid_sequences))

seqIdToPosition = {x[1]:x[0] for x in enumerate(valid_sequences)}
positionToSeqId = {x[0]:x[1] for x in enumerate(valid_sequences)}
seqs_matrix = sp.dok_matrix((len(valid_sequences), len(valid_sequences)))

print("Processing matrix with forward reads")
for kmer in seqKmers_F:
    #seqKmers_F[kmer] = list(set(seqKmers_F[kmer]).intersection(valid_sequences))

    if len(seqKmers_F[kmer]) > 0:
         id_seq_1 =  seqIdToPosition[seqKmers_F[kmer][0]]

         for posSeq in range(1, len(seqKmers_F[kmer])):
            id_seq_2 = seqIdToPosition[seqKmers_F[kmer][posSeq]]
            seqs_matrix[id_seq_1, id_seq_2] = 1

    if len(hash_overlap[kmer]) > 0:
        for ovKmer in hash_overlap[kmer]:
            id_seq_1 =  seqIdToPosition[seqKmers_F[kmer][0]]
            if len(seqKmers_F[ovKmer]) > 0:
                id_seq_2 = seqIdToPosition[seqKmers_F[ovKmer][0]]
                seqs_matrix[id_seq_1, id_seq_2] = 1

print ("Identifying clusters")
clusters = defaultdict(list)
clusters_positions = sp.csgraph.connected_components(seqs_matrix)[1]
print("Found %s clusters" % len(set(clusters_positions)))

clusterId = -1
for c_pos in range(len(clusters_positions)):

    if clusters_positions[c_pos] != clusterId:
        cluster_id = clusters_positions[c_pos]

    clusters[cluster_id].append(positionToSeqId[c_pos])

index= 0
print("Assigning sequences to files")
for k, _ in sorted({x:len(y) for x,y in clusters.items()}.iteritems(), key=lambda (k,v): (v,k), reverse=True):
    output_files[index].append(k)
    index  = index+1
    if index >= nb_output_files:
        index = 0

for idx in output_files:
    outFileName = open("%s/%s/partial_F_%s.ids" % (dir_out, pop_name, idx), "w")
    for cluster_id in output_files[idx]:
        for seqId in clusters[cluster_id]:
            outFileName.write("%s\n" % seqId)
    outFileName.close()
