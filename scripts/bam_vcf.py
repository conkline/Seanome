def writeBAM(bamFile, inList):
  
  with bamFile as outf:
  
    #write BAM entry for each read
    for j in range(len(inList)):
      locusID = inList[j][0]
      tmp_read_dict = inList[j][1]
      final_cigars = inList[j][2]
      final_readnames = inList[j][3]
      final_contig = inList[j][4]
      
      for i in range(len(final_readnames)):
      
        seq = str(tmp_read_dict[final_readnames[i]].seq)
        #wild that seqrecord doesn't have a quality string attribute but c'est la vie
        quality = tmp_read_dict[final_readnames[i]].format("fastq").split("\n")[-2]
        pop = final_readnames[i].split("_F_")[0]
        #convert cigar to bam format
        cigar_n = final_cigars[i].replace("I", "N")
        
        a = pysam.AlignedSegment()
        a.query_name = final_readnames[i]
        a.query_sequence = seq
        a.flag = 1
        a.reference_id = locusID #note - must be an integer
        a.reference_start = 1 #may need to change this in the future
        a.mapping_quality = 255 #TODO: calculate quality
        a.cigarstring = cigar_n #note - have to change 'I' to 'N'
        a.template_length = len(final_contig)
        a.query_qualities = pysam.qualitystring_to_array(quality)
        a.tags = ([("RG", pop)])
        outf.write(a)
      
  return
    
def dictToCounts(popdict, pop, ref, alts):
  #in format popname: [A, C, G, T]
  #need format [ref, alt, alt2, etc.]
  transDict = {"A":0, "C":1, "G":2, "T":3}
  refOut = [popdict[pop][transDict[ref]]]
  altOut = [popdict[pop][transDict[x]] for x in alts]
  
  #generate haplotype for this sample - ref is 0, first alt is 1, etc.
  haploList = [ref] + alts
  transDict2 = {0:"A", 1:"C", 2:"G", 3:"T"}
  haploDict = {haploList[x]:x for x in range(len(haploList))}
  haplos = []
  
  for nuc in range(4):
    if popdict[pop][nuc] != 0:
      nuc_letter = transDict2[nuc]
      haplos.append(haploDict[nuc_letter])
  
  return(refOut + altOut, sorted(haplos))
    
def writeVCF(vcfFile, vcf_header, inList):
  
  with vcfFile as outf:
    for i in range(len(inList)):
      locus = inList[i][0]
      snps = inList[i][1]
      
      for snp in snps:
      
        #TODO - calculate quality (QUAL and GQ)
        #create new record
        r = outf.new_record(contig=snp["chrom"], start=snp["pos"]-1, stop=snp["pos"], 
                            alleles=[snp["ref"]] + snp["alt"], id="_".join([snp["chrom"], str(snp["pos"])]))
        
        #add info                    
        r.info["NS"] = snp["ns"]
        r.info["DP"] = snp["depth"]
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
        
        outf.write(r)
        
  return

# def writeBAM(locusID, bamFile, tmp_read_dict, final_cigars, final_readnames, final_contig):
#   
#   with bamFile as outf:
#   
#     #write BAM entry for each read
#     for i in range(len(final_readnames)):
#       
#       seq = str(tmp_read_dict[final_readnames[i]].seq)
#       #wild that seqrecord doesn't have a quality string attribute but c'est la vie
#       quality = tmp_read_dict[final_readnames[i]].format("fastq").split("\n")[-2]
#       pop = final_readnames[i].split("_F_")[0]
#       #convert cigar to bam format
#       cigar_n = final_cigars[i].replace("I", "N")
#       
#       a = pysam.AlignedSegment()
#       a.query_name = final_readnames[i]
#       a.query_sequence = seq
#       a.flag = 1
#       a.reference_id = locusID #note - must be an integer
#       a.reference_start = 1 #may need to change this in the future
#       a.mapping_quality = 255 #TODO: calculate quality
#       a.cigarstring = cigar_n #note - have to change 'I' to 'N'
#       a.template_length = len(final_contig)
#       a.query_qualities = pysam.qualitystring_to_array(quality)
#       a.tags = ([("RG", pop)])
#       outf.write(a)
#       
#   return
#     
# def dictToCounts(popdict, pop, ref, alts):
#   #in format popname: [A, C, G, T]
#   #need format [ref, alt, alt2, etc.]
#   transDict = {"A":0, "C":1, "G":2, "T":3}
#   refOut = [popdict[pop][transDict[ref]]]
#   altOut = [popdict[pop][transDict[x]] for x in alts]
#   
#   #generate haplotype for this sample - ref is 0, first alt is 1, etc.
#   haploList = [ref] + alts
#   transDict2 = {0:"A", 1:"C", 2:"G", 3:"T"}
#   haploDict = {haploList[x]:x for x in range(len(haploList))}
#   haplos = []
#   
#   for nuc in range(4):
#     if popdict[pop][nuc] != 0:
#       nuc_letter = transDict2[nuc]
#       haplos.append(haploDict[nuc_letter])
#   
#   return(refOut + altOut, sorted(haplos))
#     
# def writeVCF(locus, vcfFile, vcf_header, snps):
#   
#   with vcfFile as outf:
#     
#     for snp in snps:
#     
#       #TODO - calculate quality (QUAL and GQ)
#       #create new record
#       vcf_header.add_meta('contig', items=[('ID', snp["chrom"])])
#       r = vcf_header.new_record(contig=snp["chrom"], start=snp["pos"]-1, stop=snp["pos"], 
#                           alleles=[snp["ref"]] + snp["alt"], id="_".join([snp["chrom"], str(snp["pos"])]))
#       
#       #add info                    
#       r.info["NS"] = snp["ns"]
#       r.info["DP"] = snp["depth"]
#       r.info["AC"] = ",".join([str(x) for x in snp["counts"]])
#       r.info["TYPE"] = snp["type"]
#       
#       #add sample info
#       for pop in snp["pops"]:
#         popcounts, haplos = dictToCounts(snp["pop_counts"], pop, snp["ref"], snp["alt"])
#         bytes_pop = pop.encode('utf-8')
#         r.samples[bytes_pop]["DP"] = snp["pop_dp"][pop]
#         r.samples[bytes_pop]["AD"] = ",".join([str(x) for x in popcounts])
#         r.samples[bytes_pop]["GT"] = "/".join([str(x) for x in haplos])
#         r.samples[bytes_pop]["GQ"] = "." 
#       
#       outf.write(r)
#       
#   return
