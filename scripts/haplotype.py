from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import re

def callHaplos(locus, alignment, read_names, snps_pos, minCount):

  #call haplotypes
  haplo_counts = dict()
  
  #build haplotypes
  haplo_dict = {}
  for row in range(len(alignment)):
    haplo_dict[read_names[row]] = tuple(alignment.iloc[row, snps_pos].tolist())
    
  #finally, count haplotypes and filter by MAF
  valid_haplos_test = [x for x in set(haplo_dict.values()) if list(haplo_dict.values()).count(x) >= minCount]
  
  #deal with NAs from overlapping reads (just for counting purposes)
  #if haplotype with NA is otherwise unique, keep it in the "valid haplos"
  #otherwise, lump it in with others and discard from list
  if any([pd.np.nan in x for x in valid_haplos_test]):
    valid_haplos = [x for x in valid_haplos_test if pd.np.nan not in x]
    has_na = [x for x in valid_haplos_test if pd.np.nan in x]
    
    for na_haplo in has_na: #for each haplotype that has an NA
      is_unique = False
      for bp_idx in range(len(na_haplo)): #for each bp in this haplotype
        if not pd.isna(na_haplo[bp_idx]): #if we don't have an NA at this location
          #check against other base pairs at this location to see if we're unique
          test_haplo = na_haplo[bp_idx]
          test_bp = [x[bp_idx] for x in valid_haplos if not pd.isna(x[bp_idx])]
          if test_haplo not in test_bp:
            is_unique = True
      if is_unique:
        valid_haplos.append(na_haplo)
    
  else:
    #if no NAs, keep everything
    valid_haplos = valid_haplos_test
  
  #link haplotypes back to sample
  sample_haplos = defaultdict(set)
  for x in haplo_dict.keys():
    tmp_sample = x.split("_")[0]
    tmp_haplo = haplo_dict[x]
    
    if tmp_haplo in valid_haplos:
      sample_haplos[tmp_sample].add(tmp_haplo)
      
  #record max # haplotypes per sample at this locus
  if len(sample_haplos) != 0:
    haplo_counts[locus] = max([len(x) for x in sample_haplos.values()])
  else:
    haplo_counts[locus] = 0
    
  return(haplo_counts)
        
