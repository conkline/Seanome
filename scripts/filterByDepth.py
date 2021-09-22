from collections import defaultdict
import sys
import os

os.chdir("/10tb_abyss/emily/snail_seanome/2021_paolo_seanome_indexed/ALL")
mindepth = int(sys.argv[1])

old_mapping = open("ALL.mapping", "r")
new_mapping = open("ALL.filtered.mapping", "w")
new_ids_file = open("ALL.filtered.ids", "w")

mapping_dict = defaultdict(int)
mapping_line = defaultdict(list)

#read in old mapping
for line in old_mapping.readlines():
  parts = line.split()
  try:
    contig = parts[1]
    mapping_dict[contig] += 1
    mapping_line[contig].append(line)
  except IndexError:
    continue

#filter and write new files
new_ids = [x for x in mapping_dict.keys() if mapping_dict[x] >= mindepth]

print("Finished filtering, found " + str(len(new_ids)) + " valid contigs")

for line in new_ids:
  for x in mapping_line[line]:
    new_mapping.write(x)
  new_ids_file.write(line + "\n")

old_mapping.close()
new_mapping.close()
new_ids_file.close()
