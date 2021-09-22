from Bio import SeqIO
import glob
import sys

filein = sys.argv[1]
fileout = sys.argv[2]
indexname = sys.argv[3].split("/")[-1]
popname = sys.argv[4]
counter = 1

with open(filein, "r") as old, open(fileout, "w") as new:
    records = SeqIO.parse(old, 'fastq')
    for r in records:
        r.id = (indexname + "_" + str(counter)).replace("partial", popname)
        r.description = ""
        counter += 1
        SeqIO.write(r, new, 'fastq')

#popname = "ACU"
#all_files = glob.glob("sample_runs_2/" + popname + "/*final.contigs")
#new_file = "sample_runs_2/" + popname + "/" + popname + "_contigs.fasta"
#counter = 1

#with open(new_file, "w") as merged:
#    for i in all_files:
#        with open(i, "r") as old:
#            records = SeqIO.parse(old, 'fasta')
#            for r in records:
#	        r.id = popname + "_" + str(counter)
#                SeqIO.write(r, merged, 'fasta')
#		counter += 1
