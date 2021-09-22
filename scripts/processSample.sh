#!/bin/bash

#TODO: clean this up
#shopt -s expand_aliases
#source ~/.profile
#alias pear='../programs/pear-0.9.11-linux-x86_64/bin/pear'
#alias fastq_quality_filter='../programs/fastx/fastq_quality_filter'
#alias fastq_rename='../programs/fastx/fastx_renamer'
#alias seqtk='../programs/seqtk/seqtk'
#alias vsearch='/home/mahdi/programs/vsearch/bin/vsearch'

#get parameters
_arg_in=$1
_arg_out=$2
_arg_sample=$3
_arg_index=$4
_arg_paired=$5
_arg_mixed=$6
_arg_threads=$7
_arg_quality=$8
_arg_percentbases=$9
_arg_clusterid=${10}
_arg_minlength=${11}
_arg_maxdiffs=${12}
_arg_nonambiguous=${13}

#if mixed file type, get whether paired or not
if [ "$_arg_mixed" != "off" ]
then
	_arg_paired=$(grep "$_arg_sample" $_arg_mixed | cut -d ',' -f4)
fi

#set up filenames and paths
index_name=$(echo "$_arg_index"| sed "s/.ids//" | sed "s/.*\///")
total_f_path_v1=$_arg_in/$_arg_sample".F.fastq"
total_f_path_v2=$_arg_in/$_arg_sample".F.fq"
new_f_path=$_arg_out/$_arg_sample/$index_name".F.fastq"
all_samples=$(ls "$_arg_out")
num_samples="${#all_samples[@]}"
THREADS=1

#1) subset fastq by index, allowing for .fq or .fastq
echo "Subsetting FASTQ files by index"
sed -i 's/_2:N:/_1:N:/g' ${_arg_index}
echo -e "\n"; date; echo -e "START seqtk subseq ${total_f_path_v1} ${_arg_index} > "${new_f_path}" || seqtk subseq ${total_f_path_v2} ${_arg_index} > "${new_f_path}""
../programs/seqtk/seqtk subseq ${total_f_path_v1} ${_arg_index} > "${new_f_path}" || ../programs/seqtk/seqtk subseq ${total_f_path_v2} ${_arg_index} > "${new_f_path}"
echo -e "\n"; date; echo -e "FINISH seqtk subseq ${total_f_path_v1} ${_arg_index} > "${new_f_path}" || seqtk subseq ${total_f_path_v2} ${_arg_index} > "${new_f_path}""

#2) if paired, pair reads
if [ "$_arg_paired" != "off" ]
then
  #paired
  #subset reverse reads
  total_r_path_v1=$_arg_in/$_arg_sample".R.fastq"
  total_r_path_v2=$_arg_in/$_arg_sample".R.fq"
  new_r_path=$_arg_out/$_arg_sample/$index_name".R.fastq"
  new_prefix=$_arg_out/$_arg_sample/$index_name
  
  sed -i 's/_1:N:/_2:N:/g' ${_arg_index}
  echo -e "\n"; date; echo -e "START seqtk subseq ${total_r_path_v1} ${_arg_index} > "${new_r_path}" || seqtk subseq ${total_r_path_v2} ${_arg_index} > "${new_r_path}""
  ../programs/seqtk/seqtk subseq ${total_r_path_v1} ${_arg_index} > "${new_r_path}" || ../programs/seqtk/seqtk subseq ${total_r_path_v2} ${_arg_index} > "${new_r_path}"
  echo -e "\n"; date; echo -e "FINISH seqtk subseq ${total_r_path_v1} ${_arg_index} > "${new_r_path}" || seqtk subseq ${total_r_path_v2} ${_arg_index} > "${new_r_path}""
  sed -i 's/_2:N:/_1:N:/g' ${new_r_path}
  
  #pair
  echo -e "\n"; date; echo -e "START pear -f ${new_f_path} -r ${new_r_path} -o ${new_prefix} -j ${THREADS} "
  pear -f ${new_f_path} -r ${new_r_path} -o ${new_prefix} -j ${THREADS}
  echo -e "\n"; date; echo -e "END pear -f ${new_f_path} -r ${new_r_path} -o ${new_prefix} -j ${THREADS} \n"
  pidArr=()

  #clean
  echo -e "\n"; date; echo -e "START: Quality filter"
  echo -e "\n"; date; echo -e "fastq_quality_filter -q ${_arg_quality} -p ${_arg_percentbases} -i ${new_prefix}.assembled.fastq -o ${new_prefix}.assembled.cleaned.fastq -Q33 "
  ../programs/fastx/fastq_quality_filter -q ${_arg_quality} -p ${_arg_percentbases} -i ${new_prefix}.assembled.fastq -o ${new_prefix}.assembled.cleaned.fastq -Q33  &
  pidArr+=($!)
  echo -e "\n"; date; echo -e "fastq_quality_filter -q ${_arg_quality} -p ${_arg_percentbases} -i ${new_prefix}.unassembled.forward.fastq -o ${new_prefix}.unassembled.forward.cleaned.fastq -Q 33 "
  ../programs/fastx/fastq_quality_filter -q ${_arg_quality} -p ${_arg_percentbases} -i ${new_prefix}.unassembled.forward.fastq -o ${new_prefix}.unassembled.forward.cleaned.fastq -Q33  &
  pidArr+=($!)
  echo -e "\n"; date; echo -e "fastq_quality_filter -q ${_arg_quality} -p ${_arg_percentbases} -i ${new_prefix}.unassembled.reverse.fastq -o ${new_prefix}.unassembled.reverse.cleaned.fastq -Q 33 "
  ../programs/fastx/fastq_quality_filter -q ${_arg_quality} -p ${_arg_percentbases} -i ${new_prefix}.unassembled.reverse.fastq -o ${new_prefix}.unassembled.reverse.cleaned.fastq -Q33  &
  pidArr+=($!)
  wait ${pidArr[@]}
  echo -e "\n"; date; echo -e "END: Quality filter"

  #merge
  echo -e "\n"; date; echo -e "START cat ${new_prefix}.assembled.cleaned.fastq ${new_prefix}.unassembled.forward.cleaned.fastq  ${new_prefix}.unassembled.reverse.cleaned.fastq > ${new_prefix}.fastq"
  cat ${new_prefix}.assembled.cleaned.fastq ${new_prefix}.unassembled.forward.cleaned.fastq  ${new_prefix}.unassembled.reverse.cleaned.fastq > ${new_prefix}.fastq
  echo -e "\n"; date; echo -e "END cat ${new_prefix}.assembled.cleaned.fastq ${new_prefix}.unassembled.forward.cleaned.fastq  ${new_prefix}.unassembled.reverse.cleaned.fastq > ${new_prefix}.fastq\n"
  merged_fastq="${new_prefix}.fastq"
  merged_fasta="${new_prefix}.fasta"

  python ./pop_contigs.py ${merged_fastq} ${merged_fastq}_ ${new_prefix} ${_arg_sample} 
  mv ${merged_fastq}_ ${merged_fastq}

else
  #not paired
  #clean
  echo -e "\n"; date; echo -e "START: Quality filter"
  echo -e "\n"; date; echo -e "fastq_quality_filter -q ${_arg_quality} -p ${_arg_percentbases} -i ${new_f_path} -o ${new_prefix}.assembled.cleaned.fastq -Q33 "
  ../programs/fastx/fastq_quality_filter -q ${_arg_quality} -p ${_arg_percentbases} -i ${new_prefix}.assembled.fastq -o ${new_prefix}.assembled.cleaned.fastq -Q33  &
  pidArr+=($!)
  wait ${pidArr[@]}
  echo -e "\n"; date; echo -e "END: Quality filter"
  
  #merge
  merged_fastq="${new_prefix}.fastq"
  merged_fasta="${new_prefix}.fasta"
  python ./pop_contigs.py ${new_prefix}.assembled.cleaned.fastq ${merged_fastq}_ ${new_prefix} ${_arg_sample}
  
fi

#echo "START renaming the fastq to add the population name"
#sed -i -E "s/^@([0-9]*)$/@${index_name}_\1/" ${merged_fastq}
#echo "END renaming the fastq to add the population name"

echo -e "\n"; date; echo -e "START seqtk seq -A ${merged_fastq} > ${merged_fasta}"
../programs/seqtk/seqtk seq -A ${merged_fastq} > ${merged_fasta}
echo -e "\n"; date; echo -e "END seqtk seq -A ${merged_fastq} > ${merged_fasta}\n"

echo -e "\n"; date; echo -e "START vsearch --cluster_fast ${merged_fasta} --strand plus --id ${_arg_clusterid} --msaout ${new_prefix}_1.msa --userout ${new_prefix}_1.out --userfields query+target+caln+qstrand+tstrand --mincols ${_arg_minlength} --maxdiffs ${_arg_maxdiffs} --threads ${THREADS}"
vsearch --cluster_fast ${merged_fasta} --strand plus --id ${_arg_clusterid} --msaout ${new_prefix}_1.msa --userout ${new_prefix}_1.out --userfields query+target+caln+qstrand+tstrand --mincols ${_arg_minlength} --maxdiffs ${_arg_maxdiffs} --threads ${THREADS}
echo -e "\n"; date; echo -e "END vsearch --cluster_fast ${merged_fasta} --strand plus --id ${_arg_clusterid} --msaout ${new_prefix}_1.msa --userout ${new_prefix}_1.out --userfields query+target+caln+qstrand+tstrand --mincols ${_arg_minlength} --maxdiffs ${_arg_maxdiffs} --threads ${THREADS}"

echo -e "\n"; date; echo -e "START muso.py  -t ${THREADS} -i1 ${new_prefix}_1.msa -o1 ${new_prefix}_mod_1.cons -i2 ${new_prefix}_1.out -o2 ${new_prefix}_mod_1.out -g 1 -m 3 -n 2000"
./muso.py  -t ${THREADS} -i1 ${new_prefix}_1.msa -o1 ${new_prefix}_mod_1.cons \
    -i2 ${new_prefix}_1.out -o2 ${new_prefix}_mod_1.out -g 1 --problemCtgs problematic_it_1
echo -e "\n"; date; echo -e "END muso.py  -t ${THREADS} -i1 ${new_prefix}_1.msa -o1 ${new_prefix}_mod_1.cons -i2 ${new_prefix}_1.out -o2 ${new_prefix}_mod_1.out -g 1 -m 3 -n 2000\n"

echo -e "\n"; date; echo -e "START vsearch --cluster_fast ${new_prefix}_mod_1.cons --strand both --id ${_arg_clusterid} --msaout ${new_prefix}_2.msa --userout ${new_prefix}_2.out --userfields query+target+caln+qstrand+tstrand --mincols ${_arg_minlength} --maxdiffs ${_arg_maxdiffs} --threads ${THREADS}"
vsearch --cluster_fast ${new_prefix}_mod_1.cons --strand both --id ${_arg_clusterid} --msaout ${new_prefix}_2.msa --userout ${new_prefix}_2.out --userfields query+target+caln+qstrand+tstrand --mincols ${_arg_minlength} --maxdiffs ${_arg_maxdiffs} --threads ${THREADS}
echo -e "\n"; date; echo -e "START vsearch --cluster_fast ${new_prefix}_mod_1.cons --strand both --id ${_arg_clusterid} --msaout ${new_prefix}_2.msa --userout ${new_prefix}_2.out --userfields query+target+caln+qstrand+tstrand --mincols ${_arg_minlength} --maxdiffs ${_arg_maxdiffs} --threads ${THREADS}"

echo -e "\n"; date; echo -e "START muso.py  -t ${THREADS} -i1 ${new_prefix}_2.msa -o1 ${new_prefix}.final.contigs -i2 ${new_prefix}_2.out -o2 ${new_prefix}.final.out -g 2 -m 3 -n 2000"
./muso.py  -t ${THREADS} -i1 ${new_prefix}_2.msa -o1 ${new_prefix}.final.contigs \
     -i2 ${new_prefix}_2.out -o2 ${new_prefix}.final.out -g 2 --problemCtgs problematic_it_2
echo -e "\n"; date; echo -e "END muso.py  -t ${THREADS} -i1 ${new_prefix}_2.msa -o1 ${new_prefix}.final.contigs -i2 ${new_prefix}_2.out -o2 ${new_prefix}.final.out -g 2 -m 3 -n 2000\n"

echo -e "\n"; date; echo -e "START trackOverlaps.py -i1 ${new_prefix}_mod_1.out  -i2 ${new_prefix}.final.out  -o ${new_prefix}.mapping_to_cons"
./x.py ${new_prefix}_mod_1.out   ${new_prefix}.final.out > ${new_prefix}.mapping_to_cons
echo -e "\n"; date; echo -e "END trackOverlaps.py -i1 ${new_prefix}_mod_1.out  -i2 ${new_prefix}.final.out  -o ${new_prefix}.mapping_to_cons\n"

####OLD BELOW

# tmpPath=$2
# SAMPLE_NAME=$1
# tmpName=$(echo "$tmpPath"| sed "s/.ids//" | sed "s/.*\///")
# addF=".F.fq"; addR=".R.fq"; newF="$tmpName$addF"; newR="$tmpName$addR"
# addNew=".tmp.ids"; newIDs="$tmpName$addNew"
# 
# #echo "$tmpPath"
# #echo "$SAMPLE_NAME"
# #echo "$tmpName"
# 
# cd "/10tb_abyss/emily/seanome_project/erika_seanome/erika_runs_12_19/"$SAMPLE_NAME
# sed 's/$/ /' $tmpPath > $newIDs
# 
# #switched from A1 to A3
# grep --color=auto -A1 --no-group-separator -F -f $newIDs "/10tb_abyss/emily/seanome_project/erika_seanome/erika_raw/adapters_removed_only/"$SAMPLE_NAME".F.fq" > $newF
# grep --color=auto -A1 --no-group-separator -F -f $newIDs "/10tb_abyss/emily/seanome_project/erika_seanome/erika_raw/adapters_removed_only/"$SAMPLE_NAME".R.fq" > $newR
# 
# input_forward_path=$newF
# input_reverse_path=$newR
# NAME=$tmpName
# THREADS=40
# 
# echo -e "\n"; date; echo -e "START pear -f ${input_forward_path} -r ${input_reverse_path} -o ${NAME} -j ${THREADS} " 
# pear -f ${input_forward_path} -r ${input_reverse_path} -o ${NAME} -j ${THREADS} 
# echo -e "\n"; date; echo -e "END pear -f ${input_forward_path} -r ${input_reverse_path} -o ${NAME} -j ${THREADS} \n" 
# pidArr=()
# 
# echo -e "\n"; date; echo -e "START: Quality filter" 
# echo -e "\n"; date; echo -e "fastq_quality_filter -q 20 -p 75 -i ${NAME}.assembled.fastq -o ${NAME}.assembled.cleaned.fastq -Q 33 " 
# fastq_quality_filter -q 20 -p 75 -i ${NAME}.assembled.fastq -o ${NAME}.assembled.cleaned.fastq -Q 33  &
# pidArr+=($!)
# echo -e "\n"; date; echo -e "fastq_quality_filter -q 20 -p 75 -i ${NAME}.unassembled.forward.fastq -o ${NAME}.unassembled.forward.cleaned.fastq -Q 33 " 
# fastq_quality_filter -q 20 -p 75 -i ${NAME}.unassembled.forward.fastq -o ${NAME}.unassembled.forward.cleaned.fastq -Q 33  &
# pidArr+=($!)
# echo -e "\n"; date; echo -e "fastq_quality_filter -q 20 -p 75 -i ${NAME}.unassembled.reverse.fastq -o ${NAME}.unassembled.reverse.cleaned.fastq -Q 33 " 
# fastq_quality_filter -q 20 -p 75 -i ${NAME}.unassembled.reverse.fastq -o ${NAME}.unassembled.reverse.cleaned.fastq -Q 33  &
# pidArr+=($!)
# wait ${pidArr[@]}
# echo -e "\n"; date; echo -e "END: Quality filter"\n
# 
# echo -e "\n"; date; echo -e "START cat ${NAME}.assembled.cleaned.fastq ${NAME}.unassembled.forward.cleaned.fastq  ${NAME}.unassembled.reverse.cleaned.fastq > ${NAME}.fastq" 
# cat ${NAME}.assembled.cleaned.fastq ${NAME}.unassembled.forward.cleaned.fastq  ${NAME}.unassembled.reverse.cleaned.fastq > ${NAME}.fastq
# echo -e "\n"; date; echo -e "END cat ${NAME}.assembled.cleaned.fastq ${NAME}.unassembled.forward.cleaned.fastq  ${NAME}.unassembled.reverse.cleaned.fastq > ${NAME}.fastq\n" 
# merged_fastq="${NAME}.fastq" 
# merged_fasta="${NAME}.fasta" 
# 
# python /10tb_abyss/emily/seanome_project/pop_contigs.py ${merged_fastq} ${merged_fastq}_ ${NAME}
# mv ${merged_fastq}_ ${merged_fastq}
# 
# echo "START renaming the fastq to add the population name"
# sed -i -r "s/^@([0-9]*)$/@${NAME}_\1/" ${merged_fastq}
# echo "END renaming the fastq to add the population name"
# 
# echo -e "\n"; date; echo -e "START seqtk seq -A ${merged_fastq} > ${merged_fasta}" 
# seqtk seq -A ${merged_fastq} > ${merged_fasta}
# echo -e "\n"; date; echo -e "END seqtk seq -A ${merged_fastq} > ${merged_fasta}\n" 
# 
# echo -e "\n"; date; echo -e "START vsearch --cluster_fast iter1"
# vsearch --cluster_fast ${merged_fasta} --strand plus --id 0.95 --msaout ${NAME}_1.msa --userout ${NAME}_1.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}
# echo -e "\n"; date; echo -e "END vsearch --cluster_fast iter 1"
# 
# echo -e "\n"; date; echo -e "START muso.py  -t ${THREADS} -i1 ${NAME}_1.msa -o1 ${NAME}_mod_1.cons -i2 ${NAME}_1.out -o2 ${NAME}_mod_1.out -g 1 -m 3 -n 2000" 
# /10tb_abyss/emily/seanome_project/programs/Seanome_BU/scripts/muso.py  -t ${THREADS} -i1 ${NAME}_1.msa -o1 ${NAME}_mod_1.cons \
#     -i2 ${NAME}_1.out -o2 ${NAME}_mod_1.out -g 1 \
#     --problemCtgs problematic_it_1
# echo -e "\n"; date; echo -e "END muso.py  -t ${THREADS} -i1 ${NAME}_1.msa -o1 ${NAME}_mod_1.cons -i2 ${NAME}_1.out -o2 ${NAME}_mod_1.out -g 1 -m 3 -n 2000\n" 
# 
# echo -e "\n"; date; echo -e "START vsearch --cluster_fast iter2"
# vsearch --cluster_fast ${NAME}_mod_1.cons --strand both --id 0.95 --msaout ${NAME}_2.msa --userout ${NAME}_2.out --userfields query+target+caln+qstrand+tstrand --mincols 80 --maxdiffs 10 --threads ${THREADS}
# echo -e "\n"; date; echo -e "START vsearch --cluster_fast iter2"
# 
# echo -e "\n"; date; echo -e "START muso.py  -t ${THREADS} -i1 ${NAME}_2.msa -o1 ${NAME}.final.contigs -i2 ${NAME}_2.out -o2 ${NAME}.final.out -g 2 -m 3 -n 2000"
#  /10tb_abyss/emily/seanome_project/programs/Seanome_BU/scripts/muso.py  -t ${THREADS} -i1 ${NAME}_2.msa -o1 ${NAME}.final.contigs \
#      -i2 ${NAME}_2.out -o2 ${NAME}.final.out -g 2 --problemCtgs problematic_it_2 
# echo -e "\n"; date; echo -e "END muso.py  -t ${THREADS} -i1 ${NAME}_2.msa -o1 ${NAME}.final.contigs -i2 ${NAME}_2.out -o2 ${NAME}.final.out -g 2 -m 3 -n 2000\n" 
# 
# echo -e "\n"; date; echo -e "START trackOverlaps.py -i1 ${NAME}_mod_1.out  -i2 ${NAME}.final.out  -o ${NAME}.mapping_to_cons" 
# /10tb_abyss/emily/seanome_project/programs/Seanome_BU/scripts/x.py ${NAME}_mod_1.out   ${NAME}.final.out > ${NAME}.mapping_to_cons
# echo -e "\n"; date; echo -e "END trackOverlaps.py -i1 ${NAME}_mod_1.out  -i2 ${NAME}.final.out  -o ${NAME}.mapping_to_cons\n" 
# 
# 
# 
