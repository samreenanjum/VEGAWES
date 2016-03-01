#!/bin/bash

#################################
## Script to generate GC files for chromosomes given the fasta file and targets list
## Required: samtools
## Input: Fasta file, Targets list, Name of the GCC folder
## Project: VEGAWES
## Author: SA
##
## command ---> ./GCScript.sh fastafile targetsfile GCCfolder
################################

## Set the parameters
fastafile=$1    ## Reference Fasta file
targetslist=$2  ## targets interval list
GCCfolder=$3    ## folder name to save the gc results

## List of chromosomes from the targets file
chrlist=`awk  '{split($0,a,":"); print a[1]}'  $targetslist | uniq`
#echo $chrlist

## For each chromosome
for chr in $chrlist
do
echo $chr
#Create the output file
echo -e "chr\tprobe_Start\tprobe_end\tGCContent" > $GCCfolder/GCContent.$chr.txt  
#Create a temp targets list for the chr
awk -v find=$chr '{split($0,a,":"); if(a[1]==find) print $0}'  $targetslist  > $GCCfolder/chr.targets.txt

  while read target
	do
    echo $target | awk  '{split($0,a,"[:|-]"); printf("%d\t%d\t%d\t",  a[1], a[2], a[3])}' >> $GCCfolder/GCContent.$chr.txt
    samtools faidx $fastafile $target | sed '1d' | grep -o -i "[G|C]"  | wc -l >> $GCCfolder/GCContent.$chr.txt
  done < $GCCfolder/chr.targets.txt 
  
##Optional: sort the file  
#sort -un -k2 GCContent.$chr.txt > sorted.gc.$chr.txt
done

# Remove the temp targets file
rm $GCCfolder/chr.targets.txt
