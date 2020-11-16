# -*- coding: utf-8 -*-
"""
Updated on Sat Nov 14 2020

@author: irina_chelysheva
"""
####################################################################
#Title:
#RNA2HLA - checking the identity of RNA-seq samples based on HLA types
#Synopsys:
#Script extracts the HLA types of I an II classes from all the files in the folder containing raw RNA-seq data (paired on single-end), 
#The alleles are then cross-compared between the samples to verify the common sourse of the samples based on HLA identities
#
#Release: 1.0
#
#Author: Irina Chelysheva, 2019-2020
#Oxford Vaccine Group, University of Oxford
#Contact: irina.chelysheva@paediatrics.ox.ac.uk
#
#Usage: 
#python RNA2HLA.py -f /folder_with_raw_RNAseq_data [-r /global_name_of_run]* [-p <int>]** [-3 <int>]*** [-c <float>]**** [-g <int]*****
#*optional: to be used as a prefix for all output files
#**optional: number of parallel search threads for bowtie optional (default: 6)
#***optional: trim int bases from the low-quality end of each read
#****optional: confidence level for HLA-typing (default: 0.05)
#*****optional: number of HLA genes to be included for typing (default: 5, maybe increased to 6 - adding DQB1)
#folder should contain raw RNA-seq samples, single- or paired-end or both types, in a compressed or not compressed formats
#
#Dependencies:
#1) RNA2HLA is a python script (python 2).
#2) Dependent python scripts single_end.py and paired_end.py must be located in the same folder.
#3) bowtie must be reachable by the command "bowtie".
#4) R must be installed.
#5) Index files must be located in the folder "references".
#6) Packages: biopython (developed with V1.58), numpy (1.3.0).
#
#Output:
#The final output - overall comparison matrix in csv format, which crosscompares all RNA-seq samples in the folder.
#Individual outputs in txt format produced for each RNA-seq sample in the folder (classes I and II are written in one file): 
#1)<sampleID>.bowtielog.txt - file with statistics of HLA mapping; 
#2) <sampleID>.ambiguity.txt - reports typing ambuigities (if more than one solution for an allele possible);
#3) <sampleID>.expression.txt - RPKM expression of HLA;
#4) <sampleID>.HLAgenotype4digits.txt - 4 digital HLA type.

####################################################################

from optparse import OptionParser
import os
import re
import sys
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns

version="1.0"

if __name__ == '__main__':
	parser = OptionParser(usage="usage: %prog -f readFolder [-r globalrunName] [-p <int>] [-3 <int>] [-c <float>] [-g <int>]", version="%prog 1")
	parser.add_option("-f",
			action="store", 
			dest="readFolder",
			help="Folder with RNA-seq data (uncompressed or gzipped fastq files)")
	parser.add_option("-r", "--globalrunName",
			action="store", 
			dest="globalrunName",
			default="HLA",             
			help="Name of this HLA typing run - the output subfolder, which will be created and where all the results will be written. Default: HLA")
	parser.add_option("-p", "--threads",
                        action="store",
                        dest="threads",
                        default="6",
                        help="Bowtie option: Launch <int> parallel search threads. Default (RNA2HLA): 6")
	parser.add_option("-3", "--trim3",
			action="store",
			dest="trim3",
			default="0",
			help="Bowtie option: -3 <int> trims <int> bases from the low quality 3' end of each read. Default: 0")
	parser.add_option("-c", "--conf",
			action="store",
			dest="conf",
			default="0.05",
			help="Confidence threshold for defining the 4digits allele of each HLA gene. Default: p<0.05")
	parser.add_option("-g", "--genes",
			action="store",
			dest="genes",
			default="5",
			help="Number of HLA genes included for typing. Default: 5, maybe increased to 6 (adding DQB1)")

	(options, args) = parser.parse_args()
	if not options.readFolder:   
		parser.error('Folder with RNA-seq data (uncompressed or gzipped fastq files) is not given.')
	readFolder=options.readFolder
	globalrunName=options.globalrunName
	trim3=str(options.trim3)
	threads=str(options.threads)
	conf=float(options.conf)
	if conf>=1:
		parser.error("Confidence level (p) has to be less than 1!")
	counter=0
	genes=int(options.genes)
	if genes!=5:
         if genes!=6:
             parser.error("Number of HLA genes to include has to be 5 or 6!") 
	genes=str(options.genes)
	if not os.path.exists(options.readFolder+"/"+options.globalrunName):
         os.mkdir(options.readFolder+"/"+options.globalrunName)
             
##1) Identifying HLAtypes for all RNA-seq samples in the directory 
for file in os.listdir(readFolder):
#first search for all the fastq files in the directory - compressed or not compressed
    if file.endswith(".fastq")|file.endswith(".fastq.gz")|file.endswith("fq")|file.endswith("fq.gz"):
        counter=counter+1 #counting the number of RNA-seq files in the folder
        #setting the variables to check sample type (single- or paired-end)
        readFile1=None
        readFile2=None
        testreadFile1=None
        testreadFile2=None
        #samplename will be used as prefix for the output files
        samplename=None
        #checking if the file is the 1st file of pair-end sample
        pair=re.search(r"._1.f|-1.f|\.1.f",os.path.join(file))
        if pair:
            readFile1=os.path.join(file)
            samplename=readFile1.partition("_1.f")[0]
            if samplename==readFile1:
                samplename=readFile1.partition("-1.f")[0]
            if samplename==readFile1:
                samplename=readFile1.partition(".1.f")[0]
            #looking for the 2nd file from paired-end RNAseq in the given folder
            for file2 in os.listdir(readFolder):
                if file2.startswith(samplename+"_2.f"):
                    readFile2=os.path.join(file2)
            if not readFile2:
                for file2 in os.listdir(readFolder):
                    if file2.startswith(samplename+"-2.f"):
                        readFile2=os.path.join(file2)
            if not readFile2:
                for file2 in os.listdir(readFolder):
                    if file2.startswith(samplename+".2.f"):
                        readFile2=os.path.join(file2)
            if readFile2:
                print "Start typing HLA from "+samplename+" - paired-end RNA-seq sample"
                #####send to paired-end RNA-seq script#####
                if globalrunName:
                    runName=(globalrunName+"/"+samplename)
                else:
                    runName=samplename
    
                os.system("python paired_end.py -1 "+readFile1+" -2 "+readFile2+" -f "+readFolder+" -r "+runName+" -p "+threads+" -3 "+trim3+ " -g "+genes)
                
            if not readFile2:
                parser.error('Input RNA-seq sample '+samplename+' seems to be paired-end, but file #2 has not been found.')
        else:
            #checking if the given file is the 2nd file of paired-end sample
            pair2=re.search(r"._2.f|-2.f|\.2.f",os.path.join(file))
            if pair2:
                #now checking the existence of 1st file in pair
                testreadFile2=os.path.join(file)
                testsamplename=testreadFile2.partition("_2.f")[0]
                if testsamplename==testreadFile2:
                    testsamplename=testreadFile2.partition(".2.f")[0]
                if testsamplename==testreadFile2:
                    testsamplename=testreadFile2.partition("-2.f")[0]               
                for testfile1 in os.listdir(readFolder):
                    if testfile1.startswith(testsamplename+"_1.f"):
                        testreadFile1=os.path.join(testfile1)
                if not testreadFile1:
                 for testfile1 in os.listdir(readFolder):
                     if testfile1.startswith(testsamplename+".1.f"):
                        testreadFile1=os.path.join(testfile1)
                if not testreadFile1:
                 for testfile1 in os.listdir(readFolder):
                     if testfile1.startswith(testsamplename+"-1.f"):
                        testreadFile1=os.path.join(testfile1)        
                if not testreadFile1:
                    parser.error('Input RNA-seq sample '+testsamplename+' seems to be paired-end, but file #1 has not been found.')                     
            if not pair2:
                #if not paired - input is single-end
                readFile1=os.path.join(file)
                #extract the samplename to use as prefix
                samplename=readFile1.partition(".fastq")[0]
                if samplename==readFile1:
                     samplename=readFile1.partition(".fq")[0]
                print "Start typing HLA from "+samplename+" - single-end RNA-seq sample"
                #####send to single-end RNA-seq script#####
                if globalrunName:
                    runName=(globalrunName+"/"+samplename)
                else:
                    runName=samplename
                os.system("python single_end.py -1 "+readFile1+" -f "+readFolder+" -r "+runName+" -p "+threads+" -3 "+trim3+ " -g "+genes)
if counter==0:
    parser.error('Given folder does not contain RNA-seq files!')


##2) When all the HLAtypes are identified for all the samples - we cross compare the samples identities pairwise
samplelist=[]
for file in os.listdir(readFolder+"/"+globalrunName):
    samplename=None
    if file.endswith(".HLAgenotype4digits.txt"):
         samplename=os.path.join(file)
         samplename=samplename.partition(".HLAgenotype4digits.txt")[0]
         samplelist.append(samplename)
df = pd.DataFrame(index=samplelist, columns=samplelist)
for file in samplelist:
    #run over each possible couple of samples 1 and 2 within the folder
    samplename1=file
    f1=open(readFolder+"/"+globalrunName+"/"+file+".HLAgenotype4digits.txt")
    #read the lines in output file containing genotypes for the 1 sample
    lines1=f1.readlines()
    for file2 in os.listdir(readFolder+"/"+globalrunName):
        if file2.endswith(".HLAgenotype4digits.txt"):
            samplename2=os.path.join(file2)
            samplename2=samplename2.partition(".HLAgenotype4digits.txt")[0]
            f2=open(readFolder+"/"+globalrunName+"/"+file2)
            lines2=f2.readlines()
            #define the lines containing the HLA-types for each HLA gene
            if genes=="5":
                mylines=[2, 3, 4, 7, 8]
            else:
                mylines=[2, 3, 4, 7, 8, 9]
            #set up the counter to count the identical alleles between two samples
            counter=0 
            #set up the counter to count total number of alleles (out of 12) considered for comparison in each couple of samples
            allalleles=0
            for i in mylines:
                #from now on: sample 1 index with x, sample 2 indexed with y
                Ax1=None
                Ay1=None
                Ax2=None
                Ay2=None
                Ax11=None
                Ax22=None
                Ay11=None
                Ay22=None
                allelesx=2
                allelesy=2
                Ax=lines1[i]
                #cut only the value of significance for each determined allele
                Ax11=Ax.split("\t")[2]  
                #compare the p-value to our confidence level (conf) - threshold
                try:
                    if float(Ax11)<=conf:
                        #Ax1=Ax.split("*",1)[1][0:5]
                        Ax1=Ax.split("\t")[1]
                #if p is more than the confidence level we exclude this allele
                    else:
                        Ax1="x1"
                        allelesx=allelesx-1
                #if p is NA, this part will be running and the allele will be excluded
                except ValueError:
                    pass
                    Ax1="x1"
                    allelesx=allelesx-1
                Ax22=Ax.split("\t")[4]
                try:
                    if float(Ax22)<=conf:
                       #Ax2=Ax.rsplit("*",1)[1][0:5]
                        Ax2=Ax.split("\t")[3]
                    else:
                        Ax2="x2"
                        allelesx=allelesx-1
                except ValueError:
                   pass
                   Ax2="x2"
                   allelesx=allelesx-1
                Ay=lines2[i]
                Ay11=Ay.split("\t")[2]
                try:
                    if float(Ay11)<=conf:
                        #Ay1=Ay.split("*",1)[1][0:5]
                        Ay1=Ay.split("\t")[1]
                    else:
                        Ay1="y1"
                        allelesy=allelesy-1
                except ValueError:
                   pass
                   Ay1="y1"
                   allelesy=allelesy-1                   
                Ay22=Ay.split("\t")[4] 
                try:
                    if float(Ay22)<=conf:
                        #Ay2=Ay.rsplit("*",1)[1][0:5]
                        Ay2=Ay.split("\t")[3]
                    else:
                        Ay2="y2"
                        allelesy=allelesy-1
                except ValueError:
                    pass
                    Ay2="y2"
                    allelesy=allelesy-1
                #print Ay2
                if int(allelesx)>0 and int(allelesy)>0:
                    alleles=(allelesx+allelesy)//2
                else:
                    alleles=0
                allalleles=allalleles+alleles
                #print allalleles
                if alleles==2:
                    if Ax1==Ay1 or Ax1==Ay2:
                        counter=counter+1
                    if Ax2==Ay2 or Ax2==Ay1:
                        counter=counter+1
                if alleles==1:
                    if Ax1==Ay1 or Ax1==Ay2:
                        counter=counter+1
                    elif Ax2==Ay2 or Ax2==Ay1:
                        counter=counter+1
            df.at[samplename1,samplename2]=round(float(counter)/float(allalleles)*100,2),counter,allalleles

print df
#Finally writing a comparison matrix of all the samples and their % of identical HLA genotypes
df.to_csv(readFolder+"/"+globalrunName+"/all_samples_comparison_matrix.csv", sep='\t')
