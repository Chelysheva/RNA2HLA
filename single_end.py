# -*- coding: utf-8 -*-
"""
Updated on Sat Nov 14 2020

@author: irina_chelysheva
"""

##########################################################################################################

#Title:
#Single_end - dependent script from RNA2HLA tool.
#
#Synopsis:
#Script for HLA typing from single-end RNA-seq sample.
#
#Author: Irina Chelysheva, 2019-2020
#Oxford Vaccine Group, Department of Paediatrics, University of Oxford
#Contact: irina.chelysheva@paediatrics.ox.ac.uk
#
#References:
#Developed based on from seq2HLA (release: 2.3) by Sebastian Boegel, 2012 - 2017 (c)
#BDOI: 10.1186/gm403.
#Usage: 
#Script is automatically called while running RNA2HLA and is not intended to be used independently.

##########################################################################################################

from operator import itemgetter
import sys
import linecache
import os
from Bio import SeqIO
import numpy
from optparse import OptionParser
import fourdigits
import gzip
import subprocess
from glob import glob
import fileinput

version="1.0"

def main_single(runName,readFile1,fastaClassI,fastaClassII,bowtiebuildClassI,bowtiebuildClassII,trim3,threads,genes,class1_list,class2_list,length_dict1,length_dict2):
	#print the sample and run name
	#print "Running RNA2HLA based on seq2HLA on a run "+runName
	#define a type of input RNA-seq sample
	#print "Input is a single-end RNA-seq sample...."
	#determine if input is uncompressed or gzipped fastq files and determine number of lines (important for expression calculation
	f=gzip.open(readFile1, 'r')
	#if first line can be read without error, input is gzipped
	try:
		first_line = f.readline()
                gzipped=1
		#process_wc=subprocess.Popen(['bash','-c','zcat '+readFile1+' | sed \'2q;d\' - | wc -L'],stdout=subprocess.PIPE)		
		cmd = "zcat %s | sed '2q;d' | wc -L" % readFile1
		process_wc=subprocess.Popen( ['bash','-c',cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		readlength=process_wc.communicate()[0]
		print "Input is a gzipped file ...."
        except Exception, e:
		cmd = "sed '2q;d' %s | wc -L" % readFile1
                process_wc=subprocess.Popen( ['bash','-c',cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                readlength=process_wc.communicate()[0]
                gzipped=0
		print "Input is uncompressed fastq file ...."
	#as shown in the publication HLA typing with RNA-Seq works best by allowing as less mismatches as necessary
	if int(readlength)<=50:
                mismatch=1
        else:
                if int(readlength)<=100:
                        mismatch=2
                else:
                        mismatch=3
	#concatenate mapping parameters
	mapopt="-p "+str(threads)+" -a -v"+str(mismatch)
	print "The read length of your input fastq was determined to be "+str(int(readlength))+", so "+str(mismatch)+" mismatches will be allowed and "+str(threads)+" threads will be used by bowtie."	
	#call HLA typing for Class I
	print str(numpy.datetime64('now'))
	mainHLA(runName+"-ClassI",readFile1,bowtiebuildClassI,fastaClassI,mapopt,trim3,genes,gzipped,class1_list,1,length_dict1)
	#call HLA typing for Class II
	print str(numpy.datetime64('now'))
	mainHLA(runName+"-ClassII",readFile1,bowtiebuildClassII,fastaClassII,mapopt,trim3,genes,gzipped,class2_list,2,length_dict2)	
	merge(runName)
#---------------Class I-------------------------------
def mainHLA(runName,readFile1,bowtiebuild,hla1fasta,mapopt,trim3,genes,gzip,locus_list,hla_class,length_dict):
	twodigits1={}
	fourdigits1={}
	fourDigit_solutions1={}
	twodigits2={}
	fourdigits2={}
	fourDigit_solutions2={}
	fourDigit_solutions3={}
	twodigits3={}
	finalAlleles={}

	#-------1st iteration-----------------------------------
	if hla_class==1:
		print "----------HLA class I------------"
	else:
		print "----------HLA class II-----------"
	sam1=runName+"-iteration1.sam"
	iteration=1
	print "First iteration starts....\nMapping ....."
	mapping(sam1,runName,readFile1,bowtiebuild,1,mapopt,trim3,gzip)
	medians={}
	for locus in locus_list:
		medians[locus]=0
	medianflag=False
	#Calculation of first digital haplotype.....
	output1=runName+".digitalhaplotype1"	
	print "Calculation of first digital haplotype....."
	(map,readcount,readspergroup)=createRefDict(hla1fasta,locus_list)
	readcount=readMapping(map,sam1,readcount)
	fourDigitString1=predictHLA(sam1,medians,output1,medianflag,locus_list,hla_class,readcount,readspergroup)
	print "1st iteration done.\nNow removing reads that mapped to the three top-scoring groups ......."
	try:
		removeReads(runName,createRemoveList(runName,map,locus_list))
	except IOError:
                print "Nothing to remove\n"
	
	#------2nd iteration------------------------------------------
	print "Second iteration starts ....\n Mapping ....."
	medians={}
	iteration=2
	sam2=runName+"-iteration2.sam"
	newReadFile1=runName+"-2nditeration.fq"
	mapping(sam2,runName,newReadFile1,bowtiebuild,2,mapopt,trim3,gzip)
	medianfile=runName+".digitalhaplotype1"
	row=2
	for locus in locus_list:
		medians[locus]=linecache.getline(medianfile, row).split('\t', 3)[2]
		row+=1
	medianflag=True
	output2=runName+".digitalhaplotype2"
	finaloutput=runName+".HLAgenotype2digits"
	#Calculation of second digital haplototype
	print "Calculation of second digital haplotype....."
	(map,readcount,readspergroup)=createRefDict(hla1fasta,locus_list)
	readcount=readMapping(map,sam2,readcount)
	fourDigitString2=predictHLA(sam2,medians,output2,medianflag,locus_list,hla_class,readcount,readspergroup)
	print "2nd iteration done."
	reportHLAgenotype(output1,output2,finaloutput,locus_list)
	print "Calculation of locus-specific expression ..."
	try:
		expression(locus_list,length_dict,map,runName,readFile1,hla_class)
	except IOError:
                tmp=""
		for locus in locus_list:
			tmp+=locus+": 0 RPKM\n"
		print tmp
	
	#-----3rd iteration in case of at least one homozygous call-----------
	f1=fourDigitString1.split(",")
	f2=fourDigitString2.split(",")
	alleleTwoDigit_index=0
	alleleFourDigit_index=1
	numberSolutions_index=2
	for locus in locus_list:
		twodigits1[locus]=f1[alleleTwoDigit_index]
		fourdigits1[locus]=f1[alleleFourDigit_index]
		if twodigits1[locus]=="no":
			fourDigit_solutions1[locus]=f1[numberSolutions_index]
		else:
			fourDigit_solutions1[locus]=int(f1[numberSolutions_index])
		fourDigit_solutions2[locus]=f2[numberSolutions_index]
		twodigits2[locus]=f2[alleleTwoDigit_index]
		fourdigits2[locus]=f2[alleleFourDigit_index]
		alleleTwoDigit_index+=3
		alleleFourDigit_index+=3
		numberSolutions_index+=3
	
	#try-catch block to prevent IO-Error in case of no expression
	try:
		if "no" in twodigits2.values():
			removeReads(runName,createRemoveListFourDigits(runName,map,fourdigits1))
			sam3=runName+"-iteration3.sam"
			mapping(sam3,runName,newReadFile1,bowtiebuild,2,mapopt,trim3,gzip)
			readcount=readMapping(map,sam3,readcount)
			output3=runName+".digitalhaplotype3"
			fourDigitString3=predictHLA(sam3,medians,output3,medianflag,locus_list,hla_class,readcount,readspergroup)
			f3=fourDigitString3.split(",")
			alleleTwoDigit_index=0
			numberSolutions_index=2
			for locus in locus_list:
				twodigits3[locus]=f3[alleleTwoDigit_index]
				fourDigit_solutions3[locus]=f3[numberSolutions_index]
				alleleTwoDigit_index+=3
				numberSolutions_index+=3
		#1st allele--------------------------------
		for locus in locus_list:
			allele1=fourdigits1[locus]
			if not fourDigit_solutions1[locus]=="":
				#if fourDigit_solutions1[locus]>1:
					#allele1+="'"
				finalAlleles[locus]=allele1+"\t"+str(getMaxP(sam1+".4digits"+locus+"1.solutions"))+"\t"
			else:
				finalAlleles[locus]="no\tNA\t"
	
		#2nd allele--------------------------------	
		alleleFourDigit_index=1
		for locus in locus_list:
			allele2=twodigits2[locus]
			if allele2=="no":
				allele3=twodigits3[locus]
				if allele3=="no" or fourDigit_solutions1[locus]==1 or not allele3.split(":")[0]==fourdigits1[locus].split(":"):
					if fourDigit_solutions1[locus]=="":
						finalAlleles[locus]+="no\tNA"
					else:
						finalAlleles[locus]+=fourdigits1[locus]+"\t"+str(getP(locus,finaloutput))
				else:
					allele3=f3[alleleFourDigit_index]
					#if int(f3[alleleFourDigit_index+1])>1:
						#allele3+="'"
					finalAlleles[locus]+=allele3+"\t"+str(getMaxP(sam3+".4digits"+locus+"2.solutions"))
					alleleFourDigit_index+=3
			else:
				allele2=fourdigits2[locus]
				#if int(fourDigit_solutions2[locus])>1:
					#allele2+="'"
				finalAlleles[locus]+=allele2+"\t"+str(getMaxP(sam2+".4digits"+locus+"2.solutions"))
				alleleFourDigit_index+=3
	
			
		finaloutput4digit=runName+".HLAgenotype4digits"
		#output the 4-digit type (stdout and to file)
		reportHLA4digitGenotype(finalAlleles,finaloutput4digit,locus_list,hla_class)
	except IOError:
		#no expression"
                print "no expression"
		finaloutput4digit=runName+".HLAgenotype4digits"
		for locus in locus_list:
			finalAlleles[locus]="no\tNA\tno\tNA"
		reportHLA4digitGenotype(finalAlleles,finaloutput4digit,locus_list,hla_class)

	cleanup(newReadFile1,runName)


def cleanup(newReadFile1,runName):
	os.remove(newReadFile1)
	for file in glob(runName+"*.sam"):
		os.remove(file)
	for file in glob(runName+"*.4digits*"):
		os.remove(file)
	for file in glob(runName+"*.aligned*"):
		os.remove(file)
	for file in glob(runName+"*.digitalhaplo*"):
		os.remove(file)
	for file in glob(runName+"*.readspergroup*"):
		os.remove(file)
	#for file in glob(runName+"*.ambiguity"):
		#os.remove(file)
  
def merge(runName):
            file_list=[runName+"-ClassI.bowtielog",runName+"-ClassII.bowtielog"]
            with open(runName+".bowtielog.txt", 'w') as outfile:
                for fname in file_list:
                        if fname==runName+"-ClassI.bowtielog":
                            with open(fname) as infile: 
                                outfile.write("----------HLA class I------------\n")
                                outfile.write(infile.read())
                        else:
                            with open(fname) as infile: 
                                outfile.write("----------HLA class II-----------\n")
                                outfile.write(infile.read())
            for file in glob(runName+"-Class*.bowtielog"):
                os.remove(file)
            file_list=[runName+"-ClassI.expression",runName+"-ClassII.expression"]
            with open(runName+".expression.txt", 'w') as outfile:
                for fname in file_list:
                    with open(fname) as infile: 
                        outfile.write(infile.read())
            for file in glob(runName+"-Class*.expression"):
                os.remove(file)
            file_list=[runName+"-ClassI.HLAgenotype4digits",runName+"-ClassII.HLAgenotype4digits"]
            with open(runName+".HLAgenotype4digits.txt", 'w') as outfile:
                for fname in file_list:
                    with open(fname) as infile: 
                        outfile.write(infile.read())
            for file in glob(runName+"-Class*.HLAgenotype*digits"):
                os.remove(file)

#performs the bowtie mapping for the 2 iterations using the given parameters
def mapping(sam,runName,readFile1,bowtiebuild,iteration,mapopt,trim3,gzip):
	if iteration==1:
		if gzip==0:
			mapping_cmd="(bowtie -3 "+trim3+" -S "+mapopt+" --al "+runName+".aligned "+bowtiebuild+" "+readFile1+" |  awk -F \'\\t\' '$3 != \"*\"{ print $0 }' > "+sam+") 2> "+runName+".bowtielog"
		else:
			mapping_cmd="(bowtie -3 "+trim3+" -S "+mapopt+" --al "+runName+".aligned "+bowtiebuild+" <(zcat "+readFile1+") |  awk -F \'\\t\' '$3 != \"*\"{ print $0 }' > "+sam+") 2> "+runName+".bowtielog"
	if iteration==2:
		mapping_cmd="bowtie -3 "+trim3+" -S "+mapopt+" "+bowtiebuild+" "+readFile1+" "+sam
	#execute bowtie
	process_mapping=subprocess.Popen(['bash','-c',mapping_cmd],stdout=subprocess.PIPE)
	out,err=process_mapping.communicate()
	
	if iteration==1:
		#print alignment stats
		printcommand="cat "+runName+".bowtielog"
		printcommand_proc=subprocess.Popen(['bash','-c',printcommand],stdout=subprocess.PIPE)
		out,err=printcommand_proc.communicate()
		print out

#create dictionary "map", that contains all IMGT/HLA-allele names as keys and the allele name (e.g. A*02:01:01) as value
#dictionary "readcount" is initialized with 0 for each allele
#dictionary "readspergroup" is initialized with 0 for each group (2digit, e.g. A*01)
#dictionary "allelesPerLocus" stores the number of alleles per locus.
def createRefDict(hlafasta,class1_list):
	map={}
	allelesPerLocus={}
	readcount={}
	readspergroup={}
	
	for locus in class1_list:
		allelesPerLocus[locus]=0
	handle=open(hlafasta,'r')
	for record in SeqIO.parse(handle, "fasta") :
		l=record.description.split(' ')
		hlapseudoname=l[0]
		hlaallele=l[1]
		locus=hlaallele.split('*')[0]
		if locus in class1_list:
			map[hlapseudoname]=hlaallele
			readcount[hlaallele]=0
			readspergroup[hlaallele.split(":")[0]]=0
			allelesPerLocus[hlaallele.split('*')[0]]+=1
	handle.close()
	return (map,readcount,readspergroup)

#Open sam-file and count mappings for each allele 
def readMapping(map,sam,readcount):
	samhandle=open(sam,'r')
	for line in samhandle:
		if line[0]!='@':
			l=line.split('\t')
	                hla=l[2]
                	readcount[map[hla]]+=1
	return readcount

#predict HLA type 
def predictHLA(sam,medians,output,medianflag,locus_list,hla_class,readcount,readspergroup): 
	maxAlleles={}
	maxkey={}
	allele_dict={}
	fourdigits_dict={}	
	alleleVector={}
	alleleCount={}
	readspergroup_dict_list={}
	readspergroup_dict={}
	fourdigits_dict={}
	fourdigits_sorted_dict={}

	for locus in locus_list:
		readspergroup_dict_list[locus]=[]
		readspergroup_dict[locus]={}
		fourdigits_dict[locus]={}
		fourdigits_sorted_dict[locus]={}
		alleleVector[locus]=""

	maxallelepergroup={}
	
	#for each allele, to which at least 1 read map, find the allele which has the most reads within a group (2-digit-level, e.g. A*02") and save
	#i) the key of this allele as ambassador for the group => maxallelepergroup holds for each group the top-scoring allele
	#ii) the number of reads mapping to the top-scoring allele => readspergroup
	for key in readcount:
		if readcount[key] > 0:
			group=key.split(":")[0]
			if readspergroup[group] <=readcount[key]:
				readspergroup[group]=readcount[key]
				maxallelepergroup[group]=key

	readspergrouphandle=open(sam+".readspergrouphandle","a")
	for group in readspergroup:
		readspergrouphandle.write(group+"\t"+str(readspergroup[group])+"\n")
	readspergrouphandle.close()
	#consider all A-,B-, and C-groups seperately
	#readspergroup<A|B|C>list = list of all reads mapping to the top-scoring groups minus the decision threshold (which is 0 in the first iteration)
	#readspergroup<A|B|C> = contains the same entries as the list, but the entries are uniquely accessible via the group-key (e.g. B*27)
	for key in readspergroup:
		locus=key.split("*")[0]
		readspergroup_dict_list[locus].append(readspergroup[key]-float(medians[locus]))
		readspergroup_dict[locus][key]=readspergroup[key]-float(medians[locus])
	
	#Determine top-scoring group of the whole locus (A,B,C) and store it
	#maxkey<A,B,C> = group (e.g. A*02) with the most reads
	#It can be that, e.g. in cancer cells A whole locus is lost. For that reason it is checked if the 
	#number of mapping reads of the top-scoring group maxkey<A,B,C> is > 0, otherwise "no" ist reported for this locus
	for locus in locus_list:
		if len(readspergroup_dict_list[locus])>0:
			maxkey[locus]=max(readspergroup_dict[locus],key=lambda a:readspergroup_dict[locus].get(a))
			if readspergroup_dict[locus][maxkey[locus]] > 0:
				maxAlleles[locus]=maxallelepergroup[maxkey[locus]]
				allele_dict[locus]=maxkey[locus]
			else:
				allele_dict[locus] = "no"
		else:
			allele_dict[locus] = "no"
	
	for key in readcount:
		for locus in locus_list:	
			if key.split(":")[0]==maxkey[locus]:
				fourdigits_dict[locus][key]=readcount[key] 
	
	for locus in locus_list:
		fourdigits_sorted_dict[locus]=sorted(fourdigits_dict[locus].items(), key=itemgetter(1),reverse=True)
	
		if medianflag:
			iteration=2
		else:
			iteration=1
	
		fourdigit_writehandle=open(sam+".4digits"+locus+str(iteration),"w")
		for key in fourdigits_sorted_dict[locus]:
			fourdigit_writehandle.write(key[0]+": "+str(key[1])+"\n")
		fourdigit_writehandle.close()

		readspergrouplocus=sorted(readspergroup_dict[locus].items(), key=itemgetter(1),reverse=True)

		if medianflag:
			#in the 2nd iteration: 
			#1.) DO NOT remove the top-scoring group from the set <a,b,c>vec as this is more strict when calculating the probability of the top scoring group being an outlier
			#The strings <a,b,c>vec are used by the R-script to calculate the probability of the top scoring group being an outlier
			for key in maxallelepergroup:
				if key.split("*")[0]==locus:
					alleleVector[locus]+=str(readcount[maxallelepergroup[key]])+","
			#2.) Add the decision thresholds to the sets <a,b,c>vec, as this enables measuring the distance of the top-scoring group to this distance
			if allele_dict[locus]=="no":
				alleleVector[locus]+=str(medians[locus])
			#3.) In case of, e.g. a loss of whole HLA-locus (A,B,C), <a,b,c>vec only contain median[<0|1|2>].
			#To avoid errors in the R-Script, set to 0
			if alleleVector[locus]=="" or alleleVector[locus]==medians[locus]:
				alleleVector[locus] = "0"
				alleleCount[locus]= "0"
			else:
				alleleCount[locus]=str(readcount[maxallelepergroup[maxkey[locus]]])

		else:
			#in the 1st iteration: remove the top-scoring group from the set <a,b,c>vec as this increases the certainty when calculating the probability of the top scoring group being an outlier
			for key in maxallelepergroup:
				if key.split("*")[0]==locus:
					if key!=maxkey[locus]:
						alleleVector[locus]+=str(readcount[maxallelepergroup[key]])+","
			#2.) DO NOT add the decision thresholds to the sets <a,b,c>vec
			if alleleVector[locus]=="":
                                alleleVector[locus] = "0"
                                alleleCount[locus]= "0"
                        else:
				alleleCount[locus]=str(readcount[maxallelepergroup[maxkey[locus]]])

	rstring=""
	numberArgs=0
	for locus in locus_list:
		numberArgs+=3
		rstring+=str(alleleCount[locus])+" "+alleleVector[locus]+" "+maxkey[locus]+" "
	#call R-script "commmand.R" to calculate the confidence of the top-scoring allele
	routput=os.popen("R --vanilla < "+os.path.abspath(os.path.dirname(sys.argv[0]))+"/command.R --args "+str(numberArgs)+" "+rstring).read()
	parseOutput=routput.split("\n")

	entries = []
	for entry in parseOutput:
		if entry[0:3]=="[1]":
			entries.append(str(entry[4:len(entry)]))
		
	fourDigitString=""
	if medianflag:
		iteration=2
	else:
		iteration=1

	entryIndex=1
	for locus in locus_list:
		if not allele_dict[locus]=="no":
			fourDigitString+=allele_dict[locus]+","+fourdigits.determine4digits_main(sam+".4digits"+locus+str(iteration),alleleVector[locus],runName,hla_class)+","
		else:
			fourDigitString+=allele_dict[locus]+",,,"
			if entries[entryIndex]!="NA":
	                        entries[entryIndex]=str(1-float(entries[entryIndex]))
		entryIndex+=2
			
	pred2File(entries,readspergroup_dict,output,allele_dict,hla_class,locus_list)
	return fourDigitString
	
#write digital haplotype into file
def pred2File(entries,readspergroup_dict,output,allele_dict,HLAclass,locus_list):
	out = open(output, 'w')
        out.write("HLA\tHLA1\tmedian-Value\talternative\tp-Value\n")
	index=0
	for locus in locus_list:
		out.write(locus+"\t"+allele_dict[locus]+"\t")
		#compute the decision threshold (median) for homozygosity vs. heterozygosity for the second iteration 
		out.write(str(numpy.median(readspergroup_dict[locus].values()))+"\t")
		out.write(entries[index]+"\t"+entries[index+1]+"\n")
		index+=2

	out.close()
	#print "The digital haplotype is written into "+output

#open mapping file and all read ids to the list "removeList", which map to one of the three groups in "alleles"
def createRemoveList(runName,map,locus_list):
	removeList={}
	alleles = []
	sam=runName+"-iteration1.sam"
	alleles_in=runName+".digitalhaplotype1"
	line=2
	for locus in locus_list:
		alleles.append(linecache.getline(alleles_in, line).split('\t', 2)[1])
		line+=1

	samhandle=open(sam,'r')
	for line in samhandle:
		if line[0]!='@':
			illuminaid=line.split("\t")[0]
			hlapseudoname = line.split("\t")[2]
			if map[hlapseudoname].split(':')[0] in alleles:
				removeList[illuminaid]=1
	samhandle.close()
	return removeList

#open mapping file and all read ids to the list "removeList", which map to one of the three four-digit alleles in "alleles"
def createRemoveListFourDigits(runName,map,fourdigits1_dict):
	removeList={}
	sam=runName+"-iteration1.sam"

	samhandle=open(sam,'r')
	for line in samhandle:
		if line[0]!='@':
			illuminaid=line.split("\t")[0]
			hlapseudoname = line.split("\t")[2]
			fourdigits=map[hlapseudoname].split(':')[0]+":"+map[hlapseudoname].split(':')[1]
			if fourdigits in fourdigits1_dict:
				removeList[illuminaid]=1
	samhandle.close()
	return removeList	
	
#Remove reads that mapped to the three top-scoring alleles and write the remaining reads into two new read files
def removeReads(runName,removeList):
	aligned1=runName+".aligned"
	newReadFile1=runName+"-2nditeration.fq"
	#r1 which is the input of bowtie in the 2nd iteration
	r1=open(newReadFile1,"w")
	#open the file, that contain the reads that mapped in the 1st iteration
	aligned_handle1=open(aligned1,"r")
	
	#One read entry consists of 4 lines: header, seq, "+", qualities.
	for record in SeqIO.parse(aligned_handle1, "fastq"):
		illuminaid=record.id.split('/')[0].split(' ')[0]#find exact id, which also appears in the mapping file
		if not illuminaid in removeList:
			SeqIO.write(record, r1, "fastq")

#write the final prediction (both digital haplotypes) to file and stdout
def reportHLAgenotype(output1,output2,finaloutput,locus_list):
	filehandle1=open(output1,'r').readlines()[1:len(locus_list)+1]
	filehandle2=open(output2,'r').readlines()[1:len(locus_list)+1]
	outfile=open(finaloutput, 'w')
	outfile.write("#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n")
	
	for i in range(len(filehandle1)):
		filehandle1[i]=filehandle1[i][0:-1]
	for i in range(len(filehandle2)):
		filehandle2[i]=filehandle2[i][0:-1]
	print "-----------2 digit typing results-------------"
        print "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence"
	line=0
	for locus in locus_list:
		allele1 = filehandle1[line].split('\t',2)[1]
		allele1_score = filehandle1[line].split('\t')[4]
		allele2 = filehandle2[line].split('\t',2)[1]
		if allele2 == "no":
			allele2 = "hoz("+filehandle2[line].split('\t')[3]+")"
		allele2_score = filehandle2[line].split('\t')[4]
		line+=1
		#write complete HLA genotype to file
		outfile.write(locus+"\t"+allele1+"\t"+allele1_score+"\t"+allele2+"\t"+allele2_score+"\n")
		#.. and print it to STDOUT
		print locus+"\t"+allele1+"\t"+allele1_score+"\t"+allele2+"\t"+allele2_score
	outfile.close()

def reportHLA4digitGenotype(finalAlleles,finaloutput4digit,locus_list,hla_class):
	#write complete HLA genotype at four-digit-level to file
	outfile=open(finaloutput4digit, 'w')
	if hla_class==1:
		outfile.write("----------HLA class I------------\n")  
	else:
		outfile.write("----------HLA class II-----------\n")
	outfile.write("#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n")
	 #.. and print it to STDOUT
        print "-----------4 digit typing results-------------"
	print "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence"
	for locus in locus_list:
		outfile.write(locus+"\t"+finalAlleles[locus]+"\n")
		print locus+"\t"+finalAlleles[locus]
	
#calculate locus-specific expression
def expression(locus_list,length_dict,map,runName,readFile1, hla_class):
	outfile=open(runName+".expression", 'w')
	if hla_class==1:
		outfile.write("----------HLA class I------------\n")  
	else:
		outfile.write("----------HLA class II-----------\n") 
	aligned1=runName+".aligned"
	sam=runName+"-iteration1.sam"
	logfile=runName+".bowtielog"
	print logfile
	totalreads=float(linecache.getline(logfile,1).split(':')[1])
	alleles_in=runName+".digitalhaplotype1"
	alleles=[]
	line=2
	for locus in locus_list:
		alleles.append(linecache.getline(alleles_in, line).split('\t', 2)[1])
		alleles.append(linecache.getline(alleles_in, line).split('\t')[3])
		line+=1
	
	#create read dictionary
	reads={}
	aligned_handle1=open(aligned1,"r")
	for record in SeqIO.parse(aligned_handle1, "fastq"):
		illuminaid=record.id.split('/')[0].split(' ')[0]#find exact id, which also appears in the mapping file
		reads[illuminaid]={}
		for locus in locus_list:
			reads[illuminaid][locus]=0
	samhandle=open(sam,'r')
	for line in samhandle:
		if line[0]!='@':
			illuminaid=line.split("\t")[0].split('/')[0].split(' ')[0]
			hlapseudoname = line.split("\t")[2]
			if map[hlapseudoname].split(':')[0] in alleles:
				reads[illuminaid][map[hlapseudoname].split('*')[0]]+=1
	count={}
	for locus in locus_list:
		count[locus]=0
	for key in reads:
		n=0
		for locus in reads[key]:
			if reads[key][locus] > 0:
				n+=1
		for locus in reads[key]:
			if reads[key][locus] > 0:
				count[locus]+=float(1.0/float(n))
	
	#Calculate RPKM and print expression values for each locus to stdout
	for locus in count:
		print locus+": "+str(round(float((1000.0/length_dict[locus]))*float((1000000.0/totalreads))*count[locus],2))+" RPKM"
		outfile.write(locus+": "+str(round(float((1000.0/length_dict[locus]))*float((1000000.0/totalreads))*count[locus],2))+" RPKM\n")
        outfile.close()

#In case of ambiguous typings, the allele(s) with the best p-value (which is actually the smallest one, so the name of the function is misleading - sorry) is reported and thus the min(p) is returned
def getMaxP(file):
	p=[]
	for line in open(file,"r"):
		if not line[0]=="#":
			if line.split("\t")[1][0:-1]=="NA":
				return "NA"
			p.append(float(line.split("\t")[1][0:-1]))
	try:
		return min(p)
	except ValueError:
		return "NA"

#Return the p value of a prediction, which is stored in the intermediate textfile
def getP(locus,finaloutput):
    if genes==5:
        	if locus=="A" or locus=="DRB1":
        		p=linecache.getline(finaloutput, 2).split('\t')[4]
        	if locus=="B" or locus=="DPB1":
        		p=linecache.getline(finaloutput, 3).split('\t')[4]
        	if locus=="C":
        		p=linecache.getline(finaloutput, 4).split('\t')[4]
        	return p[0:-1]
    else: 
        	if locus=="A" or locus=="DQB1":
        		p=linecache.getline(finaloutput, 2).split('\t')[4]
        	if locus=="B" or locus=="DRB1":
        		p=linecache.getline(finaloutput, 3).split('\t')[4]
        	if locus=="C" or locus=="DPB1":
        		p=linecache.getline(finaloutput, 4).split('\t')[4]
        	return p[0:-1]
  
if __name__ == '__main__':
	length_dict1={}
	length_dict2={}
	parser = OptionParser(usage="usage: %prog -1 readFile1 -f readFolder -r runName [-p <int>] [-3 <int>]")
	parser.add_option("-1",
			action="store", 
			dest="readFile1",
			help="File name of #1 mates (uncompressed or gzipped fastq)")
	parser.add_option("-f",
			action="store", 
			dest="readFolder",
			help="Folder with RNA-seq data (uncompressed or gzipped fastq files)")
	parser.add_option("-r", "--runName",
			action="store", 
			dest="runName",
			help="Name of this HLA typing run. Wil be used throughout this process as part of the name of the newly created files.")	
	parser.add_option("-p", "--threads",
                        action="store",
                        dest="threads",
                        default="6",
                        help="Bowtie option: Launch <int> parallel search threads. Default (from seq2HLA): 6")
	parser.add_option("-3", "--trim3",
			action="store",
			dest="trim3",
			default="0",
			help="Bowtie option: -3 <int> trims <int> bases from the low quality 3' end of each read. Default: 0")
	parser.add_option("-g", "--genes",
			action="store",
			dest="genes",
			default="5",
			help="Number of HLA genes included for typing. Default: 5, maybe increased to 6 (adding DQB1)")

	(options, args) = parser.parse_args()
	if not options.readFile1:   
		parser.error('File name is not given.')
	if not options.runName:   
		parser.error('Run name not given.')
	readFile1=(options.readFolder+"/"+options.readFile1)
	#change the output directory
	runName=(options.readFolder+"/"+options.runName)
	bowtiebuildClassI=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/ClassIWithoutNQex2-3.plus75"
	fastaClassI=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/ClassIWithoutNQex2-3.plus75.fasta"
	trim3=str(options.trim3)
	threads=str(options.threads)
	class1_list=["A","B","C"]
	for locus in class1_list:
		length_dict1[locus]=694	
  	genes=int(options.genes)
  	if genes==5:
		class2_list=["DRB1","DPB1"]
		length_dict2["DRB1"]=421
		length_dict2["DPB1"]=414
		bowtiebuildClassII=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/HLA2.ex2.plus75_short"
		fastaClassII=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/HLA2.ex2.plus75_short.fasta"	
  	else:
		class2_list=["DQB1","DRB1","DPB1"]
		length_dict2["DQB1"]=421
		length_dict2["DRB1"]=421
		length_dict2["DPB1"]=414
		bowtiebuildClassII=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/HLA2.ex2.plus75"
		fastaClassII=os.path.abspath(os.path.dirname(sys.argv[0]))+"/references/HLA2.ex2.plus75.fasta"
 
main_single(runName,readFile1,fastaClassI,fastaClassII,bowtiebuildClassI,bowtiebuildClassII,trim3,threads,genes,class1_list,class2_list,length_dict1,length_dict2)