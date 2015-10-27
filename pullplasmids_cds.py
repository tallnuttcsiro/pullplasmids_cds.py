from __future__ import print_function
from Bio import SeqIO
import sys

#Script to remove sequences from a file based on a list of sequence ids from a blast file (T. Allnutt, CSIRO, 2015)
#usage: pullplasmids_cds inputfile_cds.faa blast_format6file.txt plasmids_cds.fasta cds_no_plasmids_output.fasta
#requires module: Biopython
#python >= v2.7



inputfile1 = sys.argv[1] #input file

reads = sys.argv[2]# blast file against the plasmid reference

h=open(sys.argv[4],'w')#file of input with items from the list file removed

outputfile2 = sys.argv[3] #file plasmid reads pulled out

fmt = sys.argv[5]
g = open(outputfile2,'w')

ids = open(reads,'r')

ilist=[]
c=0
k=0
n=0
y=""

phits={}

print ("reading blast..")
#add all the hits in a scaffold then compare to total scaffold length. if >0% then call plasmid?
for line in ids:
	n = line.split("\t")[0]
	b = float(line.split("\t")[-1])
	
	if n not in phits.keys():
		phits[n]=b
	else:
		phits[n]=phits[n]+b


#print len(ilist)
x = 0
print ("reading file to memory..")
record = SeqIO.to_dict(SeqIO.parse(inputfile1, fmt))
hitnames=""
#print (record.keys())

hits=open("tmp.hits",'a')
c1=0
for y in phits.keys():
	
	if y in record:
		
		scaflen=float(len(record[y].seq))
		pc= phits[y]/scaflen *100
		
		if pc >= 0:
			c1=c1+1
			print (y, phits[y])
			print(pc,"%")
			hitnames=hitnames+y+"\t"
			x=x+1
			SeqIO.write(record[y],g,fmt) #write  contigs from list
			#remove y from dict
			del record[y]
hits.write(inputfile1+"\t"+hitnames+"\n")
print (c1, "plasmid hits")	
#write records that were not in list
for i in record.keys():

	SeqIO.write(record[i],h,fmt)
		
	
ids.close()

g.close()









	

		






