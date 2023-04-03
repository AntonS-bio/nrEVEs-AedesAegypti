from calendar import c
import subprocess
import multiprocessing as mp
import numpy as np
from os import listdir
from os.path import isfile, join
from Bio import SeqIO

assemblyDir="/mnt/storage5/anton/Mosquitoes/DENV/assemblies/"
wd="/mnt/storage5/anton/Mosquitoes/DENV/"
trimmedAssemblies="/mnt/storage5/anton/Mosquitoes/DENV/trimmedAssemblies/"
assemblies = [f for f in listdir(assemblyDir) if isfile(join(assemblyDir, f)) and f.find(".fasta")>-1]

##IF blastDB does not yet exist, make one using:
##makeblastdb -in GCF_002204515.2_AaegL5.0_genomic.fna -title AaegL5 -out ./blastDB/AaegL5 -dbtype nucl  -blastdb_version 4 

def getBlastResults(queryFile, db):
    blastResults=[]
    blastCommand=f"blastn -query "+queryFile+" -task 'blastn'  -max_target_seqs 10000 -db "+db+" -num_threads 11 -evalue 0.01 -word_size 11 -outfmt '6 delim=  qseqid qstart qend'"
    blastOutput=subprocess.run(blastCommand, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in blastOutput.stdout.decode().split("\n"):
        if len(line)>0:
            contig, start, stop=line.split("\t")
            blastResults.append( [contig, int(start), int(stop)])
    blastResults.sort(key=lambda t: (t[0], t[1], t[2]))
    bedInput=""
    for value in blastResults:
        bedInput=bedInput+'\t'.join(str(f) for f in value )+"\n"
    return bedInput

def mergeBlastResults(toInclude, toExclude):
    intervalsToInclude={}
    intervalsToExclude={}
    for bedToolsOutput, mergedIntervals in zip([toInclude, toExclude],[intervalsToInclude,intervalsToExclude]):
        if bedToolsOutput.stdout is None:
            continue
        for line in bedToolsOutput.stdout.decode().split("\n"):
            if len(line)==0:
                continue
            contig, start, end=line.split("\t")
            if contig not in mergedIntervals:
                mergedIntervals[contig]=[]
            mergedIntervals[contig].append([int(start), int(end)])
    retainedIntervals={}
    for contig in intervalsToInclude:
        hasExclusions=False
        for interval in intervalsToInclude[contig]:
            intervalRange=range(interval[0],interval[1]+1)
            if contig in intervalsToExclude:
                for excludedInterval in intervalsToExclude[contig]:
                    intervalRange=np.setdiff1d(intervalRange, range(excludedInterval[0],excludedInterval[1]+1) )
                    if len(intervalRange)< len(range(interval[0],interval[1]+1)): #exclusions only exists if two intervals overlap
                        hasExclusions=True
            #for speed, the final determination of intervals must be split in cases where contig has exclusions and where it does not
            if len(intervalRange)<25:#the remaining sequence is too short
                continue
            intervalStart=intervalRange[0]
            if not hasExclusions:
                retainedIntervals[contig]=[]
                for value in intervalsToInclude[contig]:
                    retainedIntervals[contig].append( value )
            else:
                for i in range(1, len(intervalRange) ):
                    if intervalRange[i-1]+1!=intervalRange[i] or i==len(intervalRange)-1:
                        intervalEnd=intervalRange[i] if i==len(intervalRange)-1 else intervalRange[i-1]
                        if intervalEnd-intervalStart>25: #this simply adds contig interval
                            if contig not in retainedIntervals:
                                retainedIntervals[contig]=[]
                            retainedIntervals[contig].append( [intervalStart, intervalEnd ] )
                        intervalStart=intervalRange[i]

                
    return retainedIntervals

##IF viruses blastDB does not yet exist, make one using:
##makeblastdb -in ~/Mosquitoes/DENV/InputData/RefSeq_ExEVEs.fa -title RefSeq_ExEVEs -out ~/Mosquitoes/DENV/InputData/RefSeq_ExEVEs -dbtype nucl  -blastdb_version 4 

counter=0
blastResults={}
for assembly in assemblies:
    print(assembly)
    if counter % 20==0:
        print(counter/len(assemblies))
    counter+=1
    
    #don't reinvent bedtools
    bedToInclude=getBlastResults(assemblyDir+assembly,"~/Mosquitoes/DENV/InputData/RefSeq_ExEVEs")
    with open(wd+"temp.fasta","w") as output:
        output.write(bedToInclude)
    bedtoolsIncludeOutput=subprocess.run("bedtools merge  -d 10 -i "+wd+"temp.fasta", shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bedToExclude=getBlastResults(assemblyDir+assembly,"~/DRerio/DRerio")
    if len(bedToExclude)!=0:
        with open(wd+"temp.fasta","w") as output:
            output.write(bedToExclude)
        bedtoolsExcludeOutput=subprocess.run("bedtools merge  -d 5 -i "+wd+"temp.fasta", shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        bedtoolsExcludeOutput=subprocess.CompletedProcess(args=['ls', '-l'], returncode=0)

    mergedIntervals=mergeBlastResults(bedtoolsIncludeOutput, bedtoolsExcludeOutput)

    with open(trimmedAssemblies+assembly,"w") as trimmedFasta:
        for record in SeqIO.parse(assemblyDir+assembly, "fasta"):
            if record.id in mergedIntervals:
                trimmedFasta.write(">"+assembly.replace(".fasta","")+"_"+record.id+"\n")
                for interval in mergedIntervals[record.id]:
                    trimmedFasta.write(str(record.seq)[interval[0]-1:interval[1]]+"\n")

#cat *.fasta >> ~/Mosquitoes/DENV/DB/allTrimmed.fasta

#makeblastdb -in allTrimmed.fasta -title allTrimmed -out allTrimmed -dbtype nucl  -blastdb_version 4  

#blastn -query allTrimmed.fasta -task 'blastn' -max_target_seqs 1000000000 -db allTrimmed -num_threads 8 -evalue 0.01 -word_size 11 -outfmt '6 delim=  qseqid qstart qend sseqid sstart send pident evalue' > ../allTrimmedBlast.tsv
                