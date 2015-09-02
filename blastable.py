#!/usr/bin/env python

#Script to identify the presence of a panel of query genes across a large number of whole genomes. 


import sys
import re
import os
import glob
import subprocess
from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter
import pandas as pd

dic = {}
plot = {}
dic['1-Querylength'] = {}
plot['1-Querylength'] = {}
blast_tab_list = []
scratch = os.path.basename("scratch") #specify blast output folder
suffix = ''


def main():
    
    #read in input query for blast
    queryFile = args.input
    
    
    #read dir with genomes and make blast do blast
    listGenomes = parseGenomes(args.genomes)
    for genome in listGenomes:
        files_to_blast = [queryFile, genome]
        returned = doBlast(files_to_blast)
        
        blast_tab_file = returned[0]
        suffix = returned[1]
        blast_tab_list.append(blast_tab_file)
        
    for blast_result in blast_tab_list: #iterate through blast results list
#         blast_result = blast_result.rstrip()
        genome_name = blast_result.split('.vs.')[1] #reverse format blast_result name to parse out genome name
        genome_name = re.sub('\.' + suffix, '', genome_name)
        
        if blast_result not in dic: #initialise dic for blast result
            dic[genome_name] = {}
            plot[genome_name] = {}
    
#         print "Looking for blast hit in " + scratch + "/" + blast_result
#         if not os.path.exists(scratch): os.makedirs(scratch) and warning("Output folder does not exist, creating one...")
#         resultFH = open(scratch + "/" + blast_result, 'r') #open blast result tab file with the directory prefix
#     
        resultFH = open(blast_result, 'r')
        for hit in resultFH: #go through hits in a blast_result file
            process_hit(hit, blast_result, genome_name) #feed the hitline and the name of the blast result file to process_hit()
            #print "Processing: " + blast_result
    df = pd.DataFrame(dic).T
    df.fillna(0, inplace=True)
    df.to_csv("allHits.csv")
    
    plotdf = pd.DataFrame(plot).T
    plotdf.fillna(0, inplace=True)
    plotdf.to_csv("topHit.csv")
    
    
    
    if args.plot: # do something with matplotlib here to generate a heatmap
        pass
    

# def doBlast(db):
#     print "Doing blast on %s" % (db)
#     
def doBlast(inputList):
    for file in inputList:
        if determineFileType(file) in {"genbank", "embl"}: #if file is genbank, convert to fasta
            convert2Fasta(file)
        else:
            pass
        
#    for i in range(0,len(inputList) - 1): #pair up seqs for blast
    i = 0
    
    queryName = re.sub(r"\.(\w+$)", r"", inputList[i]) #strip suffixes from filename
    subjecName = re.sub(r"\.(\w+$)", r"", inputList[i+1])
    
    queryName = os.path.basename(queryName)
    subjecName = os.path.basename(subjecName)
    
    queryFile = ''
    subjecFile = ''
                   
    if determineFileType(inputList[i]) in {"genbank", "embl"}: #Query sequence for BLAST. if genbank, use fasta file generated earlier
        queryFile = re.sub(r"\.\w+$", r".fa", inputList[i])
    else:
        queryFile = inputList[i]
        
        #check if query is is a multi fasta
#             filetype = determineFileType(inputList[i])
#             records = list(SeqIO.parse(inputList[i], filetype))
#             if len(records) > 1:
#                 if args.verbose: print "Multifasta detected, concatenating prior to BLASTing"
#                 mergedFile = mergeRecords(inputList[i])
#                 queryFile = mergedFile
#                     
#             else:
#                 queryFile = inputList[i]
#                 pass
    
    if determineFileType(inputList[i+1]) in {"genbank", "embl"}: #Subject sequence for BLAST, if genbank, use fasta file generated earlier
        subjecFile = re.sub(r"\.\w+$", r".fa", inputList[i+1])
    else:
        subjecFile = inputList[i+1]
        #check if query is is a multi fasta 
#             filetype = determineFileType(inputList[i+1])
#             records = list(SeqIO.parse(inputList[i+1], filetype))
#             if len(records) > 1:
#                 if args.verbose: print "Multifasta detected, concatenating prior to BLASTing"
#                 mergedFile = mergeRecords(inputList[i+1])
#                 subjecFile = mergedFile
#                     
#             else:
#                 subjecFile = inputList[i+1]
#                 pass
    
    blastType = args.blast
    
    if args.verbose: print "Performing blast"
    blastOptionsPre = (args.flags if args.flags else "")
    blastOptions = re.sub(r"-(\w+)\s", r"\1_", blastOptionsPre)
    blastOptions = re.sub(r"\s+", r".", blastOptions)
    if args.verbose: print "with options: %s" % (blastOptionsPre)
    
    blast_out = "%s.vs.%s.%s.%s.tab.test" % (queryName, subjecName, blastOptions, blastType)
    suffix = "%s.%s.tab.test" % (blastOptions, blastType)
    
    if os.path.exists(blast_out): #check if blast output exists
        if args.verbose: warning("Existing blast results detected, skipping...")
        pass
    else:    
        #run BLAST
        #make blastDB
        #print "making blast db for: " + subjecFile
        subprocess.Popen("makeblastdb -dbtype nucl -in %s" % (subjecFile), shell=True).wait()
    
        subprocess.Popen("%s -query %s -db %s -outfmt '6 std qlen' -out %s %s" % (blastType, queryFile, subjecFile, blast_out, blastOptionsPre), shell=True).wait()
        #print "%s -query %s -subject %s -outfmt '6 std qlen' -out %s %s" % (blastType, queryFile, subjecFile, blast_out, blastOptionsPre) #uncomment to print blast command
    
    returnList = []
    returnList.append(blast_out)
    returnList.append(suffix)
    
    return returnList #return the blast output filename AND the suffix


            

def parseGenomes(dir):
    
    print "Looking for genomes/sequences in: " + os.path.join(dir, '')
    recognizedFileTypes = ('*.fas', '*.fna','*.fa','*.fasta', '*.gb', '*.gbk')
    files_grabbed = []
    for type in recognizedFileTypes:
        files_grabbed.extend(glob.glob(os.path.join(dir, '') + type))
        
    if not files_grabbed:  error("No fasta files found. Use fa, fna, fasta, fas, gbk or gb extension")
    
    return files_grabbed

def checkHit(hitline):
    pident = float(hitline[2])
    length = int(hitline[3])
    qlen = int(hitline[12])
    tol = ((pident*length)/qlen) #threshold value/formula
    
    tolValue = float(args.tol)
    
    if tol >= tolValue:
        #print "Tol: " + str(tol)
        #print "Pident: " + str(pident)
        return True
    else:
        return False
        
def process_hit(hit, blast_result, genome_name):
    elements = hit.rstrip().split('\t')
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen) = elements #assign hit line to meaningful variable names
    
    #Check and process qseqid
    
    if not qseqid.count(',') == 3: error("Please ensure that your fasta headers contain 4 comma-separated values. Please refer to the seqfindr input format and ensure no spaces occur before the 3 commas. BLAST stops parsing header" % (qseqid)) 
    qseqid_list = qseqid.split(',') #replace comma at the end of query ID. Based on a query header formatted for Seqfindr
    qseqid = qseqid_list[1]

    dic['1-Querylength'][qseqid] = qlen
    plot['1-Querylength'][qseqid] = qlen

    if checkHit(elements):
        pident = float(pident)
        length = int(length)
        qlen = int(qlen)
        tol = ((pident*length)/qlen) #threshold value/formula
        
        if args.verbose: print "Query %s has a hit above the cutoff in %s. TOL: %s. %s mismatches and %s gaps across %s for query of length %s" % (qseqid, sseqid, tol, elements[4], elements[5], elements[3], elements[12]) #print info about hit that passed filter

#This big parses all the hits showing the presence of multiple copies/hits
        if not qseqid in dic[genome_name]: #check if the query already has a previous hit
            tol = int(tol)
            tol = str(tol)
            dic[genome_name][qseqid] = tol
        else:                               #if there IS a previous hit, stick the two numbers together. 
             tol = int(tol)
             tol = str(tol)
             existingtol = str(dic[genome_name][qseqid])
             newtol = '%s:%s' % (existingtol, tol) #stick two tol values together
             dic[genome_name][qseqid] = newtol

#This bit only parses the top hit, ignoring all other hits even if they are above the tolValue

        if not qseqid in plot[genome_name]: #check if the query already has a previous hit
            plot[genome_name][qseqid] = int(tol)
        else:                               #if there IS a previous hit, stick the two numbers together. 
            plot[genome_name][qseqid] = int(tol)
#              tol = int(tol)
#              tol = str(tol)
#              existingtol = str(plot[genome_name][qseqid])
#              newtol = '%s:%s' % (existingtol, tol) #stick two tol values together
#              plot[genome_name][qseqid] = newtol
    
    return plot
    return dic

######################### Functions from BWAST ########################

def determineFileType(file):
    if re.search("|".join(["gb$", "gbk$"]), file):
        return "genbank"
    elif re.search("|".join(["embl$", "emb$"]), file):
        return "embl"
    elif re.search("|".join(["fasta$", "fa$", "fna$", "fas$"]), file):
        return "fasta"
    else:
        return False #defaults to fasta
        
        
def mergeRecords(file): #adapted from SeqHandler by NF Alikhan (github.com/happykhan/seqhandler)

#SeqHandler is a script for merging, converting and splitting sequence files (Genbank, EMBL, fasta and others). Please use it to merge multi-Genbank files before running bwast.py

    filetype = determineFileType(file) #determine file type
    readInMultifasta = open(file, "r")
    records = list(SeqIO.parse(readInMultifasta, filetype))

    mergingFile = records[0]
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    contigs = SeqFeature(FeatureLocation(0, len(mergingFile) ), type="fasta_record",\
                strand=1)
    contigs.qualifiers["note"] = records[0].name #pull out contig number of first contig
    mergingFile.features.append(contigs) #append first contig to mergingFile
    for nextRecord in records[1:]:
        contigs = SeqFeature(FeatureLocation(len(mergingFile), len(mergingFile) + len(nextRecord)), type="fasta_record",\
                strand=1)
        contigs.qualifiers["note"] = nextRecord.name 
        mergingFile.features.append(contigs) #append subsequent contigs to mergingFile
        mergingFile += nextRecord
    mergingFile.name = records[0].name
    mergingFile.description = records[0].description
    mergingFile.annotations = records[0].annotations

    for feature in mergingFile.features:
        if feature.type == 'source':
            mergingFile.features.remove(feature)
    contigs = SeqFeature(FeatureLocation(0, len(mergingFile)), type="source", strand=1)
    mergingFile.features.insert(0,contigs)
    merged_file = re.sub(r"\.\w+$", r".merged.fa", file)
    out_handle = open(merged_file, "w")
    SeqIO.write(mergingFile, out_handle, filetype)
    return merged_file

#END of code adapted from SeqHandler.

def convert2Fasta(file):
    convertedFilename = re.sub(r"\.\w+$", r".fa", file)
    SeqIO.convert(file, determineFileType(file), convertedFilename, 'fasta')


#########################################################################



def error(message):
    sys.stderr.write("\n%s\n\n" % message)
    sys.exit(1)
    
def warning(message):
    sys.stderr.write("\n%s\n\n" % message)


    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
Script to identify the presence of a panel of query genes across a large number of whole genomes. 
  
Requires: blast on your path.

Requires the spaces to be removed after the commas in seqfindr. 
    ''', formatter_class=RawTextHelpFormatter)
        
    
    #takes in the input files
    parser.add_argument('genomes', action="store", help="Directory/Folder containing genomes (fasta format)")

    parser.add_argument("-i", "--input", action="store", required=True, help="Input blast query. E.g. panel of genes formatted for SeqFindr. REQUIRED")
    parser.add_argument("-o", "--output", action="store", required=False, help="Output directory for BLAST results")
    parser.add_argument("-t", "--tol", action="store", default="60", help="TOL cut-off value. Number of aligned bases/total query length")
    parser.add_argument("-f", "--flags", action="store", help="Custom BLAST options, enclosed in quotes. E.g. -f '-task blastn -evalue 0.001'")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Verbose mode")
    parser.add_argument("-p", "--plot", action="store_true", default=False, help="plot heatmap")
    parser.add_argument("-b", "--blast", action="store", default="blastn", choices=("blastn", "tblastx"), help="Blast program to use. Default [blastn]")
    
    args = parser.parse_args()


    main() #run main script




