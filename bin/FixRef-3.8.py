#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
from contextlib import nullcontext
from pyfaidx import Fasta
import os
import re

def get_args():
    parser = argparse.ArgumentParser(description='The program fixes the issue with VCF file that reference allele is in ALT column instead of REF.')
    parser.add_argument("-f", "--reference", help="reference sequence in fasta format ", required=True, dest="ref")
    parser.add_argument("-v", "--vcf", help="input vcf file", required=True)
    parser.add_argument("-o", "--output", help="output vcf file", required=True)
    return parser.parse_args()

def processVCF(infile, outfile, refile):
    genome = Fasta(refile)
    nRef = nFlip = nErr = 0
    input = open(infile)
    output = open(outfile, "w")
    nread = 0
    for line in input.readlines():
        lines = line.split("\t")
        if line.startswith("#"):
            output.write(line)
        else:
            Chr = lines[0]
            REF = lines[3]
            ALTs = lines[4].split(",")            
            Start = int(lines[1])-1
            End1 = int(Start) + len(REF)
            refSeq = genome[Chr][Start:End1].seq 

            NoRef = True #none of the allele can be found in reference sequence
            if REF == refSeq:
                nRef = nRef+1
                output.write(line)
                NoRef = False
                continue
            else:
                for i,allele in enumerate(ALTs):
                    End2 = int(Start) + len(allele)
                    if End1 != End2:
                        refSeq = genome[Chr][Start:End1].seq 
                    if allele == refSeq:
                        NoRef = False
                        nFlip = nFlip+1
                        lines[3] = refSeq
                        ALTs.remove(allele)
                        ALTs.insert(i, REF)
                        lines[4] = ",".join(ALTs)
                        for j in range(9, len(lines)): #genotype for sample j
                            if lines[j].rstrip() == ".":
                                continue
                            Geno = lines[j].rstrip().split(":")
                            GT = Geno[0]
                            if re.search("/", GT):
                                GTs = GT.split("/")
                                tmpGTs = GTs
                                for t in range(0,2):
                                    if tmpGTs[t] == '0':
                                        GTs[t] = str(i+1)
                                    elif tmpGTs[t] == str(i+1):
                                        GTs[t] = '0'
                                if int(GTs[0]) > int(GTs[1]):
                                    GTs = GTs[::-1]
                                GT = "/".join(GTs)
                                Geno[0] = GT
                                lines[j] = ":".join(Geno)
                            elif re.search("\|", GT):
                                GTs = GT.split("|")
                                tmpGTs = GTs
                                for t in range(0,2):
                                    if tmpGTs[t] == '0':
                                        GTs[t] = str(i+1)
                                    elif tmpGTs[t] == str(i+1):
                                        GTs[t] = '0'
                                if int(GTs[0]) < int(GTs[1]):
                                    GTs = GTs[::-1]
                                GT = "|".join(GTs)
                                Geno[0] = GT
                                lines[j] = ":".join(Geno)
                            else:
                                print ("Error: Non-missing genotype not seperated by '/' nor '|' at %s:%s") % (Chr, Start)
                        line = "\t".join(lines).rstrip()+"\n"
                        output.write(line)
                        continue #stop once found match allele
            if NoRef:
                nErr = nErr+1
                print ("Reference allele %s not found for %s:%s [%s] [%s]" % (refSeq, lines[0], End1, lines[3], lines[4])) 
        
    #nread = nread+1
    #if nread%1000 == 0:
    #    print ("********** finished %s variants, %s:%s **********" % (nread, # Chr, Start))  // This might slow down the process
            
    print ("\n*****************************************")
    print ("* Number of variants have matched allele: %s" % (nRef))
    print ("* Number of variants have flipped allele: %s" % (nFlip))
    print ("* Number of variants have no-match allele: %s (dropped)" % (nErr)) 
    print ("*****************************************\n")
    input.close()
    output.close()

if __name__=="__main__":
    args = get_args()
    print ("Please note that this software can only be used to flip REF and ALT alleles and corresponding GT calling, to make the REF allele same as reference genome.\nMake sure you have the same version of reference for variant calling.")
    processVCF(args.vcf, args.output, args.ref)
    
    
