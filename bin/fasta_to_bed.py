import sys
import gzip
import re

def fasta_to_bed(input_file, output_file):
    fasta = gzip.open(input_file, "rt") if input_file.endswith(".gz") else open(input_file)
    bed = open(output_file, "w")
    start_pos = 0
    chromosome = None
    regex = re.compile(r">(.*)")
    chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                   "11", "12", "13", "14", "15", "16", "17", "18",
                   "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "X", "Y"]
    sequence = ""
    for line in fasta:
        line = line.strip()
        if line.startswith(">"):
            if chromosome and chromosome in chromosomes:
                end_pos = start_pos + len(sequence)
                bed.write("\t".join([chromosome, str(start_pos), str(end_pos)]) + "\n")
            captures = regex.search(line)
            chromosome = captures.group(1)
            if chromosome not in chromosomes:
                chromosome = None
                continue
            start_pos = 0
            sequence = ""
        else:
            sequence += line
    if chromosome and chromosome in chromosomes:
        end_pos = start_pos + len(sequence)
        bed.write("\t".join([chromosome, str(start_pos), str(end_pos)]) + "\n")
    fasta.close()
    bed.close()

input_file = sys.argv[1]
output_file = sys.argv[2]
fasta_to_bed(input_file, output_file)