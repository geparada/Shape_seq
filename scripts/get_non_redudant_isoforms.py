import sys
import csv
from Bio import SeqIO

def main(fasta, knownCanonical):
    
    main_transcripts = set([])
    
    with open(knownCanonical) as knownCanonical_file:
        reader = csv.DictReader(knownCanonical_file, delimiter="\t")
        
        for row in reader:
            main_transcripts.add(row["transcript"])
            
    with open(fasta) as fasta_file:
        
        for t in  SeqIO.parse(fasta_file, "fasta"):
            
            ID =  t.id.split("|")[0]
            
            if ID in main_transcripts:
                
                print(">" + ID)
                print(t.seq)
            


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

