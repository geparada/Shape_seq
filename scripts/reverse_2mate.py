from Bio import SeqIO, bgzf, Seq
from io import StringIO
from gzip import open as gzopen
from Bio.Alphabet import generic_dna
import sys

def main( input_fastq, out_fastqgz):
	
	with bgzf.BgzfWriter(out_fastqgz, "wb") as outgz:

		for record in  SeqIO.parse(gzopen(input_fastq, "rt"), format="fastq"):
			new_record = record.seq.reverse_complement()
			new_record.id = record.id
			new_record.description = record.description
			
			SeqIO.write(sequences=new_record, handle=outgz, format="fastq")

		
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
