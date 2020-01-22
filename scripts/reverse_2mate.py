from Bio import SeqIO, bgzf
from io import StringIO
from gzip import open as gzopen

for records in  SeqIO.parse(gzopen("random_10.fastq.gz", "rt"), format="fastq"):

with bgzf.BgzfWriter("test.fastq.bgz", "wb") as outgz:
	SeqIO.write(sequences=records, handle=outgz, format="fastq")
