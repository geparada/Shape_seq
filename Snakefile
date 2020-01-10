
cell = ["Sh3", "Sh5"]
condition = ["DMSO", "NAI"]
rep = ["rep1", "rep2"]
rd = ["R1", "R2"]

adapter3 = {"R1":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "R2":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}


rule all:
	input:
		expand("fastqc/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html", cell=cell, condition=condition, rep=rep, rd=rd )
		#expand("fastq/trim/{cell}-{condition}-{rep}_combined_{rd}.fastq", cell=cell, condition=condition, rep=rep, rd=rd )


rule unzip_fastq:
	input:
		"../{rep}/{cell}-{condition}-{rep}_combined_{rd}.fastq.gz"
	output:
		temp("fastq/{cell}-{condition}-{rep}_combined_{rd}.fastq")
	shell:
		"zcat {input} > {output}"

rule fastqc:
    input:
        "../{rep}/{cell}-{condition}-{rep}_combined_{rd}.fastq.gz"
    output:
        html="fastqc/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html",
        zip="fastqc/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "-a adapters.tab -t 8"
    wrapper:
        "0.47.0/bio/fastqc"

#rule fastqc:
#	input:
#		"../{rep}/{cell}-{condition}-{rep}_combined_{rd}.fastq.gz"
#	output:
#		"fastqc/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html"
#	params:
#		"fastqc/{cell}-{condition}-{rep}_combined_{rd}/"
#	shell:
#		"fastqc {input} -o {params} -a adapters.tab -t 8"		


rule triming:
	input:
		"fastq/{cell}-{condition}-{rep}_combined_{rd}.fastq"
	output:
		temp("fastq/trim/{cell}-{condition}-{rep}_combined_{rd}.fastq")
	params:
		apt3 =  lambda wildcards : adapter3[wildcards.rd],
		minqual = 10
	conda:
		"env/core.yaml"
	shell:
		"python ../StructureFold2/fastq_trimmer.py  {input} -tp {params.apt3} -minqual {params.minqual}"
		

