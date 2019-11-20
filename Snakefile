
cell = ["Sh3", "Sh5"]
condition = ["DMSO", "NAI"]
rep = ["rep1", "rep2"]
rd = ["R1", "R2"]


rule all:
	input:
		expand("fastq/{cell}-{condition}-{rep}_combined_{rd}.fastq", cell=cell, condition=condition, rep=rep, rd=rd )


rule unzip_fastq:
	input:
		"../{rep}/{cell}-{condition}-{rep}_combined_{rd}.fastq.gz"
	output:
		temp("fastq/{cell}-{condition}-{rep}_combined_{rd}.fastq")
	shell:
		"zcat {input} > {output}"



rule triming:
	input:
		"../{rep}/{sample}{rd}.fastq.gz"
	output:
		"fastq/{sample}{rd}.fastq.trim"
	conda:
		"../envs/core.yaml"
	shell:
		"bwa index {input}"
		

