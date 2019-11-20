rule bwa_index:
	input:
		"../{rep}/{sample}.fastq.gz"
	output:
		"Round1/ME_TAGs.fa.amb"
	conda:
		"../envs/core.yaml"
	shell:
		"bwa index {input}"
