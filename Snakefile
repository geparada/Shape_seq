
rule unzip_fastq:
    input:
        "../{rep}/{sample}{rep}_combined_{rd}.fastq.gz"
    output:
        temp("fastq/{sample}{rep}_combined_{rd}.fastq")
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
		

