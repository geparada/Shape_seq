
cell = ["Sh3", "Sh5"]
condition = ["DMSO", "NAI"]
rep = ["rep1", "rep2"]
rd = ["R1", "R2"]

adapter3 = {"R1":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "R2":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"}

rule all:
	input:
		"multiqc/raw_multiqc.html",
		"multiqc/final_multiqc.html"
		#expand("fastqc/sickle/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html", cell=cell, condition=condition, rep=rep, rd=rd ),
		#expand("trimmed/sickle/{cell}-{condition}-{rep}_combined_R1.fastq.gz", cell=cell, condition=condition, rep=rep, rd=rd ),
		#expand("trimmed/{cell}-{condition}-{rep}_combined_{rd}.fastq.gz", cell=cell, condition=condition, rep=rep, rd=rd ),
		#expand("fastqc/raw/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html", cell=cell, condition=condition, rep=rep, rd=rd ),
		#expand("fastqc/cutadapt/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html", cell=cell, condition=condition, rep=rep, rd=rd )
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
        html="fastqc/raw/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html",
        zip="fastqc/raw/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "-a adapters.tab -t 8"
    wrapper:
        "0.47.0/bio/fastqc"


rule cutadapt:
    input: ["../{rep}/{cell}-{condition}-{rep}_combined_R1.fastq.gz", "../{rep}/{cell}-{condition}-{rep}_combined_R2.fastq.gz"]
    output:
        fastq1="trimmed/{cell}-{condition}-{rep}_combined_R1.fastq.gz",
        fastq2="trimmed/{cell}-{condition}-{rep}_combined_R2.fastq.gz",
        qc="trimmed/{cell}-{condition}-{rep}_combined.log"
    params:
        adapters = '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
        others = ""
    threads: 4 # set desired number of threads here
    wrapper:
        "0.47.0/bio/cutadapt/pe"


rule cutadapt_fastqc:
    input:
        "trimmed/{cell}-{condition}-{rep}_combined_{rd}.fastq.gz"
    output:
        html="fastqc/cutadapt/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html",
        zip="fastqc/cutadapt/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "-a adapters.tab -t 8"
    wrapper:
        "0.47.0/bio/fastqc"
	
rule skicle_fastqc:
    input:
        "trimmed/sickle/{cell}-{condition}-{rep}_combined_{rd}.fastq.gz"
    output:
        html="fastqc/sickle/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html",
        zip="fastqc/sickle/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "-a adapters.tab -t 8"
    wrapper:
        "0.47.0/bio/fastqc"
	
rule multiqc_raw:
    input:
        expand("fastqc/raw/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html", cell=cell, condition=condition, rep=rep, rd=rd)
    output:
        "multiqc/raw_multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    wrapper:
        "0.47.0/bio/multiqc"
	
rule multiqc:
    input:
        expand("fastqc/sickle/{cell}-{condition}-{rep}_combined_{rd}/{cell}-{condition}-{rep}_combined_{rd}_fastqc.html", cell=cell, condition=condition, rep=rep, rd=rd),
	expand("trimmed/{cell}-{condition}-{rep}_combined.log", cell=cell, condition=condition, rep=rep, rd=rd)
    output:
        "multiqc/final_multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    wrapper:
        "0.47.0/bio/multiqc"

	
rule sickle_pe:
  input:
    r1="trimmed/{cell}-{condition}-{rep}_combined_R1.fastq.gz",
    r2="trimmed/{cell}-{condition}-{rep}_combined_R2.fastq.gz"
  output:
    r1="trimmed/sickle/{cell}-{condition}-{rep}_combined_R1.fastq.gz",
    r2="trimmed/sickle/{cell}-{condition}-{rep}_combined_R2.fastq.gz",
    rs="trimmed/sickle/{cell}-{condition}-{rep}_combined_single.fastq.gz"
  params:
    qual_type="sanger",
    # optional extra parameters
    extra="-g -q 20 -l 40"
  log:
    # optional log file
    "trimmed/sickle/{cell}-{condition}-{rep}_combined.log"
  wrapper:
    "0.47.0/bio/sickle/pe"

    
rule get_non_redudant_transcriptome:
  input:
    config["transcript_fasta"],
    config["canonical_transcripts"]
  output:
    "data/transcriptome.canonical.fasta"
  shell:
    "python scripts/get_non_redudant_isoforms.py {input} > {output}"

    
rule bowtie2_build:
  input:
    "data/transcriptome.canonical.fasta"
  output:
    "data/transcriptome.canonical.fasta.1.bt2"
  conda: "envs/core.yaml" 
  shell:
    "bowtie2-build {input} {input}"

rule bowtie2:
    input:
        m1="trimmed/sickle/{cell}-{condition}-{rep}_combined_R1.fastq.gz",
        m2="trimmed/sickle/{cell}-{condition}-{rep}_combined_R2.fastq.gz",
        index = "data/transcriptome.canonical.fasta.1.bt2"
    output:
        "mapped/{sample}.bam"
    log:
        "mapped/{sample}.log"
    params:
        index="data/transcriptome.canonical.fasta",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: 8
    conda:
        "envs/core.yaml"
    shell:
        "bowtie2 -x {params.index} -1 {input.m1} -2 {input.m2} -p {threads}   > {output}"

    
