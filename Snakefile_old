configfile : "config.yaml"

cell = ["Sh3", "Sh5"]
condition = ["DMSO", "NAI"]
rep = ["rep1", "rep2"]
rd = ["R1", "R2"]

adapter3 = {"R1":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "R2":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"}

rule all:
	input:
		expand("{cell}" + "_" + "_".join(condition) + "_ln_nrm.norm.react_overlap." + "_".join(cell) + ".txt_310.15_transcriptome.canonical.fasta_Vienna_package-mfe_sht_0_md_99999", cell=cell),
		expand("react/norm/{cell}" + "_" + "_".join(condition) + "_ln_nrm.norm.react", cell=cell),
		expand("mapped/{cell}-{condition}-{rep}.rtsc",cell=cell,  condition=condition, rep=rep, rd=rd ),
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

rule reverse_mate:
    input:
        "trimmed/sickle/{cell}-{condition}-{rep}_combined_R2.fastq.gz"
    output:
        "trimmed/sickle/{cell}-{condition}-{rep}_combined_R2.rev.fastq.gz"
    conda: "env/core.yaml"		
    shell:
        "python scripts/reverse_2mate.py {input} {output}" 	
    
rule get_non_redudant_transcriptome:
  input:
    config["transcript_fasta"],
    config["canonical_transcripts"]
  conda: "env/core.yaml"
  output:
    "data/transcriptome.canonical.fasta"
  shell:
    "python scripts/get_non_redudant_isoforms.py {input} > {output}"

    
rule bowtie2_build:
  input:
    "data/transcriptome.canonical.fasta"
  output:
    "data/transcriptome.canonical.fasta.1.bt2"
  conda: "env/core.yaml" 
  shell:
    "bowtie2-build {input} {input}"

rule bowtie2:
    input:
        m1="trimmed/sickle/{cell}-{condition}-{rep}_combined_R1.fastq.gz",
        m2="trimmed/sickle/{cell}-{condition}-{rep}_combined_R2.rev.fastq.gz",
        index = "data/transcriptome.canonical.fasta.1.bt2"
    output:
        "mapped/{cell}-{condition}-{rep}.sam"
    params:
        index="data/transcriptome.canonical.fasta",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: 8
    conda:
        "env/core.yaml"
    shell:
        "bowtie2 -x {params.index} -1 {input.m1} -2 {input.m2} -p {threads}  > {output}"

rule sam_filter:
    input:
         "mapped/{cell}-{condition}-{rep}.sam"
    output:
         sam="mapped/{cell}-{condition}-{rep}_filtered.sam",
         log="mapped/{cell}-{condition}-{rep}_filtered.log"
    conda:
        "env/core.yaml"		
    shell:
        "python2 ../StructureFold2/sam_filter.py -sam {input} -logname {output.log}"
        
rule generate_rtsc:
    input:
        "mapped/{cell}-{condition}-{rep}_filtered.sam",
        "data/transcriptome.canonical.fasta"
    output:
        "mapped/{cell}-{condition}-{rep}.rtsc"
    conda:
        "env/core.yaml"
    shell:
        "python2 ../StructureFold2/sam_to_rtsc.py -trim _filtered -single {input} "

        
treatment_dict = dict()

for CELL in cell:
    for CONDITION in condition:
        treatment_dict[(CELL, CONDITION)] = expand("mapped/{cell}-{condition}-{rep}.rtsc",cell=CELL,  condition=CONDITION, rep=rep)
        
        
rule merge_rtsc:
    input:
        lambda w: treatment_dict[(w.cell, w.condition)]
    output:
        "rtsc/{cell}-{condition}.rtsc"
    conda:
        "env/core.yaml"
    shell:
        "python2 ../StructureFold2/rtsc_combine.py {input} -name {output}"

        
rule rtsc_coverage:
    input:
        transcripts = "data/transcriptome.canonical.fasta",   
        rtsc = "rtsc/{cell}-{condition}.rtsc"
    output:
        "coverage/{cell}-{condition}_coverage.csv"
    conda:
        "env/core.yaml"
    shell:
        "python2 ../StructureFold2/rtsc_coverage.py {input.transcripts} -f {input.rtsc} -bases AGCT -name {output}"

        
rule rtsc_react:
    input:
        "rtsc/{cell}-" + condition[0] + ".rtsc",
        "rtsc/{cell}-" + condition[1] + ".rtsc",
        "data/transcriptome.canonical.fasta"
    output:
        "react/{cell}" + "_" + "_".join(condition) + "_ln_nrm.react",
        "react/{cell}" + "_" + "_".join(condition) + "_ln_nrm.scale"
    params:
        "react/{cell}" + "_" + "_".join(condition) + "_ln_nrm"
    conda:
        "env/core.yaml"
    shell:
        "python2 ../StructureFold2/rtsc_to_react.py {input} -name {params} -bases AGCT"                


rule normalise_react:
    input:
        a = "rtsc/{cell}-" + condition[0] + ".rtsc",
        b = "rtsc/{cell}-" + condition[1] + ".rtsc",
        transcripts = "data/transcriptome.canonical.fasta",   
        scale  = "react/" + cell[0] + "_" + "_".join(condition) + "_ln_nrm.scale"
	#scale  = "react/Sh3_DMSO_NAI_ln_nrm.scale"
    output:
        "react/norm/{cell}" + "_" + "_".join(condition) + "_ln_nrm.norm.react"
    params:
        "react/norm/{cell}" + "_" + "_".join(condition) + "_ln_nrm.norm"
    conda:
        "env/core.yaml"
    shell:
        "python2 ../StructureFold2/rtsc_to_react.py {input.a} {input.b} {input.transcripts} -name {params} -bases AGCT -scale {input.scale}" 

rule coverage_overlap:
    input:
        transcripts = "data/transcriptome.canonical.fasta",   
        rtsc = expand("rtsc/{cell}-" + condition[1] + ".rtsc", cell=cell)
    output:
        "coverage/overlap." + "_".join(cell) + ".txt"
    conda:
        "env/core.yaml"
    shell:
        "python2 ../StructureFold2/rtsc_coverage.py {input.transcripts} -f {input.rtsc} -bases AGCT -name {output} -ol -on {output}"


rule fold:
    input:
        RNA_IDs = "coverage/overlap." + "_".join(cell) + ".txt",
        transcripts = "data/transcriptome.canonical.fasta",  	
        react =  "react/norm/{cell}" + "_" + "_".join(condition) + "_ln_nrm.norm.react"
    output:
        directory("{cell}" + "_" + "_".join(condition) + "_ln_nrm.norm.react_overlap." + "_".join(cell) + ".txt_310.15_transcriptome.canonical.fasta_Vienna_package-mfe_sht_0_md_99999")
    conda:
        "env/core.yaml"           
    shell:
        "python2 ../StructureFold2/batch_fold_rna.py {input.RNA_IDs} {input.transcripts} 2 short.txt long.txt -r {input.react}"
