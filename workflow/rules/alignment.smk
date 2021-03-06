## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

rule clean_alignment:
    shell:
        "rm -fr results/alignment/"

rule run_alignment:
    input:
        expand("results/alignment/bowtie2/{sample}_sorted.bam.bai", sample = samples(pep))

def get_fnas_for_bt2_index(config, reference):
    return config["reference"][reference]["fastas"]

def get_bt2_index(sample, pep):
    reference = lookup_sample_metadata(sample, "genome", pep)
    return "results/alignment/bowtie2_index/%s/%s"%(reference,reference)

def get_bt2_index_file(sample, pep):
    reference = lookup_sample_metadata(sample, "genome", pep)
    return "results/alignment/bowtie2_index/%s/%s.1.bt2"%(reference,reference)
    

rule pull_genbank:
    message: "Download genbank for {wildcards.accession}"
    output:
        "resources/genbanks/{accession}.gbk"
    log:
        stdout="results/alignment/logs/pull_genbank/{accession}.log",
        stderr="results/alignment/logs/pull_genbank/{accession}.err"
    threads: 1
    conda:
        "../envs/alignment.yaml"
    shell:
        "ncbi-acc-download {wildcards.accession} --out {output} > {log.stdout} "
        "2> {log.stderr}"

rule process_genbank:
    message: "Processing genbank for {wildcards.genome}"
    input:
        "resources/genbanks/{genome}.gbk"
    output:
        "results/alignment/process_genbank/{genome}/{genome}.{outfmt}"
    log:
        stdout="results/alignment/logs/process_genbank/{genome}_{outfmt}.log",
        stderr="results/alignment/logs/process_genbank/{genome}_{outfmt}.err"
    threads: 1
    conda:
        "../envs/alignment.yaml"
    shell:
         "python3 workflow/scripts/parse_genbank.py {input} "
         "--outfmt {wildcards.outfmt} "
         "--chrm '{wildcards.genome}'  "
         " > {output} 2> {log.stderr}"

rule bowtie2_index:
    input:
        lambda wildcards: get_fnas_for_bt2_index(config, wildcards.reference)
    output:
        "results/alignment/bowtie2_index/{reference}/{reference}.1.bt2"
    params:
        fastas = lambda wildcards, input: ",".join(input)
    threads:
        5
    log:
        stdout="results/alignment/logs/bowtie2_index/{reference}.log",
        stderr="results/alignment/logs/bowtie2_index/{reference}.err" 
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2-build --threads {threads} "
        "{params.fastas} "
        "results/alignment/bowtie2_index/{wildcards.reference}/{wildcards.reference} "
        "> {log.stdout} 2> {log.stderr}"

rule bowtie2_map:
    input:
        in1="results/preprocessing/trimmomatic/{sample}_trim_paired_R1.fastq.gz",
        in2="results/preprocessing/trimmomatic/{sample}_trim_paired_R2.fastq.gz",
        bt2_index= lambda wildcards: get_bt2_index_file(wildcards.sample,pep)
    output:
        temp("results/alignment/bowtie2/{sample}.bam")
    log:
        stderr="results/alignment/logs/bowtie2/{sample}_bt2.log" 
    params:
        bt2_index= lambda wildcards: get_bt2_index(wildcards.sample,pep)    
    threads: 
        5
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2 -x {params.bt2_index} -p {threads} "
        "-1 {input.in1} -2 {input.in2} --phred33  "
        "--end-to-end --very-sensitive 2> {log.stderr} "
        "| samtools view -b > {output}"

rule bam_sort:
    input:
        "results/alignment/bowtie2/{sample}.bam"
    output:
        "results/alignment/bowtie2/{sample}_sorted.bam"
    log:
        stderr="results/alignment/logs/bowtie2/{sample}_bt2_sort.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools sort {input} > {output} 2> {log.stderr}"

rule bam_index:
    input:
        "results/alignment/bowtie2/{sample}_sorted.bam"
    output:
        "results/alignment/bowtie2/{sample}_sorted.bam.bai"
    log:
        stdout="results/alignment/logs/bowtie2/{sample}_bt2_index.log",
        stderr="results/alignment/logs/bowtie2/{sample}_bt2_index.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools index {input} {output} > {log.stdout} 2> {log.stderr}"
