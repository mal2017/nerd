from os.path import realpath
from os.path import split as pathsplit
import subprocess
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import sys

# Block annoying warnings
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# META
__author__ = "Matt Lawlor"

# SETUP
shell.executable("/bin/bash")

# DETERMINE REMOTE OR LOCAL RESOURCE
def determine_resource(path):
    if "gs://" in path:
         return GSRemoteProvider().remote(path.replace("gs://",""))
    elif "ftp://" in path:
         return FTPRemoteProvider().remote(path)
    elif "s3://" in path:
         return S3RemoteProvider().remote(path.replace("s3://",""))
    elif "http://" in path:
         return HTTPRemoteProvider().remote(path.replace("http://",""))
    elif "https://" in path:
         return HTTPRemoteProvider().remote(path.replace("https://",""))
    else:
        return path
# pairings
SAMPLES = list(config["samples"].keys())

ruleorder:  trim_se > trim_pe

rule target:
    input:
        #expand("fastq/{s}_{e}.trimmed.fq.gz", s=SAMPLES, e=["r1","r2"]),
        #expand("idx/genome.{i}.ht2", i=[1,2,3,4,5,6,7,8]),
        #expand("aln/{s}.srt.bam", s=SAMPLES),
        #expand("aln/{s}.srt.bam.bai", s=SAMPLES),
        expand("assembly/{s}/transcripts.gtf", s=SAMPLES)

rule concat_fqs:
    input:
        lambda wc: [determine_resource(x) for x in config["samples"][wc.samp]["fastq"][wc.end]]
    output:
        temp("fastq/{samp}_{end}.fq.gz")
    shell:
        "cat {input} > {output}"

def get_fqs_for_trim(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return ["fastq/{samp}_r1.fq.gz"]
    else:
        return ["fastq/{samp}_r1.fq.gz", "fastq/{samp}_r2.fq.gz"]

def get_fqs_for_aln(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return ["fastq/{samp}_r1.trimmed.fq.gz"]
    else:
        return ["fastq/{samp}_r1.trimmed.fq.gz", "fastq/{samp}_r2.trimmed.fq.gz"]

def get_proper_ended_fastp_call(x):
    fqs = get_fqs_for_trim(x)
    if len(fqs) == 1:
        return "--in1 {r1}".format(r1=fqs[0].format(samp=x))
    else:
        return "--in1 {r1} --in2 {r2}".format(r1=fqs[0].format(samp=x), r2=fqs[1].format(samp=x))

def get_proper_ended_fastp_out(x):
    fqs = get_fqs_for_aln(x)
    if len(fqs) == 1:
        return "--out1 {r1}".format(r1=fqs[0].format(samp=x))
    else:
        return "--out1 {r1} --out2 {r2}".format(r1=fqs[0].format(samp=x), r2=fqs[1].format(samp=x))

rule trim_se:
    input:
        fq = lambda wc: get_fqs_for_trim(wc.samp),
    output:
        r1 = "fastq/{samp}_r1.trimmed.fq.gz",
        html = "fastq/{samp}_fastp.html",
        json = "fastq/{samp}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.samp),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.samp)
    conda:
        "envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.samp}_fastp"

rule trim_pe:
    input:
        fq = lambda wc: get_fqs_for_trim(wc.samp)
    output:
        r1 = "fastq/{samp}_r1.trimmed.fq.gz",
        r2 = "fastq/{samp}_r2.trimmed.fq.gz",
        html = "fastq/{samp}_fastp.html",
        json = "fastq/{samp}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.samp),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.samp)
    conda:
        "envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.samp}_fastp"


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522910/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4502638/


def get_proper_ended_hisat2_call(x):
    fqs = get_fqs_for_aln(x)
    if len(fqs) == 1:
        return "-U {r1}".format(r1=fqs[0].format(samp=x))
    else:
        return "-1 {r1} -2 {r2}".format(r1=fqs[0].format(samp=x), r2=fqs[1].format(samp=x))

rule hisat2_build:
    input:
        fa = determine_resource(config.get("GENOME_FA",None))
    output:
        expand("idx/genome.{i}.ht2",i=[1,2,3,4,5,6,7,8]),
        "idx/genome.fa"
    conda:
        "envs/hisat2.yaml"
    singularity:
        "docker://quay.io/biocontainers/hisat2:2.1.0--py36pl5.22.0_0"
    threads:
        8
    shell:
        "gunzip -c {input.fa} > idx/genome.fa || cp {input.fa} idx/genome.fa; "
        "hisat2-build -p {threads} idx/genome.fa idx/genome"

rule hisat2_aln:
    input:
        #fqs = lambda wc: get_fqs_for_aln(wc.samp),
        fqs = lambda wc: rules.trim_se.output if (len(config["samples"][wc.samp]["fastq"].keys()) == 1) else rules.trim_pe.output,

        idx = rules.hisat2_build.output,
    output:
        temp("aln/{samp}.sam")
    params:
        call = lambda wc: get_proper_ended_hisat2_call(wc.samp),
        lt = config.get("HISAT2_LIB_TYPE","RF"),
        conc = config.get("HISAT2_CONCORDANT_FLAG","--rf")
    conda:
        "envs/hisat2.yaml"
    singularity:
        "docker://quay.io/biocontainers/hisat2:2.1.0--py36pl5.22.0_0"
    threads:
        2
    shell:
        "hisat2 -x idx/genome {params.call} -S {output} -p {threads} "
        "--dta-cufflinks "
        "--rna-strandness {params.lt} --no-unal --no-mixed {params.conc}"

rule prep_alignments:
    input:
        rules.hisat2_aln.output
    output:
        "aln/{samp}.srt.bam"
    threads:
        2
    shell:
        "samtools sort -o {output} -@ {threads} {input}"


rule index_bam:
    input:
        "{file}.srt.bam"
    output:
        "{file}.srt.bam.bai"
    threads:
        2
    conda:
        "envs/samtools19.yaml"
    shell:
        "samtools index -@ {threads} {input}"
rule cufflinks_assemble:
    input:
        bam="aln/{samp}.srt.bam",
        bai="aln/{samp}.srt.bam.bai",
        gtf=determine_resource(config.get("GENOME_GTF",None))
    output:
        "assembly/{samp}/genes.fpkm_tracking",
        "assembly/{samp}/isoforms.fpkm_tracking",
        "assembly/{samp}/skipped.gtf",
        "assembly/{samp}/transcripts.gtf",
    threads:
        2
    conda:
        "envs/cufflinks.yaml"
    singularity:
        "docker://quay.io/biocontainers/cufflinks:2.2.1--py36_2"
    shell:
        "cufflinks -g {input.gtf} -o assembly/{wildcards.samp} -p {threads} {input.bam}"

# https://www.biostars.org/p/271203/
