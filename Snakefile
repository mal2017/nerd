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

rule target:
    input:
        expand("fastq/{s}_{e}.trimmed.fq.gz", s=SAMPLES, e=["r1","r2"]),
        expand("aln/{s}", s=SAMPLES)


rule concat_fqs:
    input:
        lambda wc: config["samples"][wc.samp]["fastq"][wc.end]
    output:
        temp("fastq/{samp}_{end}.fq.gz")
    shell:
        "cat {input} > {output}"

rule trim:
    input:
        r1 = "fastq/{samp}_r1.fq.gz",
        r2 = "fastq/{samp}_r2.fq.gz"
    output:
        r1 = "fastq/{samp}_r1.trimmed.fq.gz",
        r2 = "fastq/{samp}_r2.trimmed.fq.gz",
        html = "fastq/{samp}_fastp.html",
        json = "fastq/{samp}_fastp.json"
    threads:
        2
    conda:
        "envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.samp}_fastp"

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522910/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4502638/

def get_fqs_for_aln(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return "fastq/{samp}_r1.fq.gz"
    else:
        return ["fastq/{samp}_r1.fq.gz", "fastq/{samp}_r2.fq.gz"]

def get_proper_ended_tophat2_call(x):
    fqs = get_fqs_for_aln(x)
    if len(x) == 1:
        return "-U {r1}".format(r1=fqs[0])
    else:
        return "{r1} {r2}".format(r1=fqs[0].format(samp=x), r2=fqs[1].format(samp=x))

rule align:
    input:
        fqs = lambda wc: get_fqs_for_aln(wc.samp),
        bt2_idx_paths = [determine_resource(y) for y in config.get("BT2_FILES",None)]
    output:
        directory("aln/{samp}")
    singularity:
        "docker://quay.io/biocontainers/tophat:2.1.1--py35_0"
    conda:
        "envs/tophat2.yaml"
    params:
        call = lambda wc: get_proper_ended_tophat2_call(wc.samp),
        idx_pfx = config.get("BT2_IDX_PFX",None),
        lt = config.get("TOPHAT2_LIB_TYPE","fr-firststrand")
    threads:
        4
    shell:
        "tophat2 -p {threads} "
        "--library-type {params.lt} "
        "--no-discordant --no-mixed "
        "--b2-very-fast -o {output} "
        "{params.idx_pfx} {params.call}"
