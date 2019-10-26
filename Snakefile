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

configfile: "config.yaml"

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

"""
criteria:
    bidirectional
    located within a search space
    unspliced
    short transcript
    non-polyadenylated
"""




rule target:
    input:
        expand("fastq/{s}.trimmed.fq.gz", s=SAMPLES)


rule concat_fqs:




# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522910/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4502638/
rule trim:
    input:
        []
