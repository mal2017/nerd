# nerd

`nerd` is a pipeline for discovery of eRNAs/enhancers from
non-poly-a-selected rna-seq (gro-seq,total rna-seq, etc).

## criteria

bidirectional
located within a predefined search space
unspliced
short transcript
non-polyadenylated

## pipeline

## manifest

## usage

```bash
snakemake -d work/ --configfile config.yaml --cores 6 --use-conda
```

## dependencies
