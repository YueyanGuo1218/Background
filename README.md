# Strand-seq Library QC
Quality control for Strand-Seq bam files based on measures:
- Background
- Depth
- Spikiness
- Half-Depth
- Entropy

## Installation

Download the Background repo:

```bash
git clone https://github.com/YueyanGuo1218/ssqc.git
```

Install with pip:

```bash
cd ./ssqc && pip install .
```

## Usage

### Generate reports

Calculate backgroung for a bam file:

```bash
ssqc /paty/to/target.bam \
  --mask /path/to/mask.bed \
  --output /path/to/output_directory
```

Use `mask` argument to involve masked regions from a bed file to get more accurate result.

The `output` argument specifies the output directory, defaulting to the current working directory.

### Merge reports

Merge all TSV files in given directory:

```bash
ssqc-merge /path/to/tsv_directory \
  --output /path/to/output.tsv
```

In most cases, the input directory should be the same as the output directory of the `ssqc` cammand.

Clustering will be done and the result will be attached automaticly.

Note that putting the output TSV file in the input directory results in errors when running `ssqc-merge` command on the same input directory again.

### Parallelization

An example of using GNU Parallel:

```bash
ls /path/to/*.bam | parallel ssqc {} --output /path/to/output
```

## Uninstallation

To uninstall the package, run:

```bash
pip uninstall ssqc
```
