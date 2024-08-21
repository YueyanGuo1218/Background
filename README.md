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

Calculate backgroung for a bam file:

```bash
ssqc /paty/to/target.bam \
  --mask /path/to/mask.bed \
  --output /path/to/output_directory
```

Use `mask` argument to involve masked regions from a bed file to get more accurate result.

The `output` argument specifies the output directory, defaulting to the current working directory.

## Uninstallation

To uninstall the package, run:

```bash
pip uninstall ssqc
```
