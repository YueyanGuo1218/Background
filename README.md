# Background
Background calculator for Strand-Seq bam files.

## Installation

Download the Background repo:

```bash
git clone https://github.com/YueyanGuo1218/Background.git
```

Install with pip:

```bash
cd ./Background && pip install .
```

## Usage

Calculate backgroung for a bam file:

```bash
ssbg /paty/to/target.bam \
  --mask /path/to/mask.bed \
  --output /path/to/output_directory
```

Use `mask` argument to involve masked regions from a bed file to get more accurate result.

The `output` argument specifies the output directory, defaulting to the current working directory.

## Uninstallation

To uninstall the package, run:

```bash
pip uninstall Background
```
