# Background
Background calculator for Strand-Seq bam files.

## Usage

Download the Background repo:

```bash
git clone https://github.com/YueyanGuo1218/Background.git
```

Install with pip:

```bash
cd ./Background && pip install .
```

Calculate background:

```bash
ssbg /path/to/target.bam --mask /path/to/mask.bed
```

`--mask` is optional.
