# Genotype Imputation Pipeline

A cleaned, generalised Bash workflow for genotype-likelihood calling and imputation of human genomic data using **BCFtools**, **GLIMPSE**, **PLINK**, and **GNU Parallel**.

This repository is adapted from a research workflow and has been sanitised for public sharing. Internal server paths, user-specific paths, sample identifiers, and private dataset locations have been replaced with configurable variables or placeholders.

## Workflow

1. Create working directories
2. Generate genotype likelihoods from BAM files with BCFtools
3. Index genotype-likelihood BCF files
4. Run GLIMPSE imputation across autosomes
5. Ligate imputed chunks
6. Merge samples by chromosome
7. Concatenate chromosomes
8. Filter imputed calls at genotype probability (GP) >= 0.99 with PLINK
9. Optionally update FID/IID and sex metadata
10. Optionally convert variant IDs to `CHR:POS`
11. Run the same workflow for chromosome X

## Software

- BCFtools
- GLIMPSE
- PLINK
- GNU Parallel
- Bash
- tabix (for VCF indexing)

You must install/configure these tools separately.

## Reference data

The workflow expects reference files derived from the 1000 Genomes Project and an hg19/GRCh37-compatible reference genome.

Example placeholders used in the script:

```text
/path/to/hs37d5.fa
/path/to/1000G/reference/
```

Update the configuration variables at the top of the script before running.

## Input files

At minimum, the workflow expects:

- A text file containing one BAM path per line
- Indexed BAM files
- 1000 Genomes reference site files
- GLIMPSE reference BCF files
- GLIMPSE genetic maps
- GLIMPSE chunk-definition files

Optional metadata files:

```text
update-ids-glimpse.txt
update-sex.txt
```

`update-ids-glimpse.txt` should contain four columns:

```text
current_FID current_IID new_FID new_IID
```

`update-sex.txt` should contain three columns:

```text
FID IID Sex
```

where PLINK sex codes are `1 = male`, `2 = female`, and `0 = unknown`.

## Usage

Edit the configuration block at the top of:

```text
scripts/imputation_pipeline.sh
```

Then run the required stage, for example:

```bash
bash scripts/imputation_pipeline.sh genotype-likelihoods
bash scripts/imputation_pipeline.sh glimpse
bash scripts/imputation_pipeline.sh postprocess
```

For chromosome X:

```bash
bash scripts/imputation_pipeline.sh chrX
```

## Important notes

- The script is a **portfolio-ready generalisation** of the source workflow and is not intended to run unchanged on another cluster.
- Cluster paths and software locations must be configured for your environment.
- No research data, sample identifiers, or internal server paths are included.
- The source protocol mentions a later Beagle-based “double imputation” step, but the supplied command sections primarily document the GLIMPSE workflow. A Beagle implementation has therefore **not** been invented here.
- Check permissions and ownership before publishing research-derived code.

## Repository structure

```text
genotype-imputation-pipeline/
├── README.md
├── .gitignore
├── config/
│   └── example.env
└── scripts/
    └── imputation_pipeline.sh
```
