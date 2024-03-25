# Metagenomic Analysis Between Individuals Within Two Groups
Undergraduate thesis of Vincentius Suryo S. \
Disclaimer: Some information cannot be given until this paper is published in a journal.

## General Objectives
The aim of the project is to investigate whether there is any difference of **bacterial** abundance between two groups.

## Methodology
```mermaid
flowchart LR
    A["QC and Data Wrangling
    (CutAdapt & FastQC)"] --> B["Trimming
    (CutAdapt)"]
    B --> C["Dehosting
    (Bowtie2)"]
    C --> D["Taxonomy Identification
    (Kraken2)"]
    D --> E["Functional Analysis
    (HUMAnN3)"]
    E --> F["Statistical Evaluation
    and Visualisation"]
```

## Result
### How Even and Diverse Is The Bacterial Community in Each Individual?
!(/.assets/porto_alpha.svg)

### How Does The Bacterial Community Differ from One Another?
(beta div)
(genus)
