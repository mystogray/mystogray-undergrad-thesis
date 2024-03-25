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
Based on alpha diversity analysis, there is a significant difference (p<0.05) between the two groups, with group A having a higher Shannon index in comparison to group B. This indicates that the bacterial community in individuals within group A is **more** even and diverse compared to group B.
![plot](/.assets/porto_alpha.svg)

### How Does The Bacterial Community Differ from One Another?
Beta diversity analysis using Bray-Curtis index shows that each group is making their own clusters on different PCoA axises. However, there is one outlier on each group, which is A1 and B3. This pattern can be explained by looking at the genus diversity of each sample.\
![plot](/.assets/porto_betadiv.svg) \
Looking at the barplot below, it can be seen that most of the individuals in group B is dominated by _Burkholderia_ genus, except for B3. Meanwhile on group A, there is no predominant genus on all samples. If we look closely, sample B3 does not have a dominant genus, similar to group A. This may explain the findings before that B3 is in the same clusters as A2 and A3. \
![plot](/.assets/porto_genusabundance.svg)

