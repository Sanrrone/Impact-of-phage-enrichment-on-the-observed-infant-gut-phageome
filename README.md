# Impact of phage enrichment on the observed infant gut phageome

Data and code availability repository for the study:

> **Valenzuela-Diaz S, Dikareva E, Hickman B, Kiljunen S, Kolho K-L, de Vos W, Salonen A, Korpela K.**
> Impact of phage enrichment on the observed infant gut phageome.
> *Microbiology Spectrum* **14**(5):e02153-25 (2026).
> https://doi.org/10.1128/spectrum.02153-25

[![Published in Microbiology Spectrum](https://img.shields.io/badge/Microbiology%20Spectrum-10.1128%2Fspectrum.02153--25-005a9c)](https://journals.asm.org/doi/10.1128/spectrum.02153-25)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15040456.svg)](https://doi.org/10.5281/zenodo.15481575)

<img title="" src="figures/HumanGutVirome%20protocol.png" alt="pipeline" data-align="center">

---

## Overview

This repository provides the sequence data, supplementary tables, and processing
code supporting the paper, which evaluates how a phage-enrichment (VLP) protocol
affects the phageome recovered from infant faecal samples. It contains all mined
and filtered phage contigs, the complete annotation table, the analysis and HPC
pipeline scripts, and the manuscript figures. The repository is intended for data
availability and reproducibility.

Per-sample files are named `S<n>_<timepoint>` (e.g. `S1_6_months`), covering the
1-, 6-, and 12-month timepoints. In total the `contigs/` folder holds **37,450
phage contigs** across **65 per-sample FASTA files**.

---

## Repository contents

### `contigs/`
All phage contigs mined and filtered per sample (`S<n>_<age>_phages.fna`,
nucleotide FASTA).

### `complete_viral_table.zip`
The full per-contig annotation table — all information extracted from the contigs
(taxonomy, host prediction, lifestyle, quality, identity, etc.).

### `r_scripts/`
R scripts used for the downstream analyses and to generate the figures
(annotation parsing, taxonomy, identities, distributions, plotting helpers).

### `slurm_scripts/`
The full processing pipeline as run on an HPC cluster with the SLURM workload
manager (numbered `0_*` → `13_*`: read counting, host removal, assembly, phage
mining/annotation, host prediction, binning, read mapping, and statistics).

### `seqcleaner.py`
Script used to decontaminate the HumGut database (HGDB) using UHGV.

### `figures/`
Manuscript main and supplementary figures (`Figure1`–`Figure6`, `FigureS1`, …)
and the `HumanGutVirome protocol` pipeline schematic.

### Supplementary material
- `supplementary_table1.csv` — per-sample extraction metadata (protocol, age,
  buffer, feces weight, beads, viral DNA concentration).
- `supplementary_table2.csv` — per-metagenome table (1,366 rows).
- `supplementary_table3.xlsx` — additional supplementary data.
- `Supplementary figures.pdf` — compiled supplementary figures.

---

## Citation

If you use these data or code, please cite the paper above. Repository snapshots
are archived on Zenodo (see `CITATION.cff` for the archived DOI and metadata).

## License

Released under the [MIT License](LICENSE).

## Contact

Sandro Valenzuela-Diaz — https://orcid.org/0000-0002-2284-8243
