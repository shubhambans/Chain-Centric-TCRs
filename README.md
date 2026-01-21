Prevalent α-Chain–Centric TCR Metaclones Shape Specificity and Focus Early SARS-CoV-2 T cell Responses

## Repository layout

- `Figure_code/`: R scripts and notes organized by figure and dataset (see below).
- `Figure_PDFs/`: exported figure artifacts (PDFs or other outputs).
- `Supplementary_tables/`: supplementary table notes and helper scripts for supplementary figures.
- `data/`: placeholder for raw/processed datasets used by the scripts.

## TCR dataset references

Set | Source | Description
--- | --- | ---
Set 1 | VDJdb (VDJdb 26) | Human antigen-specific TCRs.
Set 2 | Bacher et al. 23 | Enriched for SARS-CoV-2 antigen-specific T cells.
Set 3 | Su et al. 23 | HC individuals.
Set 4 | Su et al. 23 | COVID-19 + HC individuals.
Set 5 | Combined 23–26 | Combined COVID-19 + HC individuals.

## Figure code organization

- `Figure_code/Figure_1/`: Figure 1 scripts for chain-centric frequencies and metaclone properties (Sets 1 & 3).
- `Figure_code/Figure_2/`: Figure 2 scripts for multi-specificity overlap analyses (Set 1).
- `Figure_code/Figure_3/`: Figure 3 scripts for COVID-19/publicity analyses and motif plots (Sets 4 & 5).
- `Figure_code/Figure_4/`: Figure 4 scripts for T cell subset analysis and visualization (Set 2).
- `Figure_code/Figure_5/`: Figure 5 scripts for pseudotime and mismatch analyses (Set 2).

## Notes

Many figure scripts are stored without file extensions. Open them in an editor as R scripts when needed.
