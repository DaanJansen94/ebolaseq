# Methodology

This document describes how ebolaseq produces the optional GP mAb escape report: from GP alignment through epitope scoring to the files written in `Escape/`.

## Pipeline

1. Enable report → force GP protein alignment.
2. Extract GP CDS, translate, MAFFT → `Alignment/MAFFT/GP/protein_aln.fasta`.
3. Map epitope positions (Makona 2014 frame) to MSA columns.
4. Score isolates vs Makona; summarise per therapy (Ebanga / Inmazeb / MBP134).
5. Write `Escape/` (HTML, Excel, verification files).

Makona 2014 (KJ660346.2) is the reference frame because the licensed mAbs were tested against this outbreak strain (watchlist sites match the identical UniProt Q05320 / Mayinga-76 numbering).

## GP alignment (step 2)

Extract GP CDS, translate, then MAFFT all proteins into one alignment. Gaps in the MSA can shift position numbers (see step 3).

## Position mapping (step 3)

Each literature epitope site has a Makona residue number (e.g. 116). Because MAFFT adds gap columns, each number is mapped to the correct MSA column via an anchor sequence. Output: `gp_epitope_q05320_msa_map.tsv`.

## Scoring (step 4)

Per isolate at each mapped site vs Makona `zaire_aa`:

| Status | Meaning |
|--------|---------|
| `same` | matches Makona |
| `escape` | matches a catalogued escape mutation |
| `changed` | different, not in catalogue |
| `gap` | missing or unmapped |

Grantham distance for substitutions. Per therapy: contact-site conservation, epitope-change count, escape matches, max Grantham. In silico concern flags per cocktail component. Catalogue: [MAB_ESCAPE_EPITOPES.md](MAB_ESCAPE_EPITOPES.md#catalogued-escape-mutations).

## Outputs (step 5)

| File | Content |
|------|---------|
| `gp_mab_escape_report.html` | Interactive report |
| `mab_escape_data.xlsx` | Tables for R/Excel |
| `gp_protein_aln.fasta` | GP MSA |
| `makona_gp_mature_reference.fasta` | Makona baseline |
| `gp_epitope_q05320_msa_map.tsv` | Position ↔ column map |

## Limitations

- In silico only; catalogue match does not prove reduced clinical efficacy.
- Only watchlist positions in [MAB_ESCAPE_EPITOPES.md](MAB_ESCAPE_EPITOPES.md); other changes appear as `changed`.
- MAFFT gaps and partial GP can mis-assign columns — verify with `gp_epitope_q05320_msa_map.tsv`.
- All isolates compared to Makona 2014 `zaire_aa`, including non-Zaire species, for consistent numbering.
