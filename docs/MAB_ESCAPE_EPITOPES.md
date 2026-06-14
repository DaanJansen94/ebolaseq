# mAb epitope positions (ebolaseq mab-escape-report)

Literature GP epitopes for **Ebanga (mAb114)**, **Inmazeb (REGN-EB3)**, and **MBP134 (ADI-15878 + ADI-23774)**, mapped to **UniProt Q05320** positions and **2014 Makona** amino acids ([KJ660346.2](https://www.ncbi.nlm.nih.gov/nuccore/KJ660346.2)). The report watchlist is [`gp_epitope_watchlist.json`](../ebolaseq/data/gp_epitope_watchlist.json); tables below match that file.

Pipeline steps (alignment → scoring → outputs): [Methodology](Methodology.md).

## Therapies covered

| Product | Monoclonal antibody(ies) |
|---------|--------------------------|
| Ebanga (ansuvimab) | mAb114 |
| Inmazeb (REGN-EB3) | mAb3470 (atoltivimab), mAb3471 (odesivimab), mAb3479 (maftivimab) |
| MBP134 / MBP134AF | ADI-15878, ADI-23774 |

Details: [Coordinate system](#coordinate-system).

## mAb114 (Ebanga / ansuvimab)

GP1 receptor-binding head and glycan cap.

| Pos (Q05320) | Makona 2014 aa | GP subunit | Region | Proven escape mutants | Ref |
|-------------:|---------:|:----------:|--------|----------------------|:---:|
| 112 | E | GP1 | GP1 head / RBS beta7-beta9 | — | 1 |
| 116 | P | GP1 | GP1 head / RBS beta7-beta9 | P116H | 1 |
| 118 | G | GP1 | GP1 head / RBS beta7-beta9 | G118E, G118R | 1 |
| 143 | G | GP1 | GP1 head / RBS beta7-beta9 | — | 1 |
| 144 | T | GP1 | GP1 head / RBS beta7-beta9 | T144M | 1 |
| 146 | P | GP1 | GP1 head / RBS beta7-beta9 | — | 1 |
| 270 | T | GP1 | Glycan cap | — | 1 |

1. Misasi et al., Science 2016 — epitope ([10.1126/science.aad6117](https://doi.org/10.1126/science.aad6117); [PDB 5FHC](https://www.rcsb.org/structure/5FHC))

## REGN3470 — atoltivimab (Inmazeb component 1)

GP1 glycan cap.

| Pos (Q05320) | Makona 2014 aa | GP subunit | Region | Proven escape mutants | Ref |
|-------------:|---------:|:----------:|--------|----------------------|:---:|
| 278 | N | GP1 | Glycan cap beta17-beta18 | — | 1 |
| 279 | P | GP1 | Glycan cap loop 279-306 | — | 1 |
| 280 | E | GP1 | Glycan cap beta17-beta18 | — | 1 |
| 282 | D | GP1 | Glycan cap beta17-beta18 | — | 1 |
| 296 | N | GP1 | Glycan cap loop 279-306 | — | 1 |

1. Rayaprolu et al., Cell Host Microbe 2023 — REGN-EB3 contacts and escape ([10.1016/j.chom.2023.01.002](https://doi.org/10.1016/j.chom.2023.01.002); [EMDB-26005](https://www.ebi.ac.uk/emdb/EMD-26005))

## REGN3471 — odesivimab (Inmazeb component 2)

GP1 head + glycan cap.

| Pos (Q05320) | Makona 2014 aa | GP subunit | Region | Proven escape mutants | Ref |
|-------------:|---------:|:----------:|--------|----------------------|:---:|
| 111 | L | GP1 | GP1 head / RBS beta7-beta9 | — | 1 |
| 112 | E | GP1 | GP1 head / RBS beta7-beta9 | — | 1 |
| 116 | P | GP1 | GP1 head / RBS beta7-beta9 | P116H | 1 |
| 118 | G | GP1 | GP1 head / RBS beta7-beta9 | G118E, G118R | 1 |
| 143 | G | GP1 | GP1 head / RBS beta7-beta9 | — | 1 |
| 144 | T | GP1 | GP1 head / RBS beta7-beta9 | T144M | 1 |
| 146 | P | GP1 | GP1 head / RBS beta7-beta9 | — | 1 |
| 270 | T | GP1 | Glycan cap | — | 1 |
| 272 | K | GP1 | Glycan cap | — | 1 |
| 274 | I | GP1 | Glycan cap | I274M | 1 |
| 275 | W | GP1 | Glycan cap | W275L | 1 |

1. Rayaprolu et al., Cell Host Microbe 2023 ([10.1016/j.chom.2023.01.002](https://doi.org/10.1016/j.chom.2023.01.002))
2. Davidson et al., J Virol 2015 — I274M, W275L ([10.1128/JVI.01490-15](https://doi.org/10.1128/JVI.01490-15))

## REGN3479 — maftivimab (Inmazeb component 3)

GP base / fusion loop (quaternary).

| Pos (Q05320) | Makona 2014 aa | GP subunit | Region | Proven escape mutants | Ref |
|-------------:|---------:|:----------:|--------|----------------------|:---:|
| 34 | P | GP1 | GP1 quaternary 2nd protomer | — | 1 |
| 43 | L | GP1 | GP1 quaternary 2nd protomer | — | 1 |
| 45 | V | GP1 | GP1 quaternary 2nd protomer | — | 1 |
| 504 | I | GP2 | Fusion loop / GP base | — | 1 |
| 505 | V | GP2 | Fusion loop / GP base | — | 1 |
| 506 | N | GP2 | Fusion loop / GP base | — | 1 |
| 507 | A | GP2 | Fusion loop / GP base | — | 1 |
| 527 | I | GP2 | Fusion loop / GP base | — | 1 |
| 528 | G | GP2 | Fusion loop / GP base | G528R | 1 |
| 529 | L | GP2 | Fusion loop / GP base | — | 1 |
| 530 | A | GP2 | Fusion loop / GP base | — | 1 |
| 535 | F | GP2 | Fusion loop / GP base | — | 1 |
| 536 | G | GP2 | Fusion loop / GP base | — | 1 |
| 560 | Q | GP2 | HR1 | — | 1 |
| 563 | N | GP2 | HR1 N-linked glycan | N563Y | 1 |
| 564 | E | GP2 | HR1 | — | 1 |
| 567 | Q | GP2 | HR1 | — | 1 |

1. Rayaprolu et al., Cell Host Microbe 2023 ([10.1016/j.chom.2023.01.002](https://doi.org/10.1016/j.chom.2023.01.002))

## ADI-15878 (MBP134 component 1)

Pan-ebolavirus GP-base mAb; quaternary IFL / HR1 epitope spanning two protomers (overlaps REGN3479 IFL; escape **G528E/S** not **G528R**). Binds the conserved **N-term pocket** beneath the GP2 N-terminal tail (502–510); tail residues **504–507 are displaced and disordered** in the ADI-15878–GP complex and are **not** contact residues (2).

| Pos (Q05320) | Makona 2014 aa | GP subunit | Region | Proven escape mutants | Ref |
|-------------:|---------:|:----------:|--------|----------------------|:---:|
| 34 | P | GP1 | GP1 base, 2nd protomer | — | 1 |
| 45 | V | GP1 | GP1 base, 2nd protomer | — | 1 |
| 47 | D | GP1 | GP1 base, 2nd protomer | — | 1 |
| 527 | I | GP2 | IFL, adjacent protomer | — | 1 |
| 528 | G | GP2 | IFL, adjacent protomer | G528E, G528S | 1, 3 |
| 529 | L | GP2 | IFL, adjacent protomer | — | 1 |
| 530 | A | GP2 | IFL, adjacent protomer | — | 1 |
| 535 | F | GP2 | IFL, adjacent protomer | — | 1 |
| 536 | G | GP2 | IFL, adjacent protomer | — | 1 |
| 537 | P | GP2 | IFL, adjacent protomer | — | 1 |
| 559 | R | GP2 | HR1, adjacent protomer | — | 1 |
| 560 | Q | GP2 | HR1, adjacent protomer | — | 1 |
| 561 | L | GP2 | HR1, adjacent protomer | — | 1 |
| 563 | N | GP2 | HR1 glycan (NAG563) | N563Y | 1 |
| 564 | E | GP2 | HR1, adjacent protomer | — | 1 |
| 567 | Q | GP2 | HR1, adjacent protomer | — | 1 |

1. Murin et al., Cell Rep 2018 — GP–ADI-15878 contacts ≤4 Å, Table S2; EBOV/Mak cryo-EM ([10.1016/j.celrep.2018.08.009](https://doi.org/10.1016/j.celrep.2018.08.009); [EMD-8935](https://www.ebi.ac.uk/emdb/EMD-8935))
2. West et al., mBio 2018 — GPCL crystal; N-term tail displaced, pocket binding ([10.1128/mBio.01674-18](https://doi.org/10.1128/mBio.01674-18); [PDB 6EA7](https://www.rcsb.org/structure/6EA7))
3. Bornholdt et al., Cell 2017 — G528E/S escape ([10.1016/j.cell.2017.04.037](https://doi.org/10.1016/j.cell.2017.04.037))

## ADI-23774 (MBP134 component 2)

310-pocket mAb. **Same GP epitope as parental ADI-15946.** **ADI-23774** is affinity-matured from ADI-15946 (yeast display) to improve SUDV GP binding—not a new GP footprint. MBP134 = **ADI-15878** + **ADI-23774** (non-overlapping sites).

**K510:** HC contacts **K510**; **K510E** abolishes binding/neutralization ([escape catalogue](#catalogued-escape-mutations)).

| Pos (Q05320) | Makona 2014 aa | GP subunit | Region | Proven escape mutants | Ref |
|-------------:|---------:|:----------:|--------|----------------------|:---:|
| 71 | E | GP1 | GP1 base / 310 helix (310 pocket) | — | 1 |
| 72 | G | GP1 | GP1 base / 310 helix (310 pocket) | — | 1 |
| 73 | N | GP1 | GP1 base / 310 helix (310 pocket) | — | 1 |
| 74 | G | GP1 | GP1 base / 310 helix (310 pocket) | — | 1 |
| 75 | V | GP1 | GP1 base / 310 helix (310 pocket) | — | 1 |
| 287 | E | GP1 | Glycan cap loop 279-306 | — | 1 |
| 288 | W | GP1 | Glycan cap loop 279-306 | — | 1 |
| 289 | A | GP1 | Glycan cap loop 279-306 | — | 1 |
| 290 | F | GP1 | Glycan cap loop 279-306 | — | 1 |
| 291 | W | GP1 | Glycan cap loop 279-306 | W291R | 1 |
| 292 | E | GP1 | Glycan cap loop 279-306 | — | 1 |
| 510 | K | GP2 | GP2 N-terminus / 310 pocket interface | K510E | 1 |

1. Bornholdt et al., Nat Struct Mol Biol 2019 — 310 pocket, GPCL contacts, **K510** / **K510E** ([10.1038/s41594-019-0191-4](https://doi.org/10.1038/s41594-019-0191-4); [PDB 6MAM](https://www.rcsb.org/structure/6MAM))
2. Bornholdt et al., Cell 2017 — ADI-15946 discovery ([10.1016/j.cell.2017.04.037](https://doi.org/10.1016/j.cell.2017.04.037))
3. Wec et al., Cell Host Microbe 2019 — **MBP134** cocktail ([10.1016/j.chom.2018.12.004](https://doi.org/10.1016/j.chom.2018.12.004))
4. Bornholdt et al., Cell Host Microbe 2019 — **MBP134AF** ([10.1016/j.chom.2018.12.005](https://doi.org/10.1016/j.chom.2018.12.005))

*Refs 3–4: product names for the two-mAb mix, not extra epitope sites. Positions 287–292: full-length GP / GPFL competition (contact=false in watchlist).*

## Catalogued escape mutations

| Mutation | Pos | Change | GP subunit | Linked mAb | Product(s) | Region | Ref |
|----------|----:|--------|:----------:|------------|------------|--------|:---:|
| P116H | 116 | P→H | GP1 | mAb114, REGN3471 | Ebanga, Inmazeb | GP1 head / RBS beta7-beta9 | 1 |
| G118E | 118 | G→E | GP1 | mAb114, REGN3471 | Ebanga, Inmazeb | GP1 head / RBS beta7-beta9 | 1 |
| G118R | 118 | G→R | GP1 | mAb114, REGN3471 | Ebanga, Inmazeb | GP1 head / RBS beta7-beta9 | 1 |
| T144M | 144 | T→M | GP1 | mAb114, REGN3471 | Ebanga, Inmazeb | GP1 head / RBS beta7-beta9 | 1 |
| I274M | 274 | I→M | GP1 | REGN3471 | Inmazeb | Glycan cap | 2 |
| W275L | 275 | W→L | GP1 | REGN3471 | Inmazeb | Glycan cap | 2 |
| W291R | 291 | W→R | GP1 | ADI-23774 | MBP134 | Glycan cap loop 279-306 | 1 |
| K510E | 510 | K→E | GP2 | ADI-23774 | MBP134 | GP2 N-terminus / 310 pocket interface | 5 |
| G528E | 528 | G→E | GP2 | ADI-15878 | MBP134 | Fusion loop / GP base | 4 |
| G528R | 528 | G→R | GP2 | REGN3479 | Inmazeb | Fusion loop / GP base | 3 |
| G528S | 528 | G→S | GP2 | ADI-15878 | MBP134 | Fusion loop / GP base | 4 |
| N563Y | 563 | N→Y | GP2 | REGN3479 | Inmazeb | HR1 N-linked glycan | 3 |

1. Rayaprolu et al., Cell Host Microbe 2023 — REGN-EB3 escapes ([10.1016/j.chom.2023.01.002](https://doi.org/10.1016/j.chom.2023.01.002))
2. Davidson et al., J Virol 2015 — I274M, W275L ([10.1128/JVI.01490-15](https://doi.org/10.1128/JVI.01490-15))
3. Rayaprolu et al., Cell Host Microbe 2023 — maftivimab escapes ([10.1016/j.chom.2023.01.002](https://doi.org/10.1016/j.chom.2023.01.002))
4. Bornholdt et al., Cell 2017 — ADI-15878 G528E/S ([10.1016/j.cell.2017.04.037](https://doi.org/10.1016/j.cell.2017.04.037))
5. Bornholdt et al., Nat Struct Mol Biol 2019 — K510E, W291R ([10.1038/s41594-019-0191-4](https://doi.org/10.1038/s41594-019-0191-4))

*Used by ebolaseq `mab-escape-report` vs Makona reference. Antibody links follow `escape_mab_map` in the watchlist.*

## Coordinate system

- **Pos:** UniProt [Q05320](https://www.uniprot.org/uniprotkb/Q05320) (mature GP = residues 33–676). **GP subunit:** GP1 (33–501) or GP2 (502–676).
- **Makona 2014 aa:** [KJ660346.2](https://www.ncbi.nlm.nih.gov/nuccore/KJ660346.2) / [AHX24649.2](https://www.ncbi.nlm.nih.gov/protein/AHX24649.2). At all listed sites, Makona matches Mayinga-76 (Q05320).
- Literature residues are **cross-matched** to this scheme (no offset for antibodies here). If a source uses mature-only numbering, add **32** for Q05320.
- Primary structures: REGN-EB3 escape on Makona; mAb114 / ADI-15946 crystals on Mayinga (= same indices); ADI-15878 contacts from Murin 2018 Table S2 on **EBOV/Mak** cryo-EM; ADI-23774 / ADI-15946 310 pocket from Bornholdt 2019 ([PDB 6MAM](https://www.rcsb.org/structure/6MAM)).
