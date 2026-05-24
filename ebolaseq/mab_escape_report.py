"""mAb escape epitope report for ebolaseq GP protein alignments."""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import resources
from typing import Callable, Dict, List, Optional, Set, Tuple
from urllib.request import urlopen

from Bio import AlignIO, SeqIO
from Bio.Align import PairwiseAligner


@dataclass
class PositionEntry:
    pos: int
    domain: str
    region: str
    zaire_aa: str
    mAb114: bool
    REGN3470: bool
    REGN3471: bool
    REGN3479: bool
    contact: bool
    proven_escape: List[str]


@dataclass
class CellResult:
    aa: str
    status: str  # same, changed, gap, escape
    mutation_label: Optional[str] = None
    grantham: Optional[int] = None
    ref_aa: Optional[str] = None


@dataclass
class IsolateSummary:
    seq_id: str
    species: str
    n_changed: int = 0
    n_escape_match: int = 0
    cocktail_concern: bool = False
    mab114_concern: bool = False
    mAb114_hits: int = 0
    regn3470_hits: int = 0
    regn3471_hits: int = 0
    regn3479_hits: int = 0
    max_grantham: int = 0




def _reference_label(metadata: dict) -> str:
    return metadata.get("baseline_label", "Reference")


def _coordinate_mature_sequence(metadata: dict) -> str:
    """UniProt Q05320 mature GP for Q05320 position -> MSA column mapping."""
    if metadata.get("coordinate_mature_seq"):
        return metadata["coordinate_mature_seq"]
    mature_start = int(metadata.get("mature_gp_start", 33))
    acc = metadata.get("coordinate_uniprot", metadata.get("uniprot_accession", "Q05320"))
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urlopen(url, timeout=30) as resp:
        lines = resp.read().decode("utf-8").splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq[mature_start - 1 :]


def _baseline_mature_sequence(metadata: dict) -> str:
    """Makona 2014 mature GP (reference amino acids; file is already residues 33-676)."""
    if metadata.get("baseline_mature_seq"):
        return metadata["baseline_mature_seq"]
    fasta_name = metadata.get("baseline_fasta", "makona_gp_mature.fasta")
    try:
        pkg = resources.files("ebolaseq").joinpath("data", fasta_name)
        lines = pkg.read_text(encoding="utf-8").splitlines()
    except Exception:
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, "data", fasta_name), encoding="utf-8") as fh:
            lines = fh.read().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def _q05320_mature_sequence(metadata: dict) -> str:
    return _coordinate_mature_sequence(metadata)

def load_watchlist(path: Optional[str] = None) -> Tuple[dict, List[PositionEntry]]:
    if path is None:
        try:
            pkg = resources.files("ebolaseq").joinpath("data/gp_epitope_watchlist.json")
            data = json.loads(pkg.read_text(encoding="utf-8"))
        except Exception:
            here = os.path.dirname(os.path.abspath(__file__))
            path = os.path.join(here, "data", "gp_epitope_watchlist.json")
            data = json.loads(open(path, encoding="utf-8").read())
    else:
        data = json.loads(open(path, encoding="utf-8").read())
    entries = [
        PositionEntry(
            pos=e["pos"], domain=e["domain"], region=e["region"], zaire_aa=e["zaire_aa"],
            mAb114=e["mAb114"], REGN3470=e["REGN3470"], REGN3471=e["REGN3471"],
            REGN3479=e["REGN3479"], contact=e["contact"], proven_escape=list(e.get("proven_escape", [])),
        )
        for e in data["positions"]
    ]
    return data.get("metadata", {}), entries


_GRANTHAM_CACHE: Optional[dict] = None


def _load_grantham_matrix() -> dict:
    global _GRANTHAM_CACHE
    if _GRANTHAM_CACHE is not None:
        return _GRANTHAM_CACHE
    try:
        pkg = resources.files("ebolaseq").joinpath("data/grantham_matrix.json")
        _GRANTHAM_CACHE = json.loads(pkg.read_text(encoding="utf-8"))
    except Exception:
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, "data", "grantham_matrix.json"), encoding="utf-8") as fh:
            _GRANTHAM_CACHE = json.load(fh)
    return _GRANTHAM_CACHE


def grantham_distance(ref_aa: str, alt_aa: str) -> Optional[int]:
    """Grantham distance (1974 Table 2); 5–215, higher = more dissimilar."""
    ref_aa, alt_aa = ref_aa.upper(), alt_aa.upper()
    if ref_aa == alt_aa:
        return 0
    m = _load_grantham_matrix()
    if ref_aa not in m or alt_aa not in m:
        return None
    v = m[ref_aa].get(alt_aa) or 0
    if v:
        return int(v)
    v = m[alt_aa].get(ref_aa) or 0
    return int(v) if v else None


def _ungapped(seq: str) -> str:
    return seq.replace("-", "").upper()


def _msa_col_for_ungapped_index(msa_row: str, uidx: int) -> Optional[int]:
    u = 0
    for col, aa in enumerate(msa_row):
        if aa != "-":
            if u == uidx:
                return col
            u += 1
    return None


def _build_q05320_to_msa_col(q05320_mature: str, ref_msa_row: str) -> Dict[int, int]:
    ref_ungapped = _ungapped(ref_msa_row)
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-0.5)
    aln = aligner.align(q05320_mature, ref_ungapped)[0]
    # alignment: q05320 along target 0, ref along query 1
    q_map = {}
    qi = ri = 0
    for a, b in zip(aln[0], aln[1]):
        if a != "-" and b != "-":
            q_map[qi] = ri
            qi += 1
            ri += 1
        elif a != "-":
            qi += 1
        elif b != "-":
            ri += 1
    msa_map = {}
    mature_start = 33
    for qpos in range(mature_start, mature_start + len(q05320_mature)):
        mi = qpos - mature_start
        if mi not in q_map:
            continue
        col = _msa_col_for_ungapped_index(ref_msa_row, q_map[mi])
        if col is not None:
            msa_map[qpos] = col
    return msa_map


def _infer_species(seq_id: str) -> str:
    s = seq_id.lower()
    if "bundibugyo" in s or "/bdbv" in s or "bdbv" in s:
        return "Bundibugyo"
    if "zaire" in s or "zebov" in s or "ebola virus" in s:
        return "Zaire"
    if "sudan" in s:
        return "Sudan"
    if "reston" in s:
        return "Reston"
    if "tai" in s:
        return "Tai Forest"
    parts = seq_id.split("/")
    if len(parts) > 1:
        return parts[1]
    return "unknown"



def _accession(seq_id: str) -> str:
    token = seq_id.split()[0]
    return token.split("/")[0].split("|")[-1]


def _work_dir_from_aln(gp_aln_path: str) -> str:
    return os.path.abspath(os.path.join(os.path.dirname(gp_aln_path), "..", "..", ".."))


def _outgroup_accessions(work_dir: str) -> Set[str]:
    accessions: Set[str] = set()
    fasta_dir = os.path.join(work_dir, "FASTA")
    if not os.path.isdir(fasta_dir):
        return accessions
    for fn in os.listdir(fasta_dir):
        if not fn.startswith("outgroup_") or not fn.endswith((".fasta", ".fa")):
            continue
        for rec in SeqIO.parse(os.path.join(fasta_dir, fn), "fasta"):
            accessions.add(_accession(rec.id))
    return accessions


def _is_outgroup(seq_id: str, outgroup_accs: Set[str]) -> bool:
    if _accession(seq_id) in outgroup_accs:
        return True
    return "outgroup" in seq_id.lower()


def _yes_no(flag: bool) -> str:
    return "Yes" if flag else "No"


def _load_country_map(work_dir: str) -> Dict[str, str]:
    """taxon header or accession -> location from FASTA/location.txt."""
    path = os.path.join(work_dir, "FASTA", "location.txt")
    mapping: Dict[str, str] = {}
    if not os.path.isfile(path):
        return mapping
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("taxon"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            taxon, location = parts[0], parts[1]
            mapping[taxon] = location
            mapping[_accession(taxon)] = location
    return mapping


def _infer_country(seq_id: str, country_map: Dict[str, str]) -> str:
    if seq_id in country_map:
        return country_map[seq_id]
    acc = _accession(seq_id)
    if acc in country_map:
        return country_map[acc]
    header_parts = seq_id.split("/")
    if len(header_parts) >= 3 and header_parts[2]:
        return header_parts[2].replace("_", " ")
    return "unknown"


def find_escape_baseline_reference_id(records, metadata: dict) -> str:
    """Prefer 2014–2016 Makona / West African outbreak sequence in the alignment."""
    preferred = [
        r"kj660346",
        r"makona",
        r"kissidougou",
        r"gueckedou",
        r"/2014/",
        r"/2015/",
        r"/2016/",
        r"west.africa",
    ]
    fallback = [r"zaire", r"zebov", r"ebolavirus", r"ebola virus"]
    legacy = [r"NC_002549", r"mayinga", r"yambuku"]
    for patterns in (preferred, fallback, legacy):
        for rec in records:
            h = rec.id.lower()
            desc = (rec.description or "").lower()
            if any(re.search(p, h) or re.search(p, desc) for p in patterns):
                return rec.id
    return records[0].id


def find_zaire_reference_id(records, metadata: dict) -> str:
    return find_escape_baseline_reference_id(records, metadata)


def _escape_label(ref_aa: str, pos: int, mut_aa: str, catalog: List[str]) -> Optional[str]:
    for lab in catalog:
        if len(lab) >= 4 and lab[0] == ref_aa and lab[-1] == mut_aa:
            try:
                p = int(lab[1:-1])
            except ValueError:
                continue
            if p == pos:
                return lab
    # generic
    if ref_aa != mut_aa and mut_aa not in ("-", "X"):
        return f"{ref_aa}{pos}{mut_aa}"
    return None


def score_alignment(gp_aln_path: str, watchlist: List[PositionEntry], metadata: dict, outgroup_accs: Optional[Set[str]] = None):
    records = list(AlignIO.read(gp_aln_path, "fasta"))
    if not records:
        raise ValueError("Empty GP protein alignment")
    ref_id = find_escape_baseline_reference_id(records, metadata)
    ref_row = None
    for rec in records:
        if rec.id == ref_id:
            ref_row = str(rec.seq)
            break
    if ref_row is None:
        ref_id = records[0].id
        ref_row = str(records[0].seq)

    try:
        q05320_mature = _coordinate_mature_sequence(metadata)
    except Exception as exc:
        print(f"mab-escape-report: could not load Q05320 coordinate sequence ({exc}); using sparse watchlist map")
        mature_start = int(metadata.get("mature_gp_start", 33))
        max_pos = max(e.pos for e in watchlist)
        q05320_mature = ["X"] * (max_pos - mature_start + 1)
        for e in watchlist:
            q05320_mature[e.pos - mature_start] = e.zaire_aa
        q05320_mature = "".join(q05320_mature)

    pos_to_col = _build_q05320_to_msa_col(q05320_mature, ref_row)
    missing = [e.pos for e in watchlist if e.pos not in pos_to_col]
    if missing:
        print(f"mab-escape-report: warning — could not map positions to MSA columns: {missing[:10]}{'...' if len(missing)>10 else ''}")

    matrix: Dict[str, Dict[int, CellResult]] = {}
    summaries: Dict[str, IsolateSummary] = {}

    if outgroup_accs is None:
        outgroup_accs = _outgroup_accessions(_work_dir_from_aln(gp_aln_path))
    excluded_outgroups: List[str] = []

    for rec in records:
        sid = rec.id
        if _is_outgroup(sid, outgroup_accs):
            excluded_outgroups.append(_accession(sid))
            continue
        row = str(rec.seq)
        species = _infer_species(sid + " " + (rec.description or ""))
        matrix[sid] = {}
        summ = IsolateSummary(seq_id=sid, species=species)
        is_ref = sid == ref_id

        for ent in watchlist:
            col = pos_to_col.get(ent.pos)
            if col is None:
                matrix[sid][ent.pos] = CellResult(aa="?", status="gap")
                continue
            aa = row[col] if col < len(row) else "-"
            if aa == "-":
                matrix[sid][ent.pos] = CellResult(aa="-", status="gap")
                if not is_ref and ent.contact:
                    summ.n_changed += 1
                continue
            aa = aa.upper()
            ref_aa = ent.zaire_aa.upper()
            if is_ref or aa == ref_aa:
                matrix[sid][ent.pos] = CellResult(aa=aa, status="same")
                continue
            label = _escape_label(ref_aa, ent.pos, aa, ent.proven_escape)
            status = "escape" if label in ent.proven_escape else "changed"
            g = grantham_distance(ref_aa, aa)
            matrix[sid][ent.pos] = CellResult(aa=aa, status=status, mutation_label=label, grantham=g, ref_aa=ref_aa)
            if ent.contact or ent.proven_escape:
                summ.n_changed += 1
            if label in ent.proven_escape:
                summ.n_escape_match += 1
            if ent.mAb114 and (ent.contact or ent.proven_escape):
                summ.mAb114_hits += 1
            if ent.REGN3470 and (ent.contact or ent.proven_escape):
                summ.regn3470_hits += 1
            if ent.REGN3471 and (ent.contact or ent.proven_escape):
                summ.regn3471_hits += 1
            if ent.REGN3479 and (ent.contact or ent.proven_escape):
                summ.regn3479_hits += 1

        for _ent in watchlist:
            _c = matrix[sid].get(_ent.pos)
            if _c and _c.grantham is not None:
                summ.max_grantham = max(summ.max_grantham, _c.grantham)
        summ.mab114_concern = summ.mAb114_hits > 0
        components_hit = sum([
            summ.regn3470_hits > 0,
            summ.regn3471_hits > 0,
            summ.regn3479_hits > 0,
        ])
        summ.cocktail_concern = components_hit >= 2
        summaries[sid] = summ

    if excluded_outgroups:
        print(
            "mab-escape-report: excluded phylogeny outgroup(s): "
            + ", ".join(sorted(set(excluded_outgroups)))
        )
    return ref_id, pos_to_col, matrix, summaries, excluded_outgroups



def _has_change_in_scope(sid, matrix, watchlist, position_filter):
    for ent in watchlist:
        if not position_filter(ent):
            continue
        if matrix[sid][ent.pos].status in ("changed", "escape"):
            return True
    return False


def _changed_isolates(matrix, summaries, ref_id, watchlist, position_filter=None):
    pf = position_filter or (lambda e: True)
    result = []
    for sid in matrix:
        if sid == ref_id:
            continue
        if summaries[sid].n_changed == 0 and summaries[sid].n_escape_match == 0:
            continue
        if _has_change_in_scope(sid, matrix, watchlist, pf):
            result.append(sid)
    return sorted(result, key=_accession)


def _canonical_isolate_order(matrix: dict, ref_id: str) -> List[str]:
    others = sorted([s for s in matrix if s != ref_id], key=_accession)
    if ref_id in matrix:
        return [ref_id] + others
    return others


def _bool_csv(flag: bool) -> str:
    return "TRUE" if flag else "FALSE"


def _proven_escape_join(ent: PositionEntry) -> str:
    return ";".join(ent.proven_escape) if ent.proven_escape else ""


def _cohort_label(seq_id: str) -> str:
    return "consensus" if _is_user_consensus(seq_id) else "downloaded"


def write_csv_exports(
    output_dir: str,
    gp_aln_path: str,
    metadata: dict,
    watchlist: List[PositionEntry],
    ref_id: str,
    matrix: dict,
    summaries: dict,
) -> List[str]:
    """Write R-friendly CSV tables (tidy long + wide matrix + summaries)."""
    import csv

    work_dir = _work_dir_from_aln(gp_aln_path)
    country_map = _load_country_map(work_dir)
    contact_positions = [ent for ent in watchlist if ent.contact]
    n_contact = len(contact_positions)
    ref_label = _reference_label(metadata)
    ref_accession = metadata.get("reference_accession", metadata.get("baseline_accession", ""))
    isolates = _canonical_isolate_order(matrix, ref_id)
    n_nonref = len([s for s in isolates if s != ref_id])
    written: List[str] = []

    def write_rows(filename: str, header: List[str], rows: List[List]) -> None:
        path = os.path.join(output_dir, filename)
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
        written.append(path)

    run_rows = [[
        ref_label,
        ref_accession,
        _accession(ref_id),
        ref_id,
        gp_aln_path,
        len(watchlist),
        n_nonref,
        n_contact,
        datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    ]]
    write_rows(
        "run_info.csv",
        [
            "baseline_label",
            "baseline_accession",
            "reference_accession",
            "reference_seq_id",
            "gp_alignment",
            "n_watchlist_positions",
            "n_isolates",
            "n_contact_positions",
            "generated_utc",
        ],
        run_rows,
    )

    meta_header = [
        "pos",
        "domain",
        "region",
        "ref_aa",
        "contact",
        "mAb114",
        "REGN3470",
        "REGN3471",
        "REGN3479",
        "proven_escape",
    ]
    wide_rows = []
    for ent in watchlist:
        meta = [
            ent.pos,
            ent.domain,
            ent.region,
            ent.zaire_aa,
            _bool_csv(ent.contact),
            _bool_csv(ent.mAb114),
            _bool_csv(ent.REGN3470),
            _bool_csv(ent.REGN3471),
            _bool_csv(ent.REGN3479),
            _proven_escape_join(ent),
        ]
        cells = []
        for sid in isolates:
            c = matrix[sid][ent.pos]
            if c.status == "escape" and c.mutation_label:
                cells.append(c.mutation_label)
            else:
                cells.append(c.aa)
        wide_rows.append(meta + cells)
    write_rows(
        "epitope_matrix.csv",
        meta_header + [_accession(s) for s in isolates],
        wide_rows,
    )

    cell_header = [
        "accession",
        "seq_id",
        "species",
        "country",
        "cohort",
        "is_reference",
        "pos",
        "domain",
        "region",
        "ref_aa",
        "aa",
        "status",
        "mutation_label",
        "grantham",
        "contact",
        "mAb114",
        "REGN3470",
        "REGN3471",
        "REGN3479",
        "proven_escape",
    ]
    cell_rows = []
    for sid in isolates:
        summ = summaries.get(sid)
        species = summ.species if summ else _infer_species(sid)
        country = _infer_country(sid, country_map)
        for ent in watchlist:
            c = matrix[sid][ent.pos]
            cell_rows.append([
                _accession(sid),
                sid,
                species,
                country,
                _cohort_label(sid),
                _bool_csv(sid == ref_id),
                ent.pos,
                ent.domain,
                ent.region,
                ent.zaire_aa,
                c.aa,
                c.status,
                c.mutation_label or "",
                "" if c.grantham is None else c.grantham,
                _bool_csv(ent.contact),
                _bool_csv(ent.mAb114),
                _bool_csv(ent.REGN3470),
                _bool_csv(ent.REGN3471),
                _bool_csv(ent.REGN3479),
                _proven_escape_join(ent),
            ])
    write_rows("epitope_cells.csv", cell_header, cell_rows)

    sum_header = [
        "accession",
        "seq_id",
        "species",
        "country",
        "cohort",
        "is_reference",
        "n_conserved_contact",
        "n_contact_sites",
        "pct_conserved_contact",
        "epitope_changes",
        "escape_matches",
        "max_grantham",
        "mAb114_hits",
        "REGN3470_hits",
        "REGN3471_hits",
        "REGN3479_hits",
        "mAb114_concern",
        "REGN_EB3_cocktail_concern",
    ]
    sum_rows = []
    for sid in sorted(matrix.keys(), key=_accession):
        summ = summaries[sid]
        n_same = sum(
            1 for ent in contact_positions if matrix[sid][ent.pos].status == "same"
        )
        pct = round(100 * n_same / n_contact, 1) if n_contact else ""
        sum_rows.append([
            _accession(sid),
            sid,
            summ.species,
            _infer_country(sid, country_map),
            _cohort_label(sid),
            _bool_csv(sid == ref_id),
            n_same,
            n_contact,
            pct,
            summ.n_changed,
            summ.n_escape_match,
            summ.max_grantham,
            summ.mAb114_hits,
            summ.regn3470_hits,
            summ.regn3471_hits,
            summ.regn3479_hits,
            _bool_csv(summ.mab114_concern),
            _bool_csv(summ.cocktail_concern),
        ])
    write_rows("isolate_summary.csv", sum_header, sum_rows)

    catalog = _collect_proven_escape_catalog(watchlist)
    nonref = sorted([sid for sid in matrix if sid != ref_id], key=_accession)
    cat_header = [
        "mutation",
        "pos",
        "ref_aa",
        "region",
        "linked_mAbs",
        "in_dataset",
        "n_isolates",
        "accessions",
    ]
    cat_rows = []
    for item in catalog:
        hits = _isolates_with_escape(matrix, item["pos"], item["label"], nonref)
        cat_rows.append([
            item["label"],
            item["pos"],
            item["zaire_aa"],
            item["region"],
            ";".join(item["mabs"]),
            _bool_csv(bool(hits)),
            len(hits),
            ";".join(hits),
        ])
    write_rows("proven_escape_catalog.csv", cat_header, cat_rows)

    return written


def write_csv(path: str, watchlist: List[PositionEntry], ref_id: str, matrix: dict):
    """Backward-compatible wrapper: writes only the wide matrix to *path*."""
    import csv

    isolates = _canonical_isolate_order(matrix, ref_id)
    meta_header = [
        "pos",
        "domain",
        "region",
        "ref_aa",
        "contact",
        "mAb114",
        "REGN3470",
        "REGN3471",
        "REGN3479",
        "proven_escape",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(meta_header + [_accession(s) for s in isolates])
        for ent in watchlist:
            row = [
                ent.pos,
                ent.domain,
                ent.region,
                ent.zaire_aa,
                _bool_csv(ent.contact),
                _bool_csv(ent.mAb114),
                _bool_csv(ent.REGN3470),
                _bool_csv(ent.REGN3471),
                _bool_csv(ent.REGN3479),
                _proven_escape_join(ent),
            ]
            for sid in isolates:
                c = matrix[sid][ent.pos]
                if c.status == "escape" and c.mutation_label:
                    row.append(c.mutation_label)
                else:
                    row.append(c.aa)
            w.writerow(row)




def _cell_class(c: CellResult) -> str:
    return {"same": "same", "changed": "changed", "escape": "escape", "gap": "gap"}.get(c.status, "")


def _format_cell_html(c: CellResult) -> str:
    if c.status in ("changed", "escape"):
        label = c.mutation_label if c.mutation_label else c.aa
        gs = f' <span class="gscore">({c.grantham})</span>' if c.grantham is not None else ""
        return f'<span class="mut">{label}</span>{gs}'
    return c.aa




def _collect_proven_escape_catalog(watchlist: List[PositionEntry]) -> List[dict]:
    """Unique published escape mutations from watchlist with linked antibodies."""
    seen: Set[str] = set()
    items: List[dict] = []
    for ent in watchlist:
        for lab in ent.proven_escape:
            if lab in seen:
                continue
            seen.add(lab)
            mabs: List[str] = []
            for e in watchlist:
                if lab not in e.proven_escape:
                    continue
                if e.mAb114 and "mAb114" not in mabs:
                    mabs.append("mAb114")
                if e.REGN3470 and "mAb3470" not in mabs:
                    mabs.append("mAb3470")
                if e.REGN3471 and "mAb3471" not in mabs:
                    mabs.append("mAb3471")
                if e.REGN3479 and "mAb3479" not in mabs:
                    mabs.append("mAb3479")
            items.append({
                "label": lab,
                "pos": ent.pos,
                "zaire_aa": ent.zaire_aa,
                "region": ent.region,
                "mabs": mabs,
            })
    return sorted(items, key=lambda x: (x["pos"], x["label"]))


def _isolates_with_escape(matrix: dict, pos: int, label: str, isolates: List[str]) -> List[str]:
    found = []
    for sid in isolates:
        c = matrix[sid].get(pos)
        if c and c.status == "escape" and c.mutation_label == label:
            found.append(_accession(sid))
    return sorted(found)


def _proven_escape_section(
    matrix: dict, ref_id: str, watchlist: List[PositionEntry], metadata: Optional[dict] = None
) -> str:
    catalog = _collect_proven_escape_catalog(watchlist)
    if not catalog:
        return ""
    isolates = sorted([sid for sid in matrix.keys() if sid != ref_id], key=_accession)
    rows = []
    n_found_any = 0
    for item in catalog:
        hits = _isolates_with_escape(matrix, item["pos"], item["label"], isolates)
        if hits:
            n_found_any += 1
        in_data = f'<strong>Yes</strong> ({len(hits)} isolate(s))' if hits else "No"
        acc_cell = ", ".join(hits) if hits else "&mdash;"
        mab_str = ", ".join(item["mabs"])
        row_class = "escape" if hits else ""
        rows.append(
            f"<tr><td><strong>{item['label']}</strong></td><td>{item['pos']}</td>"
            f"<td>{item['zaire_aa']}&rarr;{item['label'][-1]}</td><td>{mab_str}</td>"
            f'<td>{item["region"]}</td><td class="{row_class}">{in_data}</td>'
            f"<td>{acc_cell}</td></tr>"
        )
    ref_label = _reference_label(metadata or {})
    return f"""<h2>Catalogued escape mutations (literature)</h2>
<p>Published in vitro escape variants from the watchlist, checked against all analyzed isolates (vs {ref_label}). This is not an exhaustive list of all possible escapes.</p>
<p><strong>Detected in this run:</strong> {n_found_any}/{len(catalog)} catalogued mutation(s).</p>
<table class="proven-escape-table"><thead><tr>
<th>Mutation</th><th>Pos</th><th>Change</th><th>Linked mAb</th><th>Region</th><th>In dataset?</th><th>Accessions</th>
</tr></thead><tbody>{"".join(rows)}</tbody></table>"""



def _is_user_consensus(seq_id: str) -> bool:
    """Sequences added via ebolaseq --c_z, --c_b, etc. (ID contains /consensus)."""
    return "/consensus" in (seq_id or "").lower()


def _cohort_subset(matrix: dict, summaries: dict, cohort: str) -> Tuple[dict, dict]:
    if cohort == "consensus":
        keys = [k for k in matrix if _is_user_consensus(k)]
    else:
        keys = [k for k in matrix if not _is_user_consensus(k)]
    return {k: matrix[k] for k in keys}, {k: summaries[k] for k in keys if k in summaries}


def _render_report_panel(
    matrix: dict,
    summaries: dict,
    ref_id: str,
    watchlist: List[PositionEntry],
    country_map: Optional[Dict[str, str]] = None,
    panel_prefix: str = "",
    metadata: Optional[dict] = None,
) -> str:
    if not matrix:
        return "<p><em>No sequences in this set.</em></p>"
    country_map = country_map or {}
    metadata = metadata or {}
    ref_label = _reference_label(metadata)
    ref_col = f"{ref_label} aa"
    contact_positions = [ent for ent in watchlist if ent.contact]
    n_contact = len(contact_positions)
    all_changed = _changed_isolates(matrix, summaries, ref_id, watchlist)

    sum_rows = []
    for sid in sorted(matrix.keys(), key=_accession):
        s = summaries[sid]
        n_same_contact = sum(
            1 for ent in contact_positions if matrix[sid][ent.pos].status == "same"
        )
        conserved_pct = int(round(100 * n_same_contact / n_contact)) if n_contact else 0
        country = _infer_country(sid, country_map)
        sum_rows.append(
            f'<tr data-accession="{_accession(sid)}" data-species="{s.species}" data-country="{country}"'
            f' data-conserved="{n_same_contact}" data-conserved-pct="{conserved_pct}"'
            f' data-epitope-changes="{s.n_changed}" data-escape="{s.n_escape_match}"'
            f' data-grantham="{s.max_grantham}"'
            f' data-mab114-concern="{1 if s.mab114_concern else 0}"'
            f' data-cocktail-concern="{1 if s.cocktail_concern else 0}">'
            f"<td>{_accession(sid)}</td><td>{s.species}</td>"
            f"<td>{country}</td>"
            f"<td>{n_same_contact}/{n_contact} ({conserved_pct}%)</td>"
            f"<td>{s.n_changed}</td><td>{s.n_escape_match}</td><td>{s.max_grantham}</td>"
            f"<td>{_yes_no(s.mab114_concern)}</td>"
            f"<td>{_yes_no(s.cocktail_concern)}</td></tr>"
        )

    def analysis_isolates():
        return sorted([sid for sid in matrix.keys() if sid != ref_id], key=_accession)

    def antibody_table(title: str, position_filter):
        positions = [ent for ent in watchlist if position_filter(ent)]
        if not positions:
            return f"<h2>{title}</h2><p><em>No positions in watchlist for this antibody.</em></p>"
        slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
        table_id = f"{panel_prefix}-{slug}" if panel_prefix else slug
        isolates = _canonical_isolate_order(matrix, ref_id)
        isolates = [s for s in isolates if s != ref_id]
        changed_isolates = _changed_isolates(matrix, summaries, ref_id, watchlist, position_filter)
        n_pos = len(positions)
        unchanged_rows = changed_rows = 0
        rows = []
        for ent in positions:
            row_max = row_has_change = 0
            for sid in isolates:
                c = matrix[sid][ent.pos]
                if c.status in ("changed", "escape"):
                    row_has_change = 1
                if c.grantham is not None:
                    row_max = max(row_max, c.grantham)
            if row_has_change:
                changed_rows += 1
            else:
                unchanged_rows += 1
            tds = [f"<td>{ent.pos}</td>", f"<td>{ent.domain}</td>", f"<td>{ent.zaire_aa}</td>"]
            for sid in isolates:
                c = matrix[sid][ent.pos]
                g = c.grantham if c.grantham is not None else 0
                tds.append(f'<td class="{_cell_class(c)}" data-grantham="{g}">{_format_cell_html(c)}</td>')
            conserved_attr = "" if row_has_change else ' class="conserved-row conserved-hidden"'
            rows.append(
                f'<tr data-row-max-grantham="{row_max}"{conserved_attr}>'
                + "".join(tds)
                + "</tr>"
            )
        headers = "".join(f"<th>{_accession(sid)}</th>" for sid in isolates)
        conserved_pct = int(round(100 * unchanged_rows / n_pos)) if n_pos else 0
        iso_note = (
            f"{len(changed_isolates)} isolate(s) with &ge;1 change; "
            f"{len(isolates) - len(changed_isolates)} fully conserved at these sites."
        )
        return f"""<h2>{title}</h2>
<p><strong>Conservation vs {ref_label}:</strong> {unchanged_rows}/{n_pos} epitope positions unchanged in all isolates ({conserved_pct}%); {changed_rows} position(s) with &ge;1 change. {iso_note}</p>
<button type="button" class="toggle-conserved" data-table="{table_id}">Show conserved rows</button>
<table class="epitope-table" id="{table_id}"><thead><tr><th>Pos</th><th>Domain</th><th>{ref_col}</th>{headers}</tr></thead>
<tbody>{"".join(rows)}</tbody></table>"""

    no_changes = ""
    if not all_changed:
        no_changes = f"<p><strong>All epitope positions match {ref_label}</strong> in this sequence set.</p>"

    summary_table_id = f"{panel_prefix}-isolate-summary" if panel_prefix else "isolate-summary"
    parts = [
        '<div class="isolate-summary-block">',
        '<h2 class="isolate-summary-heading">Isolate summary '
        '<button type="button" class="toggle-isolate-summary" aria-expanded="false">'
        "Show isolate summary</button></h2>",
        '<div class="isolate-summary-body isolate-summary-collapsed">',
        f"<p>Conserved contact sites = epitope contact positions identical to {ref_label} "
        "(denominator matches epitope-changes scope). Click a column header to sort.</p>",
        f'<table class="summary-table sortable-table" id="{summary_table_id}"><thead><tr>',
        '<th class="sortable" data-sort-key="accession" data-sort-type="text">Accession</th>',
        '<th class="sortable" data-sort-key="species" data-sort-type="text">Species</th>',
        '<th class="sortable" data-sort-key="country" data-sort-type="text">Country</th>',
        '<th class="sortable" data-sort-key="conserved" data-sort-type="number">Conserved contact sites</th>',
        '<th class="sortable" data-sort-key="epitope-changes" data-sort-type="number">Epitope changes</th>',
        '<th class="sortable" data-sort-key="escape" data-sort-type="number">Escape matches</th>',
        '<th class="sortable" data-sort-key="grantham" data-sort-type="number">Max Grantham</th>',
        '<th class="sortable" data-sort-key="mab114-concern" data-sort-type="number">mAb114 (Ebanga) concern</th>',
        '<th class="sortable" data-sort-key="cocktail-concern" data-sort-type="number">REGN-EB3 cocktail concern</th>',
        "</tr></thead><tbody>",
        "".join(sum_rows) if sum_rows else '<tr><td colspan="9"><em>No isolates</em></td></tr>',
        "</tbody></table>",
        "</div></div>",
        no_changes,
        antibody_table("mAb114 (Ebanga / ansuvimab)", lambda e: e.mAb114),
        antibody_table("mAb3470 — atoltivimab (Inmazeb)", lambda e: e.REGN3470),
        antibody_table("mAb3471 — odesivimab (Inmazeb)", lambda e: e.REGN3471),
        antibody_table("mAb3479 — maftivimab (Inmazeb)", lambda e: e.REGN3479),
        _proven_escape_section(matrix, ref_id, watchlist, metadata),
    ]
    return "\n".join(parts)


def write_html(
    path: str,
    gp_aln_path: str,
    metadata: dict,
    watchlist: List[PositionEntry],
    ref_id: str,
    matrix: dict,
    summaries: dict,
    excluded_outgroups: List[str],
):
    work_dir = _work_dir_from_aln(gp_aln_path)
    country_map = _load_country_map(work_dir)

    m_dl, s_dl = _cohort_subset(matrix, summaries, "downloaded")
    m_cs, s_cs = _cohort_subset(matrix, summaries, "consensus")
    n_dl = len([k for k in m_dl if k != ref_id])
    n_cs = len([k for k in m_cs if k != ref_id])
    has_tabs = n_cs > 0

    ref_label = _reference_label(metadata)
    print(
        "mab-escape-report: comparison baseline %s (%s); alignment anchor %s"
        % (ref_label, metadata.get("reference_accession", ""), _accession(ref_id))
    )
    panel_downloaded = _render_report_panel(
        m_dl, s_dl, ref_id, watchlist, country_map, panel_prefix="downloaded", metadata=metadata
    )
    panel_consensus = (
        _render_report_panel(
            m_cs, s_cs, ref_id, watchlist, country_map, panel_prefix="consensus", metadata=metadata
        )
        if has_tabs
        else ""
    )

    tabs_html = ""
    panels_html = ""
    if has_tabs:
        tabs_html = """<div class="tabs" role="tablist">
<button type="button" class="tab-btn active" data-panel="panel-downloaded" role="tab">Downloaded (GenBank)</button>
<button type="button" class="tab-btn" data-panel="panel-consensus" role="tab">Added sequences (consensus)</button>
</div>"""
        panels_html = f"""<div id="panel-downloaded" class="tab-panel active" role="tabpanel">{panel_downloaded}</div>
<div id="panel-consensus" class="tab-panel" role="tabpanel">{panel_consensus}</div>"""
    else:
        panels_html = f'<div class="tab-panel active">{panel_downloaded}</div>'

    html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>GP mAb escape report</title>
<style>
:root {{ --bg:#f8f9fb; --card:#fff; --accent:#1a5fb4; --border:#d0d7de; }}
body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 0; background: var(--bg); color: #222; line-height: 1.5; }}
.wrap {{ max-width: 1280px; margin: 0 auto; padding: 1.25rem 1.5rem 2rem; }}
h1 {{ font-size: 1.5rem; margin: 0 0 0.5rem; }}
h2 {{ font-size: 1.15rem; margin-top: 1.75rem; border-bottom: 2px solid var(--accent); padding-bottom: 0.25rem; }}
.card {{ background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 1rem 1.15rem; margin: 1rem 0; box-shadow: 0 1px 2px rgba(0,0,0,.04); }}
.card h3 {{ margin: 0 0 0.5rem; font-size: 1rem; color: var(--accent); }}
.tabs {{ display: flex; flex-wrap: wrap; gap: 0.35rem; margin: 1.25rem 0 0; }}
.tab-btn {{ padding: 0.45rem 1rem; border: 1px solid var(--border); background: #fff; cursor: pointer; border-radius: 6px; font-size: 0.9rem; }}
.tab-btn.active {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
.tab-panel {{ display: none; margin-top: 0.5rem; }}
.tab-panel.active {{ display: block; }}
.controls {{ display: flex; flex-wrap: wrap; gap: 1rem 2rem; align-items: center; background: #eef2f7; padding: 0.75rem 1rem; border-radius: 6px; margin: 1rem 0; }}
.controls label {{ font-size: 0.9rem; }}
table {{ border-collapse: collapse; font-size: 0.84rem; margin: 0.75rem 0; background: #fff; }}
th, td {{ border: 1px solid var(--border); padding: 0.4rem 0.55rem; text-align: left; vertical-align: top; }}
th {{ background: #e8ecf1; white-space: nowrap; }}
.same {{ background: #d4edda; }} .changed {{ background: #fff3cd; }}
.escape {{ background: #f8d7da; font-weight: 600; }} .gap {{ background: #e9ecef; color: #666; }}
.gscore {{ color: #555; font-size: 0.8em; font-weight: normal; }}
.mut {{ font-weight: 600; }}
.toggle-conserved {{ margin: 0.35rem 0 0.5rem; padding: 0.35rem 0.75rem; font-size: 0.85rem; cursor: pointer; }}
.toggle-isolate-summary {{ margin-left: 0.5rem; padding: 0.3rem 0.65rem; font-size: 0.8rem; font-weight: normal; cursor: pointer; }}
.isolate-summary-collapsed {{ display: none; }}
.isolate-summary-heading {{ display: flex; flex-wrap: wrap; align-items: center; gap: 0.35rem; }}
th.sortable {{ cursor: pointer; user-select: none; }}
th.sortable:hover {{ background: #dce3eb; }}
th.sortable.sorted-asc::after {{ content: " ▲"; font-size: 0.75em; }}
th.sortable.sorted-desc::after {{ content: " ▼"; font-size: 0.75em; }}
.epitope-table tr.conserved-hidden {{ display: none; }}
</style></head><body>
<div class="wrap">
<h1>GP mAb escape epitope report</h1>
<p><strong>In silico only</strong> &mdash; not proof of neutralization.</p>

<div class="card">
<h3>Approved therapies in this report</h3>
<p><strong>Ebanga (ansuvimab)</strong> &mdash; single monoclonal antibody <strong>mAb114</strong>.</p>
<p><strong>Inmazeb (REGN-EB3)</strong> &mdash; <strong>cocktail of three</strong> monoclonal antibodies:
<strong>mAb3470</strong> (atoltivimab), <strong>mAb3471</strong> (odesivimab), and <strong>mAb3479</strong> (maftivimab).</p>
<p><strong>mAb114 vs mAb3471:</strong> epitope tables list only watchlist positions per antibody.
Shared footprint: 112, 116, 118, 143, 144, 146, 270. Odesivimab additionally tracks 111, 272, 274&ndash;275 (I274M, W275L).</p>
</div>

<div class="controls">
<label>Min Grantham (hide rows below): <input type="range" id="minGrantham" min="0" max="215" value="0" step="5"> <span id="minGranthamVal">0</span></label>
<label>Sort epitope rows: <select id="sortRows">
<option value="grantham-desc" selected>Max Grantham (high &rarr; low)</option>
<option value="grantham-asc">Max Grantham (low &rarr; high)</option>
<option value="pos">Position (low &rarr; high)</option>
</select></label>
<button type="button" id="applyBtn">Apply to epitope tables</button>
</div>

{tabs_html}
{panels_html}

</div>
<script>
(function() {{
  document.querySelectorAll('.tab-btn').forEach(btn => {{
    btn.addEventListener('click', () => {{
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
      btn.classList.add('active');
      const panel = document.getElementById(btn.dataset.panel);
      if (panel) panel.classList.add('active');
    }});
  }});
  const slider = document.getElementById('minGrantham');
  const valSpan = document.getElementById('minGranthamVal');
  const sortSel = document.getElementById('sortRows');
  function updateTableRows(tbl, minG) {{
    const tbody = tbl.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    rows.forEach(tr => {{
      if (tr.classList.contains('conserved-row') && tr.classList.contains('conserved-hidden')) {{
        tr.style.display = 'none';
        return;
      }}
      const rowMax = parseInt(tr.dataset.rowMaxGrantham || '0', 10);
      tr.style.display = rowMax >= minG ? '' : 'none';
    }});
    if (sortSel.value === 'grantham-desc') {{
      rows.sort((a,b) => parseInt(b.dataset.rowMaxGrantham||0,10) - parseInt(a.dataset.rowMaxGrantham||0,10));
    }} else if (sortSel.value === 'grantham-asc') {{
      rows.sort((a,b) => parseInt(a.dataset.rowMaxGrantham||0,10) - parseInt(b.dataset.rowMaxGrantham||0,10));
    }} else {{
      rows.sort((a,b) => parseInt(a.cells[0].textContent,10) - parseInt(b.cells[0].textContent,10));
    }}
    rows.forEach(r => tbody.appendChild(r));
  }}
  function apply() {{
    const minG = parseInt(slider.value, 10) || 0;
    valSpan.textContent = minG;
    document.querySelectorAll('.tab-panel.active table.epitope-table').forEach(tbl => updateTableRows(tbl, minG));
  }}
  slider.addEventListener('input', apply);
  sortSel.addEventListener('change', apply);
  document.getElementById('applyBtn').addEventListener('click', apply);
  document.querySelectorAll('.toggle-conserved').forEach(btn => {{
    btn.addEventListener('click', () => {{
      const tbl = document.getElementById(btn.dataset.table);
      if (!tbl) return;
      const show = !!tbl.querySelector('tbody tr.conserved-row.conserved-hidden');
      tbl.querySelectorAll('tbody tr.conserved-row').forEach(tr => {{
        if (show) tr.classList.remove('conserved-hidden');
        else tr.classList.add('conserved-hidden');
      }});
      btn.textContent = show ? 'Hide conserved rows' : 'Show conserved rows';
      updateTableRows(tbl, parseInt(slider.value, 10) || 0);
    }});
  }});
  document.querySelectorAll('.toggle-isolate-summary').forEach(btn => {{
    btn.addEventListener('click', () => {{
      const body = btn.closest('.isolate-summary-block').querySelector('.isolate-summary-body');
      const collapsed = body.classList.toggle('isolate-summary-collapsed');
      btn.setAttribute('aria-expanded', collapsed ? 'false' : 'true');
      btn.textContent = collapsed ? 'Show isolate summary' : 'Hide isolate summary';
    }});
  }});
  function summaryCellVal(tr, sortKey) {{
    return tr.getAttribute('data-' + sortKey) ?? '';
  }}
  function sortSummaryRows(tbl, sortKey, sortType, asc) {{
    const tbody = tbl.querySelector('tbody');
    if (!tbody) return;
    const rows = Array.from(tbody.querySelectorAll('tr')).filter(tr => summaryCellVal(tr, sortKey) !== '');
    rows.sort((a, b) => {{
      let va = summaryCellVal(a, sortKey);
      let vb = summaryCellVal(b, sortKey);
      if (sortType === 'number') {{
        va = parseFloat(va) || 0;
        vb = parseFloat(vb) || 0;
        return asc ? va - vb : vb - va;
      }}
      va = String(va).toLowerCase();
      vb = String(vb).toLowerCase();
      return asc ? va.localeCompare(vb) : vb.localeCompare(va);
    }});
    rows.forEach(r => tbody.appendChild(r));
  }}
  function markSummarySortHeader(tbl, sortKey, asc) {{
    tbl.querySelectorAll('thead th.sortable').forEach(th => {{
      th.classList.remove('sorted-asc', 'sorted-desc');
      if (th.dataset.sortKey === sortKey) {{
        th.classList.add(asc ? 'sorted-asc' : 'sorted-desc');
      }}
    }});
  }}
  document.querySelectorAll('table.summary-table.sortable-table').forEach(tbl => {{
    const defaultKey = 'conserved';
    sortSummaryRows(tbl, defaultKey, 'number', true);
    markSummarySortHeader(tbl, defaultKey, true);
    tbl.querySelectorAll('thead th.sortable').forEach(th => {{
      th.addEventListener('click', () => {{
        const sortKey = th.dataset.sortKey;
        const sortType = th.dataset.sortType || 'number';
        const asc = !th.classList.contains('sorted-asc');
        sortSummaryRows(tbl, sortKey, sortType, asc);
        markSummarySortHeader(tbl, sortKey, asc);
      }});
    }});
  }});
  apply();
}})();
</script>
</body></html>"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)



def run_mab_escape_report(gp_aln_path: str, output_dir: str = "Escape", watchlist_path: Optional[str] = None):
    os.makedirs(output_dir, exist_ok=True)
    metadata, watchlist = load_watchlist(watchlist_path)
    ref_id, pos_to_col, matrix, summaries, excluded = score_alignment(gp_aln_path, watchlist, metadata)
    html_path = os.path.join(output_dir, "gp_mab_escape_report.html")
    write_html(html_path, gp_aln_path, metadata, watchlist, ref_id, matrix, summaries, excluded)
    csv_paths = write_csv_exports(
        output_dir, gp_aln_path, metadata, watchlist, ref_id, matrix, summaries
    )
    print(f"mab-escape-report: wrote {html_path}")
    for p in csv_paths:
        print(f"mab-escape-report: wrote {p}")
    return html_path
