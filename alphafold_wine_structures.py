#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║         AlphaFold Wine Structure Pipeline                                    ║
║         3D structure retrieval for wine lees peptides & yeast proteins       ║
║                                                                              ║
║  STEP 2 of the Wine Peptidome series                                         ║
║  STEP 1 → github.com/314Olamda/WinePeptidome  (UniProt accession retrieval) ║
║                                                                              ║
║  Input  : accession list from WinePeptidome output                           ║
║           (uniprot_kb_entries.tsv) or a plain .txt file, one accession/line  ║
║  Output : AFDB metadata TSV · per-residue pLDDT TSV · PDB structure files   ║
╚══════════════════════════════════════════════════════════════════════════════╝

AlphaFold Database (AFDB) API  —  programmatic access
──────────────────────────────────────────────────────
Source: EBI AlphaFold Database REST API
        https://alphafold.ebi.ac.uk/api/

Endpoints used
  /prediction/{accession}   Metadata + model URLs for a UniProt accession
  PDB download URL          Coordinate file (B-factor column = pLDDT)
  PAE JSON URL              Predicted Aligned Error matrix (optional)

Scientific context
──────────────────
Wine lees released during yeast autolysis contain small proteins and large
peptides (500 Da – 100 kDa) with antifungal and antioxidant bioactivity.
Knowing the predicted 3D structure of these peptides is the first step toward:
  • Understanding surface-exposed antifungal epitopes
  • Guiding in silico docking against grapevine pathogen targets
    (Phaeoacremonium minimum / Phaeomoniella chlamydospora cell wall proteins)
  • Identifying disordered vs. structured regions that correlate with activity

Confidence interpretation (pLDDT)
  ≥ 90   Very high  — likely matches experimental structure
  70–90  Confident  — correct backbone, side chains uncertain
  50–70  Low        — treat with caution
  < 50   Very low   — intrinsically disordered or not modelled

Requirements
────────────
  pip install requests
  Python >= 3.8

Author  : Pol Giménez-Gil
ORCID   : 0000-0002-7720-3733
Affil.  : ISVV – Université de Bordeaux / AUTH / UNIWA
GitHub  : github.com/314Olamda
"""

import csv
import sys
import time
import json
import re
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("ERROR: 'requests' is not installed.  Run: pip install requests")


# ═══════════════════════════════════════════════════════════════════════════════
#  CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

AFDB_API       = "https://alphafold.ebi.ac.uk/api"

# Input: point to the WinePeptidome output TSV, or a plain .txt (one acc/line)
INPUT_FILE     = Path("../WinePeptidome/output/uniprot_kb_entries.tsv")

# Accession column name if using the WinePeptidome TSV
ACC_COLUMN     = "Entry"

# Confidence threshold — only download PDB for entries above this pLDDT mean
PLDDT_MIN      = 50.0

# Set to True to download the PAE (predicted aligned error) JSON as well
DOWNLOAD_PAE   = False

# Set to True to download PDB coordinate files (can be large for many entries)
DOWNLOAD_PDB   = True

OUTPUT_DIR     = Path(__file__).parent / "output"
PDB_DIR        = OUTPUT_DIR / "pdb_structures"

SLEEP_S        = 0.35    # polite delay between requests


# ═══════════════════════════════════════════════════════════════════════════════
#  UTILITY
# ═══════════════════════════════════════════════════════════════════════════════

def _get(url: str, params: dict = None, binary: bool = False,
         retries: int = 4) -> requests.Response | None:
    headers = {} if binary else {"Accept": "application/json"}
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=60)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None     # no AFDB entry for this accession — expected
            if r.status_code == 429:
                wait = 2 ** attempt
                print(f"  [429] rate-limited — waiting {wait}s")
                time.sleep(wait)
                continue
            print(f"  [HTTP {r.status_code}] {url[:80]}")
            time.sleep(2 ** attempt)
        except requests.RequestException as exc:
            print(f"  [connection error] {exc}")
            time.sleep(2 ** attempt)
    return None


def save_tsv(rows: list, fields: list, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields,
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  ✓ {len(rows):,} rows → {path.name}")


# ═══════════════════════════════════════════════════════════════════════════════
#  1.  LOAD ACCESSIONS
# ═══════════════════════════════════════════════════════════════════════════════

def load_accessions(path: Path) -> list:
    """
    Accept either:
      • WinePeptidome TSV (has a header with column 'Entry')
      • Plain .txt file, one UniProt accession per line
    Returns a deduplicated list of accession strings.
    """
    accs = []
    suffix = path.suffix.lower()

    if suffix == ".tsv":
        with open(path, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                acc = row.get(ACC_COLUMN, "").strip()
                if acc:
                    accs.append(acc)
        print(f"  Loaded {len(accs):,} accessions from TSV column '{ACC_COLUMN}'")
    else:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                acc = line.strip()
                if acc and not acc.startswith("#"):
                    accs.append(acc)
        print(f"  Loaded {len(accs):,} accessions from plain text file")

    # Deduplicate while preserving order
    seen, unique = set(), []
    for a in accs:
        if a not in seen:
            seen.add(a)
            unique.append(a)
    print(f"  → {len(unique):,} unique accessions after deduplication")
    return unique


# ═══════════════════════════════════════════════════════════════════════════════
#  2.  ALPHAFOLD DB METADATA
# ═══════════════════════════════════════════════════════════════════════════════

def fetch_afdb_metadata(accession: str) -> dict | None:
    """
    Query AFDB REST API for one UniProt accession.
    Returns the first prediction entry (most recent version) or None.

    Response fields used:
      entryId, gene, uniprotAccession, uniprotDescription,
      taxId, organismScientificName,
      pdbUrl, cifUrl, paeDocUrl,
      modelCreatedDate, latestVersion,
      isReviewed, isReferenceProteome,
      uniprotStart, uniprotEnd
    """
    r = _get(f"{AFDB_API}/prediction/{accession}")
    if r is None:
        return None
    data = r.json()
    if not data:
        return None
    # API returns a list; take the latest version entry
    entry = data[0] if isinstance(data, list) else data
    return entry


# ═══════════════════════════════════════════════════════════════════════════════
#  3.  pLDDT EXTRACTION FROM PDB B-FACTOR COLUMN
# ═══════════════════════════════════════════════════════════════════════════════

def extract_plddt_from_pdb(pdb_text: str) -> dict:
    """
    In AlphaFold PDB files the B-factor column encodes per-residue pLDDT.
    Parses ATOM records and returns:
      {
        'mean_plddt':   float,
        'min_plddt':    float,
        'max_plddt':    float,
        'n_residues':   int,
        'high_conf_pct': float,   # % residues >= 70
        'per_residue':  list of (res_num, plddt)
      }
    """
    residue_scores = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            res_num  = int(line[22:26].strip())
            b_factor = float(line[60:66].strip())
            # Keep one score per residue (first ATOM per residue = CA typically)
            if res_num not in residue_scores:
                residue_scores[res_num] = b_factor
        except (ValueError, IndexError):
            continue

    if not residue_scores:
        return {}

    scores = list(residue_scores.values())
    per_residue = sorted(residue_scores.items())
    high_conf   = sum(1 for s in scores if s >= 70.0)

    return {
        "mean_plddt":    round(sum(scores) / len(scores), 2),
        "min_plddt":     round(min(scores), 2),
        "max_plddt":     round(max(scores), 2),
        "n_residues":    len(scores),
        "high_conf_pct": round(100 * high_conf / len(scores), 1),
        "per_residue":   per_residue,
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  4.  PDB DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════════════

def download_pdb(accession: str, pdb_url: str) -> str | None:
    """Download PDB file and save to PDB_DIR. Returns file text or None."""
    r = _get(pdb_url, binary=True)
    if r is None:
        return None
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    fname = PDB_DIR / f"{accession}.pdb"
    fname.write_bytes(r.content)
    return r.text


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("\n" + "═" * 65)
    print("  AlphaFold Wine Structure Pipeline")
    print(f"  Input  : {INPUT_FILE}")
    print(f"  pLDDT threshold for download : ≥ {PLDDT_MIN}")
    print("═" * 65)

    if not INPUT_FILE.exists():
        sys.exit(
            f"\nERROR: Input file not found:\n  {INPUT_FILE.resolve()}\n\n"
            "Run the WinePeptidome pipeline first, or point INPUT_FILE to\n"
            "a plain .txt file with one UniProt accession per line."
        )

    accessions = load_accessions(INPUT_FILE)

    metadata_rows    = []
    plddt_rows       = []
    per_residue_rows = []
    not_in_afdb      = []

    print(f"\n  Querying AFDB for {len(accessions):,} accessions …\n")

    for i, acc in enumerate(accessions, 1):
        print(f"  [{i:>5}/{len(accessions)}]  {acc}", end="  ", flush=True)

        entry = fetch_afdb_metadata(acc)
        time.sleep(SLEEP_S)

        if entry is None:
            print("✗ not in AFDB")
            not_in_afdb.append(acc)
            continue

        pdb_url  = entry.get("pdbUrl",    "")
        pae_url  = entry.get("paeDocUrl", "")
        entry_id = entry.get("entryId",   acc)

        # ── Download PDB & extract pLDDT ────────────────────────────────────
        plddt_stats = {}
        if DOWNLOAD_PDB and pdb_url:
            pdb_text = download_pdb(acc, pdb_url)
            time.sleep(SLEEP_S)
            if pdb_text:
                plddt_stats = extract_plddt_from_pdb(pdb_text)
                # Store per-residue data for detailed TSV
                for res_num, score in plddt_stats.get("per_residue", []):
                    per_residue_rows.append({
                        "accession":  acc,
                        "entry_id":   entry_id,
                        "residue":    res_num,
                        "plddt":      score,
                        "confident":  score >= 70.0,
                    })

        mean_plddt = plddt_stats.get("mean_plddt", "")

        # ── Filter by confidence ─────────────────────────────────────────────
        if mean_plddt != "" and mean_plddt < PLDDT_MIN:
            print(f"⚠  pLDDT={mean_plddt} below threshold — skipped download")
        else:
            print(f"✓  pLDDT={mean_plddt}" if mean_plddt != "" else "✓")

        # ── Metadata row ─────────────────────────────────────────────────────
        metadata_rows.append({
            "accession":              acc,
            "entry_id":               entry_id,
            "gene":                   entry.get("gene", ""),
            "protein_description":    entry.get("uniprotDescription", ""),
            "organism":               entry.get("organismScientificName", ""),
            "tax_id":                 entry.get("taxId", ""),
            "uniprot_start":          entry.get("uniprotStart", ""),
            "uniprot_end":            entry.get("uniprotEnd", ""),
            "model_version":          entry.get("latestVersion", ""),
            "model_created_date":     entry.get("modelCreatedDate", ""),
            "is_reviewed":            entry.get("isReviewed", ""),
            "is_reference_proteome":  entry.get("isReferenceProteome", ""),
            "mean_plddt":             mean_plddt,
            "min_plddt":              plddt_stats.get("min_plddt", ""),
            "max_plddt":              plddt_stats.get("max_plddt", ""),
            "n_residues":             plddt_stats.get("n_residues", ""),
            "high_conf_pct":          plddt_stats.get("high_conf_pct", ""),
            "pdb_url":                pdb_url,
            "pae_url":                pae_url,
        })

        # ── pLDDT summary row ────────────────────────────────────────────────
        if plddt_stats:
            plddt_rows.append({
                "accession":      acc,
                "entry_id":       entry_id,
                "n_residues":     plddt_stats["n_residues"],
                "mean_plddt":     plddt_stats["mean_plddt"],
                "min_plddt":      plddt_stats["min_plddt"],
                "max_plddt":      plddt_stats["max_plddt"],
                "high_conf_pct":  plddt_stats["high_conf_pct"],
                "confidence_tier": (
                    "very_high" if plddt_stats["mean_plddt"] >= 90 else
                    "confident" if plddt_stats["mean_plddt"] >= 70 else
                    "low"       if plddt_stats["mean_plddt"] >= 50 else
                    "very_low"
                ),
            })

        # Optional PAE download
        if DOWNLOAD_PAE and pae_url:
            r_pae = _get(pae_url)
            time.sleep(SLEEP_S)
            if r_pae:
                pae_path = OUTPUT_DIR / "pae" / f"{acc}_pae.json"
                pae_path.parent.mkdir(parents=True, exist_ok=True)
                pae_path.write_text(r_pae.text, encoding="utf-8")

    # ── Save outputs ──────────────────────────────────────────────────────────
    print(f"\n{'═'*65}")
    print("  Saving outputs …")

    meta_fields = [
        "accession", "entry_id", "gene", "protein_description",
        "organism", "tax_id", "uniprot_start", "uniprot_end",
        "model_version", "model_created_date",
        "is_reviewed", "is_reference_proteome",
        "mean_plddt", "min_plddt", "max_plddt",
        "n_residues", "high_conf_pct",
        "pdb_url", "pae_url",
    ]
    save_tsv(metadata_rows,  meta_fields,
             OUTPUT_DIR / "afdb_metadata.tsv")

    plddt_fields = [
        "accession", "entry_id", "n_residues",
        "mean_plddt", "min_plddt", "max_plddt",
        "high_conf_pct", "confidence_tier",
    ]
    save_tsv(plddt_rows, plddt_fields,
             OUTPUT_DIR / "plddt_summary.tsv")

    per_res_fields = ["accession", "entry_id", "residue", "plddt", "confident"]
    save_tsv(per_residue_rows, per_res_fields,
             OUTPUT_DIR / "plddt_per_residue.tsv")

    # Summary
    high  = sum(1 for r in plddt_rows if r["confidence_tier"] == "very_high")
    conf  = sum(1 for r in plddt_rows if r["confidence_tier"] == "confident")
    low   = sum(1 for r in plddt_rows if r["confidence_tier"] == "low")
    vlow  = sum(1 for r in plddt_rows if r["confidence_tier"] == "very_low")
    pdbs  = len(list(PDB_DIR.glob("*.pdb"))) if PDB_DIR.exists() else 0

    summary = [
        "AlphaFold Wine Structure Pipeline — Run Summary",
        f"Input accessions    : {len(accessions):,}",
        f"Found in AFDB       : {len(metadata_rows):,}",
        f"Not in AFDB         : {len(not_in_afdb):,}",
        f"PDB files saved     : {pdbs:,}",
        "",
        "pLDDT confidence distribution (structures with PDB download)",
        f"  Very high (≥90)   : {high:,}",
        f"  Confident (70–90) : {conf:,}",
        f"  Low (50–70)       : {low:,}",
        f"  Very low (<50)    : {vlow:,}",
        "",
        "Not found in AFDB:",
    ] + (not_in_afdb if not_in_afdb else ["  (none)"])

    summary_path = OUTPUT_DIR / "pipeline_summary.txt"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("\n".join(summary), encoding="utf-8")
    print(f"  ✓ Summary → pipeline_summary.txt")

    print("\n" + "═" * 65)
    print(f"  DONE")
    print(f"  Accessions queried    : {len(accessions):,}")
    print(f"  Found in AFDB         : {len(metadata_rows):,}")
    print(f"  Not in AFDB           : {len(not_in_afdb):,}")
    print(f"  PDB files saved       : {pdbs:,}")
    print(f"  Per-residue pLDDT rows: {len(per_residue_rows):,}")
    print("═" * 65)
    print(f"\n  All files → {OUTPUT_DIR.resolve()}\n")


if __name__ == "__main__":
    main()
