#!/usr/bin/env python3
"""
correct_cds.py
--------------
Corrects CDS start positions across isoforms of the same gene using a shared
ATG motif strategy, then outputs an enriched FASTA with updated headers and
a gene-level summary TSV.

Strategy:
    For each gene that has AT LEAST ONE isoform with an annotated CDS:
      1. Find the reference isoform (earliest annotated CDS start).
      2. Extract a ~23 nt motif centred on its ATG.
      3. For every isoform in that gene (including those with NO annotated CDS):
           - Search for the motif in the isoform sequence.
           - If found, use that ATG as the CDS start.
           - Scan downstream in-frame for the first stop codon → CDS end.

Header changes per isoform:
    - Had CDS, coords changed   →  CDS_original=X-Y  CDS=new  CDS_method=Corrected_sharedATG  GeneType=Coding
    - Had CDS, coords unchanged →  CDS=X-Y  GeneType=Coding
    - Had NO CDS, motif found   →  CDS=new  CDS_method=Assigned_sharedATG  GeneType=Coding
    - Had NO CDS, motif absent  →  CDS=unknown  CDS_method=NotFound  GeneType=Coding
    - Gene has no CDS at all    →  GeneType=NonCoding  (rest of header untouched)

Outputs:
    <output>.fa               Corrected/enriched FASTA
    <output>_gene_summary.tsv Gene-level summary with CDS stats and GeneType

Usage:
    python correct_cds.py <input.fa> [output.fa]

Example:
    python correct_cds.py AtRTD2_19April2016_enriched.fa
    python correct_cds.py AtRTD2_19April2016_enriched.fa AtRTD2_corrected_CDS.fa
"""

import sys
import os
import argparse
from collections import defaultdict

from typing import Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STOP_CODONS      = {"TAA", "TAG", "TGA"}
MOTIF_UPSTREAM   = 11   # nt upstream of ATG in context motif
MOTIF_DOWNSTREAM = 12   # nt downstream of ATG in context motif


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

def parse_header(record: SeqRecord) -> dict:
    """Extract gene, CDS, and junction info from a FASTA header."""
    desc = record.description
    if desc.startswith(record.id):
        desc = desc[len(record.id):].lstrip()

    info = {
        "gene_id":   None,
        "cds_start": None,   # 1-based inclusive, or None if not annotated
        "cds_end":   None,   # 1-based inclusive, or None if not annotated
        "junctions": [],
    }

    for token in desc.split():
        if "=" not in token:
            continue
        key, val = token.split("=", 1)

        if key == "gene":
            info["gene_id"] = val

        elif key == "CDS":
            if val.lower() not in ("unknown", "none", ""):
                try:
                    s, e = val.split("-")
                    info["cds_start"] = int(s)
                    info["cds_end"]   = int(e)
                except ValueError:
                    pass

        elif key == "Junctions":
            val = val.strip("[]")
            if val:
                try:
                    info["junctions"] = [int(x) for x in val.split(",")]
                except ValueError:
                    pass

    if info["gene_id"] is None:
        info["gene_id"] = record.id.rsplit(".", 1)[0]

    return info


# ---------------------------------------------------------------------------
# Motif and stop-codon helpers
# ---------------------------------------------------------------------------

def extract_motif(seq: str, cds_start_1: int) -> str:
    atg_0 = cds_start_1 - 1
    start = max(0, atg_0 - MOTIF_UPSTREAM)
    end   = min(len(seq), atg_0 + 3 + MOTIF_DOWNSTREAM)
    return seq[start:end].upper()


def find_inframe_stop(seq: str, cds_start_1: int) -> Optional[int]:
    start_0 = cds_start_1 - 1
    for i in range(start_0, len(seq) - 2, 3):
        if seq[i:i+3].upper() in STOP_CODONS:
            return i + 3   # 1-based inclusive end
    return None


def apply_motif(seq: str, motif: str) -> Optional[Tuple[int, int]]:
    idx = seq.upper().find(motif.upper())
    if idx == -1:
        return None
    atg_in_motif = motif.upper().find("ATG")
    if atg_in_motif == -1:
        return None
    start_1 = idx + atg_in_motif + 1
    end_1   = find_inframe_stop(seq, start_1)
    if end_1 is None:
        return None
    return start_1, end_1


# ---------------------------------------------------------------------------
# Header rewriting
# ---------------------------------------------------------------------------

def rewrite_header(record: SeqRecord, info: dict, result: dict,
                   gene_type: str) -> str:
    """
    Rebuild the description according to the correction result and gene type.
    All existing tokens are preserved; CDS= is replaced/added as needed.
    GeneType= is appended.
    """
    desc = record.description
    if desc.startswith(record.id):
        desc = desc[len(record.id):].lstrip()

    method = result["method"]

    # Strip any pre-existing GeneType tag so we don't duplicate it
    tokens = [t for t in desc.split() if not t.startswith("GeneType=")]

    new_tokens = []
    cds_token_found = False

    for token in tokens:
        if token.startswith("CDS="):
            cds_token_found = True
            if method == "Corrected_sharedATG":
                new_tokens.append(f"CDS_original={info['cds_start']}-{info['cds_end']}")
                new_tokens.append(f"CDS={result['start']}-{result['end']}")
                new_tokens.append(f"CDS_method={method}")
            else:
                # unchanged or NoReferenceInGene — keep as-is
                new_tokens.append(token)
        else:
            new_tokens.append(token)

    if not cds_token_found:
        if method == "Assigned_sharedATG":
            new_tokens.append(f"CDS={result['start']}-{result['end']}")
            new_tokens.append(f"CDS_method={method}")
        elif method == "NotFound":
            new_tokens.append("CDS=unknown")
            new_tokens.append("CDS_method=NotFound")
        # NoReferenceInGene: no CDS token added

    # Always append GeneType at the end
    new_tokens.append(f"GeneType={gene_type}")

    return " ".join(new_tokens)


# ---------------------------------------------------------------------------
# Gene summary
# ---------------------------------------------------------------------------

def build_gene_summary(genes: dict, info_map: dict, result_map: dict,
                       gene_type_map: dict) -> pd.DataFrame:
    rows = []
    for gene_id, isoform_ids in genes.items():
        gene_type    = gene_type_map[gene_id]
        n_total      = len(isoform_ids)
        n_annotated  = sum(1 for tid in isoform_ids if info_map[tid]["cds_start"] is not None)
        n_assigned   = sum(1 for tid in isoform_ids
                           if result_map.get(tid, {}).get("method") == "Assigned_sharedATG")
        n_corrected  = sum(1 for tid in isoform_ids
                           if result_map.get(tid, {}).get("method") == "Corrected_sharedATG")
        n_not_found  = sum(1 for tid in isoform_ids
                           if result_map.get(tid, {}).get("method") == "NotFound")
        n_unchanged  = sum(1 for tid in isoform_ids
                           if result_map.get(tid, {}).get("method") == "unchanged")

        rows.append({
            "gene_id":                    gene_id,
            "gene_type":                  gene_type,
            "n_isoforms_total":           n_total,
            "n_isoforms_cds_annotated":   n_annotated,
            "n_isoforms_cds_assigned":    n_assigned,
            "n_isoforms_cds_corrected":   n_corrected,
            "n_isoforms_cds_unchanged":   n_unchanged,
            "n_isoforms_cds_not_found":   n_not_found,
        })

    df = pd.DataFrame(rows)
    df.sort_values("gene_id", inplace=True)
    return df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run(fasta_in: str, fasta_out: str) -> None:

    summary_out = os.path.splitext(fasta_out)[0] + "_gene_summary.tsv"

    print(f"[1/5] Parsing FASTA: {fasta_in}")
    records  = list(SeqIO.parse(fasta_in, "fasta"))
    seq_map  = {rec.id: str(rec.seq).upper() for rec in records}
    info_map = {rec.id: parse_header(rec) for rec in records}
    print(f"      {len(records):,} records loaded.")

    # Group isoforms by gene
    genes: dict = defaultdict(list)
    for rec in records:
        genes[info_map[rec.id]["gene_id"]].append(rec.id)
    print(f"      {len(genes):,} genes found.")

    # -----------------------------------------------------------------------
    # Determine GeneType per gene
    # Coding = at least one isoform has an annotated CDS
    # -----------------------------------------------------------------------
    gene_type_map = {
        gene_id: (
            "Coding" if any(info_map[tid]["cds_start"] is not None for tid in isoform_ids)
            else "NonCoding"
        )
        for gene_id, isoform_ids in genes.items()
    }

    n_coding    = sum(1 for t in gene_type_map.values() if t == "Coding")
    n_noncoding = sum(1 for t in gene_type_map.values() if t == "NonCoding")
    print(f"      Coding genes: {n_coding:,}  |  NonCoding genes: {n_noncoding:,}")

    # -----------------------------------------------------------------------
    # Per-gene correction / CDS assignment
    # -----------------------------------------------------------------------
    print(f"[2/5] Correcting / assigning CDS via shared ATG motif...")

    result_map: dict = {}
    n_corrected = 0
    n_assigned  = 0
    n_not_found = 0

    for gene_id, isoform_ids in genes.items():
        if gene_type_map[gene_id] == "NonCoding":
            for tid in isoform_ids:
                result_map[tid] = {"method": "NoReferenceInGene"}
            continue

        isoforms_with_cds = [
            tid for tid in isoform_ids
            if info_map[tid]["cds_start"] is not None
        ]

        ref_id = min(isoforms_with_cds, key=lambda tid: info_map[tid]["cds_start"])
        motif  = extract_motif(seq_map[ref_id], info_map[ref_id]["cds_start"])

        for tid in isoform_ids:
            info    = info_map[tid]
            had_cds = info["cds_start"] is not None
            coords  = apply_motif(seq_map[tid], motif)

            if coords is None:
                if had_cds:
                    result_map[tid] = {"method": "unchanged"}
                else:
                    result_map[tid] = {"method": "NotFound"}
                    n_not_found += 1
            else:
                corr_start, corr_end = coords
                if had_cds:
                    if corr_start != info["cds_start"] or corr_end != info["cds_end"]:
                        result_map[tid] = {
                            "method": "Corrected_sharedATG",
                            "start":  corr_start,
                            "end":    corr_end,
                        }
                        n_corrected += 1
                    else:
                        result_map[tid] = {"method": "unchanged"}
                else:
                    result_map[tid] = {
                        "method": "Assigned_sharedATG",
                        "start":  corr_start,
                        "end":    corr_end,
                    }
                    n_assigned += 1

    print(f"      CDS corrected (coords changed) : {n_corrected:,} isoforms")
    print(f"      CDS assigned  (was missing)    : {n_assigned:,} isoforms")
    print(f"      CDS not found (motif absent)   : {n_not_found:,} isoforms")

    # -----------------------------------------------------------------------
    # Build output FASTA
    # -----------------------------------------------------------------------
    print(f"[3/5] Building output FASTA...")
    new_records = []
    for rec in records:
        result    = result_map.get(rec.id, {"method": "unchanged"})
        info      = info_map[rec.id]
        gene_type = gene_type_map[info["gene_id"]]
        new_desc  = rewrite_header(rec, info, result, gene_type)
        new_records.append(SeqRecord(rec.seq, id=rec.id, description=new_desc))

    print(f"[4/5] Writing FASTA: {fasta_out}")
    with open(fasta_out, "w") as fh:
        SeqIO.write(new_records, fh, "fasta")

    # -----------------------------------------------------------------------
    # Gene summary TSV
    # -----------------------------------------------------------------------
    print(f"[5/5] Writing gene summary: {summary_out}")
    gene_df = build_gene_summary(genes, info_map, result_map, gene_type_map)
    gene_df.to_csv(summary_out, sep="\t", index=False)

    print(f"\nDone!")
    print(f"  Records written : {len(new_records):,}")
    print(f"  Coding genes    : {n_coding:,}")
    print(f"  NonCoding genes : {n_noncoding:,}")
    print(f"  Outputs:")
    print(f"    {fasta_out}")
    print(f"    {summary_out}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def default_output(fasta_in: str) -> str:
    base, ext = os.path.splitext(fasta_in)
    return f"{base}_corrected_CDS{ext}"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Correct/assign CDS via shared ATG motif, add GeneType flag, output corrected FASTA + gene summary."
    )
    parser.add_argument("fasta", help="Input enriched FASTA (from add_junctions.py)")
    parser.add_argument("output", nargs="?",
                        help="Output FASTA (default: <input>_corrected_CDS.fa)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: file not found: {args.fasta}")
        sys.exit(1)

    run(args.fasta, args.output or default_output(args.fasta))
