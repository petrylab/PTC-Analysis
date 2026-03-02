#!/usr/bin/env python3
"""
detect_ptcs.py
--------------
Detects Premature Termination Codons (PTCs) in transcripts from an enriched
FASTA file (produced by add_junctions.py).

PTC Rule (canonical NMD rule):
    The annotated stop codon (= CDS end) is a PTC if it lies more than
    50-55 nt UPSTREAM of the last exon-exon junction (EJC).

    In transcript coordinates (all 0-based):
        PTC = (cds_end_0based) < (last_junction - ptc_distance)

    Single-exon transcripts (no junctions) → never a PTC.

Since CDS coordinates and junction positions are both already in transcript
space inside the FASTA header, NO sequence scanning is required.

Usage:
    python detect_ptcs.py <enriched_fasta> [options]

Example:
    python detect_ptcs.py AtRTD2_19April2016_enriched.fa
    python detect_ptcs.py AtRTD2_19April2016_enriched.fa --ptc-distance 55 --workers 8

Outputs:
    <basename>_ptc_annotated.fa         Enriched FASTA with PTC info in headers
    <basename>_transcript_summary.tsv   Per-transcript PTC summary
    <basename>_gene_summary.tsv         Per-gene PTC summary
"""

import sys
import os
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_PTC_DISTANCE = 50   # nt upstream of last EJC


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

def parse_header(record_id: str, record_description: str) -> dict:
    """Extract key=value pairs from a FASTA header."""
    desc = record_description
    if desc.startswith(record_id):
        desc = desc[len(record_id):].lstrip()

    info = {
        "gene_id":    "unknown",
        "cds_start":  None,   # 1-based inclusive
        "cds_end":    None,   # 1-based inclusive
        "junctions":  [],     # 0-based transcript positions
    }

    for token in desc.split():
        if "=" not in token:
            continue
        key, val = token.split("=", 1)

        if key == "gene":
            info["gene_id"] = val

        elif key == "CDS":
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

    return info


# ---------------------------------------------------------------------------
# Core PTC logic — pure coordinate arithmetic, no sequence scanning
# ---------------------------------------------------------------------------

def detect_ptc(cds_end_1: int, junctions: list, ptc_distance: int) -> dict:
    """
    Parameters
    ----------
    cds_end_1    : CDS end, 1-based inclusive (from header, includes stop codon)
    junctions    : list of 0-based EJC positions in transcript coordinates
    ptc_distance : nt threshold for NMD rule

    Returns
    -------
    dict with has_ptc, cds_end_0based, last_junction, distance_stop_to_last_ejc,
    is_single_exon
    """
    # Convert to 0-based (last nt of stop codon)
    cds_end_0 = cds_end_1 - 1

    if not junctions:
        return {
            "has_ptc":                    False,
            "cds_end_0based":             cds_end_0,
            "last_junction":              None,
            "distance_stop_to_last_ejc":  None,
            "is_single_exon":             True,
        }

    last_junction = junctions[-1]

    # How many nt is the stop codon upstream of the last EJC?
    # Positive = stop is 5' of the junction (upstream in transcript)
    distance = last_junction - cds_end_0

    return {
        "has_ptc":                    distance > ptc_distance,
        "cds_end_0based":             cds_end_0,
        "last_junction":              last_junction,
        "distance_stop_to_last_ejc":  distance,
        "is_single_exon":             False,
    }


# ---------------------------------------------------------------------------
# Worker (called in subprocess)
# ---------------------------------------------------------------------------

def process_transcript(args: tuple) -> dict:
    rec_id, rec_desc, ptc_distance = args

    info = parse_header(rec_id, rec_desc)
    cds_start_1 = info["cds_start"]
    cds_end_1   = info["cds_end"]
    junctions   = info["junctions"]

    base = {
        "transcript_id": rec_id,
        "gene_id":       info["gene_id"],
        "cds_start":     cds_start_1,
        "cds_end":       cds_end_1,
        "n_junctions":   len(junctions),
        "warning":       "",
    }

    if cds_end_1 is None:
        base.update({
            "has_ptc":                   False,
            "cds_end_0based":            None,
            "last_junction":             None,
            "distance_stop_to_last_ejc": None,
            "is_single_exon":            len(junctions) == 0,
            "warning":                   "NO_CDS_ANNOTATION",
        })
        return base

    base.update(detect_ptc(cds_end_1, junctions, ptc_distance))
    return base


# ---------------------------------------------------------------------------
# Header builder
# ---------------------------------------------------------------------------

def build_ptc_header(record_id: str, record_description: str, result: dict) -> str:
    desc = record_description
    if desc.startswith(record_id):
        desc = desc[len(record_id):].lstrip()

    dist = result["distance_stop_to_last_ejc"]
    ptc_tag = (
        f"PTC={'YES' if result['has_ptc'] else 'NO'} "
        f"dist_stop_to_last_EJC={dist if dist is not None else 'NA'}"
    )
    if result["warning"]:
        ptc_tag += f" Warning={result['warning']}"

    return f"{desc} {ptc_tag}".strip()


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def build_transcript_summary(results: list) -> pd.DataFrame:
    rows = [{
        "transcript_id":             r["transcript_id"],
        "gene_id":                   r["gene_id"],
        "cds_start_1based":          r["cds_start"],
        "cds_end_1based":            r["cds_end"],
        "n_junctions":               r["n_junctions"],
        "last_junction_0based":      r["last_junction"],
        "cds_end_0based":            r["cds_end_0based"],
        "distance_stop_to_last_EJC": r["distance_stop_to_last_ejc"],
        "is_single_exon":            r["is_single_exon"],
        "has_ptc":                   r["has_ptc"],
        "warning":                   r["warning"],
    } for r in results]

    df = pd.DataFrame(rows)
    df.sort_values(["gene_id", "transcript_id"], inplace=True)
    return df


def build_gene_summary(transcript_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene_id, grp in transcript_df.groupby("gene_id"):
        n_total       = len(grp)
        n_ptc         = int(grp["has_ptc"].sum())
        n_single_exon = int(grp["is_single_exon"].sum())
        n_no_cds      = int((grp["warning"] == "NO_CDS_ANNOTATION").sum())
        pct_ptc       = round(100 * n_ptc / n_total, 2) if n_total > 0 else 0.0

        if n_ptc == 0:
            status = "no_PTC"
        elif n_ptc == n_total:
            status = "all_PTC"
        else:
            status = "mixed"

        rows.append({
            "gene_id":                gene_id,
            "n_isoforms_total":       n_total,
            "n_isoforms_with_ptc":    n_ptc,
            "n_isoforms_without_ptc": n_total - n_ptc,
            "n_single_exon":          n_single_exon,
            "n_no_cds_annotation":    n_no_cds,
            "pct_isoforms_with_ptc":  pct_ptc,
            "gene_ptc_status":        status,
            "ptc_isoforms":           ";".join(grp[grp["has_ptc"]]["transcript_id"].tolist()),
        })

    df = pd.DataFrame(rows)
    df.sort_values("gene_id", inplace=True)
    return df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pipeline(fasta_file: str, ptc_distance: int, workers: int,
                 out_fasta: str, out_transcript_tsv: str, out_gene_tsv: str) -> None:

    print(f"[1/5] Parsing enriched FASTA: {fasta_file}")
    transcripts = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"      {len(transcripts):,} records found.")

    print(f"[2/5] Detecting PTCs (stop codon > {ptc_distance} nt upstream of last EJC)...")
    tasks = [(r.id, r.description, ptc_distance) for r in transcripts]

    results_map = {}
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_id = {executor.submit(process_transcript, t): t[0] for t in tasks}
        for i, future in enumerate(as_completed(future_to_id), 1):
            tid = future_to_id[future]
            try:
                res = future.result()
                results_map[res["transcript_id"]] = res
            except Exception as exc:
                print(f"      ERROR processing {tid}: {exc}")
            if i % 5000 == 0:
                print(f"      {i:,}/{len(tasks):,} processed...")

    results_ordered = [results_map[r.id] for r in transcripts if r.id in results_map]

    n_ptc    = sum(1 for r in results_ordered if r["has_ptc"])
    n_single = sum(1 for r in results_ordered if r["is_single_exon"])
    print(f"      PTCs found: {n_ptc:,} / {len(results_ordered):,} transcripts "
          f"({n_single:,} single-exon, not eligible for PTC).")

    print(f"[3/5] Writing annotated FASTA: {out_fasta}")
    enriched_records = []
    for record in transcripts:
        result = results_map.get(record.id)
        if result is None:
            enriched_records.append(record)
            continue
        new_desc = build_ptc_header(record.id, record.description, result)
        enriched_records.append(SeqRecord(record.seq, id=record.id, description=new_desc))

    with open(out_fasta, "w") as fh:
        SeqIO.write(enriched_records, fh, "fasta")

    print(f"[4/5] Writing transcript summary: {out_transcript_tsv}")
    transcript_df = build_transcript_summary(results_ordered)
    transcript_df.to_csv(out_transcript_tsv, sep="\t", index=False)

    print(f"[5/5] Writing gene summary: {out_gene_tsv}")
    gene_df = build_gene_summary(transcript_df)
    gene_df.to_csv(out_gene_tsv, sep="\t", index=False)

    # Console summary
    n_genes       = len(gene_df)
    n_genes_any   = int((gene_df["n_isoforms_with_ptc"] > 0).sum())
    n_genes_all   = int((gene_df["gene_ptc_status"] == "all_PTC").sum())
    n_genes_mixed = int((gene_df["gene_ptc_status"] == "mixed").sum())
    n_genes_none  = int((gene_df["gene_ptc_status"] == "no_PTC").sum())

    print(f"\n{'='*55}")
    print(f"  SUMMARY")
    print(f"{'='*55}")
    print(f"  Transcripts analyzed  : {len(results_ordered):,}")
    print(f"  Transcripts with PTC  : {n_ptc:,} ({100*n_ptc/max(len(results_ordered),1):.1f}%)")
    print(f"  Single-exon (no EJCs) : {n_single:,}")
    print(f"  Genes analyzed        : {n_genes:,}")
    print(f"  Genes with any PTC    : {n_genes_any:,}")
    print(f"    - all isoforms PTC  : {n_genes_all:,}")
    print(f"    - mixed             : {n_genes_mixed:,}")
    print(f"    - no PTC            : {n_genes_none:,}")
    print(f"{'='*55}")
    print(f"\nOutputs:")
    print(f"  {out_fasta}")
    print(f"  {out_transcript_tsv}")
    print(f"  {out_gene_tsv}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def default_name(fasta_path: str, suffix: str) -> str:
    return os.path.splitext(fasta_path)[0] + suffix


def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect PTCs using CDS end vs last EJC position from enriched FASTA headers."
    )
    parser.add_argument("fasta",
        help="Enriched FASTA (output of add_junctions.py)")
    parser.add_argument("--ptc-distance", type=int, default=DEFAULT_PTC_DISTANCE,
        help=f"Stop codon must be > N nt upstream of last EJC to call PTC (default: {DEFAULT_PTC_DISTANCE})")
    parser.add_argument("--workers", type=int, default=4,
        help="Parallel workers (default: 4)")
    parser.add_argument("--out-fasta",
        help="Output annotated FASTA (default: <input>_ptc_annotated.fa)")
    parser.add_argument("--out-transcript-tsv",
        help="Output per-transcript TSV (default: <input>_transcript_summary.tsv)")
    parser.add_argument("--out-gene-tsv",
        help="Output per-gene TSV (default: <input>_gene_summary.tsv)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: FASTA file not found: {args.fasta}")
        sys.exit(1)

    run_pipeline(
        fasta_file=args.fasta,
        ptc_distance=args.ptc_distance,
        workers=args.workers,
        out_fasta=args.out_fasta or default_name(args.fasta, "_ptc_annotated.fa"),
        out_transcript_tsv=args.out_transcript_tsv or default_name(args.fasta, "_transcript_summary.tsv"),
        out_gene_tsv=args.out_gene_tsv or default_name(args.fasta, "_gene_summary.tsv"),
    )
