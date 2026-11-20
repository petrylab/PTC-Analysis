#!/usr/bin/env python3
"""
reannotate_CDS.py
-----------------
Gives every isoform of a gene the same START codon and re-scans the STOP codon
from it. The AtRTD2 CDS end is never used; AtRTD2 is only used to identify the
START (ATG) of the canonical isoform.

Strategy, for each gene with AT LEAST ONE isoform with an annotated CDS:
  1. Canonical (reference) isoform: the isoform with the longest ORF, scanning in
     frame from its own AtRTD2 START to the first stop codon.
  2. The genomic position of its ATG is taken as the gene's START.
  3. For every isoform of the gene (including those with NO annotated CDS):
       - If the gene's START is contained in the isoform (and is an ATG in its
         sequence), that ATG is used and the first in-frame stop codon gives
         the CDS end.
       - Otherwise (the isoform cannot use the gene's START, e.g. the start exon
         is skipped), the isoform keeps its own AtRTD2 START and the stop codon
         is re-scanned from it. Isoforms with no AtRTD2 START either get no CDS.
Isoforms are matched by genomic position (exon structure from the GTF), not by
sequence context, so a START is still found when splicing changes the sequence
around it.

Header changes per isoform:
    - Had CDS, coords changed   ->  CDS_original=X-Y  CDS=new  CDS_method=Corrected_sharedATG
                                    (or Corrected_ownSTART if it keeps its own START)  GeneType=Coding
    - Had CDS, coords unchanged ->  CDS=X-Y  GeneType=Coding
    - Had NO CDS, START found   ->  CDS=new  CDS_method=Assigned_sharedATG  GeneType=Coding
    - No START / no stop found  ->  CDS=unknown  CDS_method=NotFound  GeneType=Coding
    - Gene has no CDS at all    ->  GeneType=NonCoding  (rest of header untouched)

Outputs:
    <output>.fa               Reannotated FASTA
    <output>_gene_summary.tsv Gene-level summary with CDS stats and GeneType

Usage:
    python reannotate_CDS.py <enriched.fa> <annotation.gtf> [output.fa]

Example:
    python reannotate_CDS.py AtRTD2_19April2016_enriched.fa AtRTD2_19April2016.gtf
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

def find_inframe_stop(seq: str, cds_start_1: int) -> Optional[int]:
    start_0 = cds_start_1 - 1
    for i in range(start_0, len(seq) - 2, 3):
        if seq[i:i+3].upper() in STOP_CODONS:
            return i + 3   # 1-based inclusive end
    return None


def load_exons(gtf_path: str) -> dict:
    """transcript_id -> (strand, [(exon_start, exon_end, bases_before_exon), ...]) in transcript order."""
    cols = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gtf = pd.read_csv(gtf_path, sep="\t", header=None, comment="#", names=cols,
                      dtype={"start": int, "end": int})
    gtf = gtf[gtf["feature"] == "exon"].copy()
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"')
    exons = {}
    for tid, g in gtf.groupby("transcript_id"):
        strand = g["strand"].iloc[0]
        g = g.sort_values("start", ascending=(strand != "-"))
        lengths = (g["end"] - g["start"] + 1).values
        before = [0] + list(lengths.cumsum()[:-1])
        exons[tid] = (strand, list(zip(g["start"].values, g["end"].values, before)))
    return exons


def tx_to_genomic(exons: dict, tid: str, pos_1: int) -> Optional[int]:
    """1-based transcript position -> genomic coordinate."""
    if tid not in exons:
        return None
    strand, exon_list = exons[tid]
    for s, e, before in exon_list:
        if before < pos_1 <= before + (e - s + 1):
            offset = pos_1 - before - 1
            return int(e - offset) if strand == "-" else int(s + offset)
    return None


def genomic_to_tx(exons: dict, tid: str, g_pos: int) -> Optional[int]:
    """Genomic coordinate -> 1-based transcript position (None if not in an exon)."""
    if tid not in exons:
        return None
    strand, exon_list = exons[tid]
    for s, e, before in exon_list:
        if s <= g_pos <= e:
            offset = (e - g_pos) if strand == "-" else (g_pos - s)
            return int(before + offset + 1)
    return None


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
            if method in ("Corrected_sharedATG", "Corrected_ownSTART"):
                new_tokens.append(f"CDS_original={info['cds_start']}-{info['cds_end']}")
                new_tokens.append(f"CDS={result['start']}-{result['end']}")
                new_tokens.append(f"CDS_method={method}")
            elif method == "NotFound":
                new_tokens.append(f"CDS_original={info['cds_start']}-{info['cds_end']}")
                new_tokens.append("CDS=unknown")
                new_tokens.append("CDS_method=NotFound")
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
        n_own_start  = sum(1 for tid in isoform_ids
                           if result_map.get(tid, {}).get("method") == "Corrected_ownSTART")

        rows.append({
            "gene_id":                    gene_id,
            "gene_type":                  gene_type,
            "n_isoforms_total":           n_total,
            "n_isoforms_cds_annotated":   n_annotated,
            "n_isoforms_cds_assigned":    n_assigned,
            "n_isoforms_cds_corrected":   n_corrected,
            "n_isoforms_cds_own_start":   n_own_start,
            "n_isoforms_cds_unchanged":   n_unchanged,
            "n_isoforms_cds_not_found":   n_not_found,
        })

    df = pd.DataFrame(rows)
    df.sort_values("gene_id", inplace=True)
    return df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run(fasta_in: str, gtf_in: str, fasta_out: str) -> None:

    summary_out = os.path.splitext(fasta_out)[0] + "_gene_summary.tsv"

    print(f"[1/5] Parsing FASTA: {fasta_in}")
    records  = list(SeqIO.parse(fasta_in, "fasta"))
    seq_map  = {rec.id: str(rec.seq).upper() for rec in records}
    info_map = {rec.id: parse_header(rec) for rec in records}
    print(f"      {len(records):,} records loaded.")
    print(f"      Loading exon structures: {gtf_in}")
    exons = load_exons(gtf_in)

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
    print(f"[2/5] Assigning the gene's START to every isoform and re-scanning STOP codons...")

    result_map: dict = {}
    n_corrected = 0
    n_assigned  = 0
    n_not_found = 0
    n_own_start = 0

    for gene_id, isoform_ids in genes.items():
        if gene_type_map[gene_id] == "NonCoding":
            for tid in isoform_ids:
                result_map[tid] = {"method": "NoReferenceInGene"}
            continue

        # Canonical isoform: longest ORF from its own AtRTD2 START (AtRTD2 CDS end is ignored)
        best_len, ref_id = -1, None
        for tid in isoform_ids:
            cs = info_map[tid]["cds_start"]
            if cs is None:
                continue
            end = find_inframe_stop(seq_map[tid], cs)
            orf_len = (end - cs + 1) if end else 0
            if orf_len > best_len:
                best_len, ref_id = orf_len, tid
        ref_atg_genomic = tx_to_genomic(exons, ref_id, info_map[ref_id]["cds_start"])

        for tid in isoform_ids:
            info    = info_map[tid]
            had_cds = info["cds_start"] is not None
            seq     = seq_map[tid]

            pos = genomic_to_tx(exons, tid, ref_atg_genomic) if ref_atg_genomic else None
            if pos is not None and seq[pos - 1:pos + 2] == "ATG":
                start, method_new = pos, "shared"
            elif had_cds:
                start, method_new = info["cds_start"], "own"
            else:
                result_map[tid] = {"method": "NotFound"}
                n_not_found += 1
                continue

            end = find_inframe_stop(seq, start)
            if end is None:
                result_map[tid] = {"method": "NotFound"}
                n_not_found += 1
                continue

            if had_cds:
                if start == info["cds_start"] and end == info["cds_end"]:
                    result_map[tid] = {"method": "unchanged"}
                else:
                    result_map[tid] = {
                        "method": "Corrected_sharedATG" if method_new == "shared" else "Corrected_ownSTART",
                        "start": start, "end": end,
                    }
                    if method_new == "shared":
                        n_corrected += 1
                    else:
                        n_own_start += 1
            else:
                result_map[tid] = {"method": "Assigned_sharedATG", "start": start, "end": end}
                n_assigned += 1

    print(f"      CDS corrected, gene's START (coords changed) : {n_corrected:,} isoforms")
    print(f"      CDS corrected, own START (gene's START absent): {n_own_start:,} isoforms")
    print(f"      CDS assigned  (was missing)                   : {n_assigned:,} isoforms")
    print(f"      No CDS (no START or no in-frame stop)         : {n_not_found:,} isoforms")

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
        description="Give all isoforms of a gene the same START (matched by genomic position), re-scan STOP codons, add GeneType flag, output reannotated FASTA + gene summary."
    )
    parser.add_argument("fasta", help="Input enriched FASTA (from add_junction.py)")
    parser.add_argument("gtf", help="GTF annotation used for the FASTA (exon lines are used)")
    parser.add_argument("output", nargs="?",
                        help="Output FASTA (default: <input>_corrected_CDS.fa)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not os.path.exists(args.fasta):
        print(f"Error: file not found: {args.fasta}")
        sys.exit(1)

    if not os.path.exists(args.gtf):
        print(f"Error: file not found: {args.gtf}")
        sys.exit(1)

    run(args.fasta, args.gtf, args.output or default_output(args.fasta))
