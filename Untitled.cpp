// AtRTD2 PTC Detection
// PTC detection is based on STOP codons presence upstream of a splicing junction (>50 nt rule). 
// 
// 1) Run add_junction.py to generate a file containing CDS and Junctions information.
// Using the original AtRTD2_19April2016 fasta (.fa) and gtf (.gtf) files we generated an enriched FASTA with CDSs and Junctions. 

AtRTD2 has many isoforms from coding genes lacking any associated CDS. 
On the other hand, some isoforms have shifted START codons (compare to others in the same gene) to avoid short ORFs (and PTCs), 
so it is important to analyze the transcriptome gene-by-gene to assign proper and common START codons to all the isoforms.
    To do so, run reannotate_CDS.py:

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


3) The final step is to assign PTCs to all the isoforms. Run detect_PTCs.py:

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


The results can be easily analyzed in Excel by loading the "_gene_summary.tsv" table. 

