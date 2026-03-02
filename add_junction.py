#!/usr/bin/env python3
"""
add_junctions.py
----------------
Enriches a transcriptome FASTA with exon-exon junction coordinates
(in transcript space, 0-based) derived from a GTF file.

Usage:
    python add_junctions.py <fasta_file> <gtf_file> [output_fasta]

Example:
    python add_junctions.py AtRTD2_19April2016.fa AtRTD2_19April2016.gtf

Output:
    AtRTD2_19April2016_enriched.fa  (default, or specify as 3rd argument)

Junction format added to each header:
    Junctions=[341,523,681]   <- 0-based transcript positions of exon-exon boundaries
    Junctions=[]              <- single-exon transcript (no junctions)
"""

import sys
import os
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Suppress pandas SettingWithCopyWarning
pd.options.mode.chained_assignment = None


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

def parse_gtf_to_df(gtf_file: str) -> pd.DataFrame:
    """Parse a GTF file into a DataFrame, extracting transcript_id and gene_id."""
    columns = ['chr', 'source', 'feature', 'start', 'end',
               'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(
        gtf_file, sep='\t', header=None, comment='#', names=columns,
        dtype={'start': int, 'end': int}
    )
    df['transcript_id'] = df['attribute'].str.extract(r'transcript_id "([^"]+)";')
    df['gene_id']       = df['attribute'].str.extract(r'gene_id "([^"]+)";')
    return df


def build_junction_lookup(gtf_df: pd.DataFrame) -> dict:
    """
    Pre-compute junction coordinates (0-based, transcript space) for every
    transcript.  Returns a dict: transcript_id -> list[int]
    """
    exons = gtf_df[gtf_df['feature'] == 'exon'].copy()

    lookup = {}
    for tid, group in exons.groupby('transcript_id'):
        strand = group['strand'].iloc[0]

        # Sort exons in transcript order
        ascending = (strand != '-')
        group = group.sort_values('start', ascending=ascending)

        lengths = (group['end'] - group['start'] + 1).values
        cum_lengths = lengths.cumsum()

        # Junction positions are the cumulative ends of all exons except the last
        junctions = list(cum_lengths[:-1].astype(int))
        lookup[tid] = junctions

    return lookup


# ---------------------------------------------------------------------------
# Per-transcript processing
# ---------------------------------------------------------------------------

def build_enriched_description(record: SeqRecord, junctions: list) -> str:
    """
    Build a clean header string.
    Original description may start with the record ID — strip it to avoid duplication.
    """
    original_desc = record.description
    # Biopython sets description = "<id> <rest>", so strip the leading id
    if original_desc.startswith(record.id):
        original_desc = original_desc[len(record.id):].lstrip()

    junction_str = f"Junctions=[{','.join(map(str, junctions))}]"
    return f"{original_desc} {junction_str}".strip()


def process_record(args):
    """Worker function: takes (record_id, record_seq, record_desc, junctions) tuple."""
    record_id, record_seq, record_desc, junctions = args

    # Reconstruct a minimal SeqRecord for building the description
    class _Stub:
        id = record_id
        description = record_desc

    new_desc = build_enriched_description(_Stub(), junctions)
    return record_id, record_seq, new_desc


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def enrich_fasta(fasta_file: str, gtf_file: str, output_fasta: str,
                 workers: int = 4) -> None:

    print(f"[1/4] Parsing GTF: {gtf_file}")
    gtf_df = parse_gtf_to_df(gtf_file)

    print(f"[2/4] Building junction lookup table...")
    junction_lookup = build_junction_lookup(gtf_df)
    print(f"      {len(junction_lookup):,} transcripts found in GTF.")

    print(f"[3/4] Parsing FASTA: {fasta_file}")
    transcripts = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"      {len(transcripts):,} records found in FASTA.")

    # Prepare args for parallel processing
    tasks = []
    skipped = 0
    for record in transcripts:
        junctions = junction_lookup.get(record.id)
        if junctions is None:
            # Transcript not in GTF — keep record unchanged, add empty junctions
            junctions = []
            skipped += 1
        tasks.append((record.id, str(record.seq), record.description, junctions))

    if skipped:
        print(f"      Warning: {skipped} FASTA records had no GTF entry "
              f"(Junctions=[] assigned).")

    print(f"[4/4] Enriching records (workers={workers})...")
    results = {}

    # Use ProcessPoolExecutor for true parallelism (CPU-bound work)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_id = {executor.submit(process_record, task): task[0]
                        for task in tasks}
        for i, future in enumerate(as_completed(future_to_id), 1):
            tid = future_to_id[future]
            try:
                rec_id, rec_seq, new_desc = future.result()
                results[rec_id] = (rec_seq, new_desc)
            except Exception as exc:
                print(f"      ERROR processing {tid}: {exc}")
            if i % 5000 == 0:
                print(f"      {i:,}/{len(tasks):,} records processed...")

    # Reconstruct SeqRecords in original FASTA order
    enriched_records = []
    for record in transcripts:
        if record.id in results:
            seq, desc = results[record.id]
            enriched_records.append(
                SeqRecord(record.seq, id=record.id, description=desc)
            )
        else:
            enriched_records.append(record)

    print(f"      Writing to: {output_fasta}")
    with open(output_fasta, "w") as out_handle:
        SeqIO.write(enriched_records, out_handle, "fasta")

    print(f"\nDone! Enriched FASTA saved to '{output_fasta}'.")
    print(f"      Total records written: {len(enriched_records):,}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Enrich a transcriptome FASTA with exon-exon junction "
                    "coordinates from a GTF file."
    )
    parser.add_argument("fasta", help="Input transcriptome FASTA file")
    parser.add_argument("gtf",   help="Input GTF annotation file")
    parser.add_argument(
        "output", nargs="?",
        help="Output enriched FASTA file (default: <fasta_basename>_enriched.fa)"
    )
    parser.add_argument(
        "--workers", type=int, default=4,
        help="Number of parallel workers (default: 4)"
    )
    return parser.parse_args()


def default_output_name(fasta_path: str) -> str:
    base, ext = os.path.splitext(fasta_path)
    return f"{base}_enriched{ext}"


if __name__ == "__main__":
    args = parse_args()

    output_path = args.output or default_output_name(args.fasta)

    if not os.path.exists(args.fasta):
        print(f"Error: FASTA file not found: {args.fasta}")
        sys.exit(1)
    if not os.path.exists(args.gtf):
        print(f"Error: GTF file not found: {args.gtf}")
        sys.exit(1)

    enrich_fasta(args.fasta, args.gtf, output_path, workers=args.workers)
