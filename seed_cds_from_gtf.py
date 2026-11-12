#!/usr/bin/env python3
"""
Step 1b: seed CDS start positions into the enriched FASTA.

add_junction.py extracts junctions from the QUASI GTF, which has exon lines only.
The AtRTD2 GTF (AtRTD2_19April2016.gtf, not QUASI) carries the CDS lines. This
script reads the CDS start of each transcript from that GTF, converts it from
genomic to transcript coordinates using the QUASI exon structure, and writes a
CDS=<start>-<start+2> token into the FASTA header. reannotate_CDS.py then
re-scans the real stop codon and propagates the shared start codon across
isoforms of the same gene.

Usage:
    seed_cds_from_gtf.py enriched.fa quasi.gtf atrtd2_with_cds.gtf out_seeded.fa
"""
import sys
import pandas as pd
from Bio import SeqIO

COLS = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]


def read_gtf(path):
    df = pd.read_csv(path, sep="\t", header=None, comment="#", names=COLS,
                     dtype={"start": int, "end": int})
    df["transcript_id"] = df["attribute"].str.extract(r'transcript_id "([^"]+)";')
    return df


def build_exon_structure(quasi):
    exons = quasi[quasi["feature"] == "exon"]
    struct = {}
    for tid, g in exons.groupby("transcript_id"):
        strand = g["strand"].iloc[0]
        g = g.sort_values("start", ascending=(strand != "-"))
        lengths = (g["end"] - g["start"] + 1).values
        cum_before = [0] + list(lengths.cumsum()[:-1])
        struct[tid] = (strand, list(zip(g["start"].values, g["end"].values, cum_before)))
    return struct


def genomic_to_transcript(struct, tid, g_pos):
    if tid not in struct:
        return None
    strand, exon_list = struct[tid]
    for s, e, cum_before in exon_list:
        if s <= g_pos <= e:
            offset = (e - g_pos) if strand == "-" else (g_pos - s)
            return int(cum_before + offset + 1)  # 1-based
    return None


def main(fasta_in, quasi_gtf, cds_gtf, fasta_out):
    struct = build_exon_structure(read_gtf(quasi_gtf))
    cds = read_gtf(cds_gtf)
    cds = cds[cds["feature"] == "CDS"]

    def cds_start(g):
        return int(g["end"].max()) if g["strand"].iloc[0] == "-" else int(g["start"].min())

    starts = cds.groupby("transcript_id").apply(cds_start, include_groups=False)

    seed, missing = {}, 0
    for tid, gpos in starts.items():
        tpos = genomic_to_transcript(struct, tid, gpos)
        if tpos is None:
            missing += 1
        else:
            seed[tid] = tpos
    print(f"CDS-annotated transcripts: {len(starts)} | seeded: {len(seed)} | not mapped: {missing}")

    out, n = [], 0
    for rec in SeqIO.parse(fasta_in, "fasta"):
        rest = rec.description[len(rec.id):].lstrip() if rec.description.startswith(rec.id) else rec.description
        if rec.id in seed:
            c = seed[rec.id]
            rest = rest.replace(" Junctions=", f" CDS={c}-{c + 2} Junctions=")
            n += 1
        rec.description = f"{rec.id} {rest}"
        out.append(rec)
    SeqIO.write(out, fasta_out, "fasta-2line")
    print(f"Wrote {fasta_out} ({n} headers seeded)")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
