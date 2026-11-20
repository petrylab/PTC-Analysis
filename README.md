# AtRTD2 PTC Detection

PTC detection is based on the presence of STOP codons upstream of the last exon-exon junction (>50 nt rule).

We implemented a streamlined three-step pipeline in Python to perform the analyses.

## Pipeline

**1) Enriched FASTA with junctions.**
Using the original AtRTD2_19April2016 fasta (`.fa`) and gtf (`.gtf`) files (available at <https://ics.hutton.ac.uk/atRTD/>) we generated an enriched FASTA with CDSs and Junctions.
The AtRTD2 fasta headers already carry the annotated CDS (`CDS=start-end`, transcript coordinates); `add_junction.py` adds the exon-exon junctions from the gtf.

```
python add_junction.py AtRTD2_19April2016.fa AtRTD2_19April2016.gtf step1.fa
```

**2) Reannotate CDSs gene by gene.**
AtRTD2 has many isoforms from coding genes lacking any associated CDS.
On the other hand, some isoforms have shifted START codons (compared to others in the same gene) to avoid short ORFs (and PTCs),
so it is important to analyze the transcriptome gene by gene to assign proper and common START codons to all the isoforms.

The AtRTD2 CDS is only used to identify the START codon of the canonical isoform of each gene (the isoform with the longest ORF, scanning in frame from its annotated START).
That START is then used for all the isoforms of the gene that contain it, matched by genomic position using the exon structure in the gtf, and the STOP codon is found by scanning in frame from it.
The AtRTD2 CDS end is never used. An isoform that does not contain the gene's START (for example, because the start exon is skipped) keeps its own annotated START and the STOP is re-scanned from it.

```
python reannotate_CDS.py step1.fa AtRTD2_19April2016.gtf step2.fa
```

**3) Assign PTCs to all the isoforms.**

```
python detect_ptcs.py step2.fa
```

The results can be easily visualized and further analyzed by loading the `_gene_summary.tsv` table.
