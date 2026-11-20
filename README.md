# AtRTD2 PTC Detection

PTC detection is based on the presence of STOP codons upstream of an exon-exon junction (>50 nt rule).

We implemented a streamlined three-step pipeline in Python to perform the analyses.

## Pipeline

**1) Enriched FASTA with junctions.**

The original AtRTD2_19April2016 fasta (.fa) and gtf (.gtf) files are used (available at https://ics.hutton.ac.uk/atRTD/). The AtRTD2 fasta headers already carry the gene (gene=) and the annotated CDS (CDS=start-end, in transcript coordinates). add_junction.py reads the exons of each transcript from the gtf and adds the exon-exon junctions to the header (Junctions=[...], 0-based positions in transcript space, empty for single-exon transcripts).

```
python add_junction.py AtRTD2_19April2016.fa AtRTD2_19April2016.gtf step1.fa
```

**2) Reannotate CDSs gene by gene.**
AtRTD2 has many isoforms from coding genes lacking any associated CDS. In addition, some isoforms have shifted START codons compared to the other isoforms of the same gene, which avoids short ORFs and PTCs. To prevent isoforms of the same gene from using different ATGs, reannotate_CDS.py works gene by gene. AtRTD2 is used only to identify the START codon of the canonical isoform, and the AtRTD2 CDS end is never used.

  -The canonical isoform of a gene is the one with the longest ORF, scanning in frame from its own annotated START to the first stop codon.
  -The genomic position of the canonical ATG is then looked up in every isoform of the gene using the exon structure in the gtf. Every isoform that contains it, including those with no annotated CDS, uses that ATG as its START.
  -The STOP codon of each isoform is the first in-frame stop codon after its START.
  -An isoform that does not contain the gene’s START (for example, because the start exon is skipped) keeps its own annotated START, and the stop codon is re-scanned from it. Isoforms with no START available, or with no in-frame stop, are left without a CDS (CDS=unknown).
  -Genes with no CDS annotated in any isoform are labeled non-coding and left unchanged.
  -The output FASTA records the change in each header (CDS_original, CDS, CDS_method, GeneType), and a _gene_summary.tsv summarizes the CDS changes per gene.

```
python reannotate_CDS.py step1.fa AtRTD2_19April2016.gtf step2.fa
```

**3) Assign PTCs to all the isoforms.**
detect_ptcs.py calls an isoform as PTC-containing when its stop codon lies more than 50 nt upstream of the last exon-exon junction (or any junction). Single-exon transcripts are never called PTC, and isoforms without a CDS are not called. It writes an annotated FASTA plus transcript-level and gene-level summary tables.
```
python detect_ptcs.py step2.fa
```

The results can be easily visualized and further analyzed by loading the `_gene_summary.tsv` table.
