# AtRTD2 PTC Detection

PTC detection is based on the presence of STOP codons upstream of a splicing junction (>50 nt rule).

We implemented a streamlined four-step pipeline in Python to perform the analyses.

## Pipeline

**1) Enriched FASTA with junctions.**
Using the original AtRTD2_19April2016 fasta (`.fa`) and gtf (`.gtf`) files (available at <https://ics.hutton.ac.uk/atRTD/>) we generated an enriched FASTA with CDSs and Junctions.
To generate a file containing CDS and Junctions information:

```
python add_junction.py AtRTDv2_QUASI_19April2016.fa AtRTDv2_QUASI_19April2016.gtf step1.fa
```

**2) Seed the CDS start of each transcript.**
The QUASI GTF contains exon lines only, so `add_junction.py` cannot provide a CDS.
The non-QUASI `AtRTD2_19April2016.gtf` carries the CDS lines. `seed_cds_from_gtf.py` reads the CDS start of each
transcript from it, converts it to transcript coordinates using the QUASI exon structure, and writes a `CDS=` token in the header.
Without this step `reannotate_CDS.py` finds no coding genes.

```
python seed_cds_from_gtf.py step1.fa AtRTDv2_QUASI_19April2016.gtf AtRTD2_19April2016.gtf step1_seeded.fa
```

**3) Reannotate CDSs gene by gene.**
AtRTD2 has many isoforms from coding genes lacking any associated CDS.
On the other hand, some isoforms have shifted START codons (compared to others in the same gene) to avoid short ORFs (and PTCs),
so it is important to analyze the transcriptome gene by gene to assign proper and common START codons to all the isoforms.

```
python reannotate_CDS.py step1_seeded.fa step2.fa
```

**4) Assign PTCs to all the isoforms.**

```
python detect_ptcs.py step2.fa
```

The results can be easily visualized and further analyzed by loading the `_gene_summary.tsv` table.
