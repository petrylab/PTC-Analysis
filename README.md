AtRTD2 PTC Detection

PTC detection is based on STOP codons presence upstream of a splicing junction (>50 nt rule). 

We implemented a streamlined three-step pipeline in Python to perform the analyses.

1) Using the original AtRTD2_19April2016 fasta (.fa) and gtf (.gtf) files (available at https://ics.hutton.ac.uk/atRTD/) we generated an enriched FASTA with CDSs and Junctions.
To generate a file containing CDS and Junctions information...

    Run add_junction.py

2) The QUASI GTF contains exon lines only, so add_junction.py cannot provide a CDS.
The non-QUASI AtRTD2_19April2016.gtf carries the CDS lines. seed_cds_from_gtf.py reads the CDS start of each
transcript from it, converts it to transcript coordinates using the QUASI exon structure, and writes a CDS= token in the header.
Without this step reannotate_CDS.py finds no coding genes.

    Run seed_cds_from_gtf.py enriched.fa AtRTDv2_QUASI_19April2016.gtf AtRTD2_19April2016.gtf seeded.fa

3) AtRTD2 has many isoforms from coding genes lacking any associated CDS. 
On the other hand, some isoforms have shifted START codons (compare to others in the same gene) to avoid short ORFs (and PTCs), 
so it is important to analyze the transcriptome gene-by-gene to assign proper and common START codons to all the isoforms.

    Run reannotate_CDS.py

4) The final step is to assign PTCs to all the isoforms.

    Run detect_PTCs.py


The results can be easily visualized and further analyzed by loading the "_gene_summary.tsv" table. 


Typical run:

    python add_junction.py AtRTDv2_QUASI_19April2016.fa AtRTDv2_QUASI_19April2016.gtf step1.fa
    python seed_cds_from_gtf.py step1.fa AtRTDv2_QUASI_19April2016.gtf AtRTD2_19April2016.gtf step1_seeded.fa
    python reannotate_CDS.py step1_seeded.fa step2.fa
    python detect_ptcs.py step2.fa

Figure 1 analysis (figure1_analysis/)

Scripts, output tables and figures for the isoform/PTC analysis (panels A to D). The scripts read
AtRTD2_reannotated_CDS_transcript_summary.tsv from the working directory (a copy is in figure1_analysis/tables/).

- gene_ptc_status_pies.py: genes with 1 transcript vs genes with alternative isoforms, split into All PTC / All No_PTC / Mixed (panel A)
- gene_ptc_status_pies_2v3plus.py: the same comparison for genes with 2 isoforms vs 3 or more
- conditional_probability_PTC.py: conditional probability of a PTC isoform given non-PTC isoforms already present, and the reverse (panel B)
- productive_fraction_vs_isoforms_graphs_dodged.py: fraction of non-PTC isoforms vs isoforms per gene, and total vs non-PTC isoforms per bin (panels C and D)
- tables/: the numbers behind each plot, plus the per-gene table (table_per_gene_isoform_ptc_status.csv)
