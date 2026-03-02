AtRTD2 PTC Detection

PTC detection is based on STOP codons presence upstream of a splicing junction (>50 nt rule). 

We implemented a streamlined three-step pipeline in Python to perform the analyses.

1) Using the original AtRTD2_19April2016 fasta (.fa) and gtf (.gtf) files (available at https://ics.hutton.ac.uk/atRTD/) we generated an enriched FASTA with CDSs and Junctions.
To generate a file containing CDS and Junctions information...

    Run add_junction.py

3) AtRTD2 has many isoforms from coding genes lacking any associated CDS. 
On the other hand, some isoforms have shifted START codons (compare to others in the same gene) to avoid short ORFs (and PTCs), 
so it is important to analyze the transcriptome gene-by-gene to assign proper and common START codons to all the isoforms.

    Run reannotate_CDS.py

4) The final step is to assign PTCs to all the isoforms.

    Run detect_PTCs.py


The results can be easily visualized and further analyzed by loading the "_gene_summary.tsv" table. 
