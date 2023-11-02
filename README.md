# KABOOM
KABOOM (pro**KA**ryote BUS**CO** phylogen**OM**ics)

Perform phylogenomic reconstruction using BUSCO genes

Usage
`KABOOM.py [-h] -i INPUT -o OUTPUT -t THREADS [-l LINEAGE] [-bp BUSCO_EXTRA_PARAMETERS]
[-psc PSC] [-ap ALIGNMENT_PROGRAM] [-ma MAFFT_ALGORITHM] [-am ANALYSIS_MODE]
[--trimal_strategy TRIMAL_STRATEGY] [--missing_character MISSING_CHARACTER]
[--phylogenetics_only] [-j JACKKNIFE] [-bb UFBOOTSTRAP]`

optional arguments:
`-h, --help show this help message and exit
-i INPUT, --input INPUT Input directory containing assembled genomes or completed BUSCO runs
-o OUTPUT, --output OUTPUT Output directory to store results
-t THREADS, --threads THREADS Number of threads to use
-l LINEAGE, --lineage LINEAGE Choose witch lineage to use by BUSCO. Default: bacteria_odb10)
-bp BUSCO_EXTRA_PARAMETERS, --busco_extra_parameters BUSCO_EXTRA_PARAMETERS Extra parameters for BUSCO. You can use this to restart or force a busco analysis. 
-psc PSC, --percent_single_copy PSC BUSCO presence cut-off. BUSCOs that are complete and single-copy in at least [-psc] percent of species will be included in the concatenated alignment [default=100.0]
-ap ALIGNMENT_PROGRAM, --alignment_program ALIGNMENT_PROGRAM Choose witch program to use to align
-ma MAFFT_ALGORITHM, --mafft_algorithm MAFFT_ALGORITHM Choose witch algorithm to use with mafft: fast, genafpair (E-INS-i), localpair (L-INS-i), globalpair (G-INS-i) [default=globalpair] 
-am ANALYSIS_MODE, --analysis_mode ANALYSIS_MODE Choose to perform phylogenetic reconstruction using nucleotide (DNA) sequences or aminoacids sequences (AA) or both (it will perform both analysis one after the other). Default: AA
--trimal_strategy TRIMAL_STRATEGY trimal trimming strategy (automated1, gappyout, strict, strictplus) [default=automated1]
--missing_character MISSING_CHARACTER Character to represent missing data [default='?']
--phylogenetics_only  It flags to only run the phylogenetic analysis. It expects to encounter BUSCO files
-j JACKKNIFE, --jackknife JACKKNIFE Choose the number of jackknife replicates to by done by IQ-TREE. Default: 1000
-bb UFBOOTSTRAP, --ufbootstrap UFBOOTSTRAP Choose the number of bootstrap replicates to by done by IQ-TREE. Default: 1000`

See an example of this work in:
