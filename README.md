# KABOOM
KABOOM (pro**KA**ryote BUS**CO** phylogen**OM**ics)

Perform phylogenomic reconstruction using BUSCO genes.

Usage

```
usage: KABOOM.py [-h] -i INPUT -o OUTPUT [-t THREADS] [-m MODE] [-l LINEAGE]
                 [-bp BUSCO_EXTRA_PARAMETERS] [-psc PSC]
                 [-as ALIGNMENT_SOFTWARE] [-ma MAFFT_ALGORITHM]
                 [-am ANALYSIS_MODE] [--trimal_strategy TRIMAL_STRATEGY]
                 [--missing_character MISSING_CHARACTER] [-j JACKKNIFE]
                 [-bb UFBOOTSTRAP]

```

optional arguments:

```
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory containing assembled genomes or
                        completed BUSCO runs
  -o OUTPUT, --output OUTPUT
                        Output directory to store results
  -t THREADS, --threads THREADS
                        Number of threads to use
  -m MODE, --mode MODE  Choose witch mode to run: all, concatenate, phylogeny
                        [default=all]
  -l LINEAGE, --lineage LINEAGE
                        Choose witch lineage to use by BUSCO. Default:
                        bacteria_odb10)
  -bp BUSCO_EXTRA_PARAMETERS, --busco_extra_parameters BUSCO_EXTRA_PARAMETERS
                        Extra parameters for BUSCO. You can use this to
                        restart or force a busco analysis.
  -psc PSC, --percent_single_copy PSC
                        BUSCO presence cut-off. BUSCOs that are complete and
                        single-copy in at least [-psc] percent of species will
                        be included in the concatenated alignment
                        [default=100.0]
  -as ALIGNMENT_SOFTWARE, --alignment_software ALIGNMENT_SOFTWARE
                        Choose witch software to use to align (muscle, mafft)
  -ma MAFFT_ALGORITHM, --mafft_algorithm MAFFT_ALGORITHM
                        Choose witch algorithm to use with mafft: fast,
                        genafpair (E-INS-i), localpair (L-INS-i), globalpair
                        (G-INS-i) [default=globalpair]
  -am ANALYSIS_MODE, --analysis_mode ANALYSIS_MODE
                        Choose to perform phylogenetic reconstruction using
                        DNA sequences (nucleotides), Amino Acids sequences
                        (aminoacids) or both (it will perform both analysis
                        one after the other). Default: aminoacids
  --trimal_strategy TRIMAL_STRATEGY
                        trimal trimming strategy (automated1, gappyout,
                        strict, strictplus) [default=automated1]
  --missing_character MISSING_CHARACTER
                        Character to represent missing data [default='?']
  -j JACKKNIFE, --jackknife JACKKNIFE
                        Choose the number of jackknife replicates to by done
                        by IQ-TREE. Default: 1000
  -bb UFBOOTSTRAP, --ufbootstrap UFBOOTSTRAP
                        Choose the number of bootstrap replicates to by done
                        by IQ-TREE. Default: 1000
```

See an example of the usage of this script in:
