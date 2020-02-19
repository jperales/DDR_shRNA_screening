# DDR_shRNA_screening
Javier Perales-Patón, jperales@ext.cnio.es

## Bioinformatics analysis of shRNA screenings
Raw sequencing reads were trimmed using Trimmomatic software (v. 0.32) with default parameters, by removing the adaptor from the miR30 construct from the raw sequencing reads. Only processed reads with at least 17 nucleotide lengths were retrieved. Processed reads were aligned to the shRNA library reference using Bowtie2 (v. 2.1.0) [ref2], using the options for a very sensitive alignment with no mismatches following manual instructions (-D 20 -R 3 -N 0 -L 17 -i S,1,0.50). The read counts per hairpin were obtained with a custom Perl script that processes the Sequence Alignment/Map (SAM) file. Within-sample normalization was performed by calculating the read Counts Per Million per hairpin for each independent pool from the shRNA library design. Contrasts of hairpin representation upon treatment were obtained by the average log2-fold-change between treated and untreated cell lines from three replicates.

# References
[ref1] Anthony M. Bolger, Marc Lohse, Bjoern Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data, Bioinformatics, Volume 30, Issue 15, 1 August 2014, Pages 2114–2120

[ref2] Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

