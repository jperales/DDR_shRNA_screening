In general terms, shRNA libraries were brought as Excel files. Each pool were indicated with colors. Each specific shRNA sequence were indicated as 2 columns: the id and the duplex sequence. The duplex sequence contained 7p_primer-stem-loop-stem-5p_primer, while we expected to sequence the reverse complement of TruSeq_ADAPTOR-loop-stem...

# Rules
	- Three columns were copy-pasted into txt files from Excel as TSV format. These files ennamed [lib name].original.tsv.
	- Identifiers with a dot '.' where transform to '_'.
	- We extract the antisense specific of the stem 3' by using REXPR and Perl.
