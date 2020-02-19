#!/bin/bash


# The .original.tsv file was created by copy-paste from Excel into a txt file.

# Get .original.fa : Fasta sequences from the excel. It is useful to discern what sequences are
#       primers and loop (common) from specific sequences (specific shRNA sequence)
# For that purpose, we will use jalview and Percentaje of identity: only base position with 100% PI will be primers and adaptors, the same with the loop.

awk -F"\t" '{print ">"$1;print $2;}' DDR_miRE.original.tsv > DDR_miRE.original.fa

# The seq to eliminate are 3'->5':
# The starting 3' seq:
#       TGCTGTTGACAGTGAGCG
#       TGCTGTTGACAGTGAGCG
# Follow of one of the stem:
#       N{22}
# Follow of the loop:
#       TAGTGAAGCCACAGATGTA
#       TAGTGAAGCCACAGATGTA
# The 5' end (Backbone):
#       TGCCTACTGCCTCGGA
#       TGCCTACTGCCTCGGA

# so,
awk -F"\t" 'BEGIN{OFS="\t";}{gsub(/^TGCTGTTGACAGTGAGCG......................TAGTGAAGCCACAGATGTA/,"",$2); gsub(/TGCCTACTGCCTCGGA$/,"",$2); print $1,$2,$3;}' DDR_miRE.original.tsv

# Only to check that everything was right recognized:
awk -F"\t" 'BEGIN{OFS="\t";}{gsub(/^TGCTGTTGACAGTGAGCG......................TAGTGAAGCCACAGATGTA/,"",$2); gsub(/TGCCTACTGCCTCGGA$/,"",$2); print length($2);}' DDR_miRE.original.tsv | sort -u
# 22nt

# All right, we run it:
awk -F"\t" 'BEGIN{OFS="\t";}{gsub(/\./,"_",$1); gsub(/^TGCTGTTGACAGTGAGCG......................TAGTGAAGCCACAGATGTA/,"",$2); gsub(/TGCCTACTGCCTCGGA$/,"",$2); print $1,$2,$3;}' DDR_miRE.original.tsv > DDR_miRE.SPECIFIC_STEMs.tsv

# Now we have to get the reverse complement seq of these specific stems
perl -e 'my $fl="./DDR_miRE.SPECIFIC_STEMs.tsv"; %mapping_nt=("A"=>"T","T"=>"A","C"=>"G","G"=>"C"); open(IN,"$fl"); while(my $line=<IN>){ chomp($line); my @fields=split("\t",$line); my @seq=split("",$fields[1]); my @comp_seq=(); for(my $i=0;$i<scalar(@seq);$i++){ $comp_seq[$i]=$mapping_nt{$seq[$i]};} my @rev_comp=reverse(@comp_seq); print $fields[0],"\t",join("",@rev_comp),"\t",$fields[2],"\n";} close(IN);' > DDR_miRE.rev_compSPECIFIC_STEMs.tsv

# Get the uniq sequences: There are some hairpins duplicated across pools.
cut -f 1,2 DDR_miRE.rev_compSPECIFIC_STEMs.tsv | sort -u > DDR_miRE.uniq_rev_compSPECIFIC_STEMs.tsv
# Get the fasta of that
awk -F"\t" '{print ">"$1;print $2;}' DDR_miRE.uniq_rev_compSPECIFIC_STEMs.tsv > DDR_miRE.rev_compSPECIFIC_STEMs.fa

# Bowtie2
module add NGS/Bowtie/2.1.0
bowtie2-build  DDR_miRE.rev_compSPECIFIC_STEMs.fa bowtie2_DDR_miRE_rc
