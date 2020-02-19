
## Intern variables: Defined by the user
module add NGS/bowtie/2.1.0

## Intern variables: Defined by the user
PROJ_DIR="/drives/slave/jperales/Projects/SquatritoM_U251_sgMSH6_May17/HighTMZt1w_HighDMSOt1w";
#--- DDR Library
prefixREF="${PROJ_DIR}/../REFERENCES/shRNA_libs/DDR_miRE/bowtie2_DDR_miRE_rc";

for sname in `ls -1 $PROJ_DIR"/FastQ/" | sed -e "s/\.fastq//g"`;do
echo -e "Trimming raw reads of $sname at $(date)";
java -jar /local/jperales/Soft/NGS/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 8 -phred33 -trimlog $PROJ_DIR"/Processed_FastQ/"$sname".log" $PROJ_DIR"/FastQ/"$sname".fastq" $PROJ_DIR"/Processed_FastQ/"$sname".fastq" ILLUMINACLIP:$PROJ_DIR"/../REFERENCES/Primers_and_adapters/shorter_adapter.txt":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:17 &> $PROJ_DIR"/Processed_FastQ/"$sname".log2";

echo -e "Alignment of ${sname}";
bowtie2 -x ${prefixREF} -p 12 --phred33 --norc -D 20 -R 3 -N 0 -L 17 -i S,1,0.50 -U "./Processed_FastQ/"${sname}".fastq" -S "./Aln/"${sname}".sam" &> "./Aln/"${sname}".alnlog";

echo -e "Counting read counts for ${sname}"; 
perl ${PROJ_DIR}/Build/tag_counter.pl ${PROJ_DIR}"/../REFERENCES/shRNA_libs/DDR_miRE/DDR_miRE.rev_compSPECIFIC_STEMs.fa" "./Aln/"${sname}".sam" > "./Count/"${sname}".counts.tsv";
done
