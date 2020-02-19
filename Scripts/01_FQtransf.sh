
## Intern variables: Defined by the user
PROJ_DIR="/drives/slave/jperales/Projects/SquatritoM_U251_sgMSH6_May17/HighTMZt1w_HighDMSOt1w";


# CASE
gunzip -c ${PROJ_DIR}"/../ajimenez_PCRseq_170503/AJ220317_10_1.fastq.gz" > ${PROJ_DIR}"/FastQ/U251_sgMSH6_hDDR_High_TMZ_t1w_E1.fastq";
gunzip -c ${PROJ_DIR}"/../ajimenez_PCRseq_170503/AJ220317_11_1.fastq.gz" > ${PROJ_DIR}"/FastQ/U251_sgMSH6_hDDR_High_TMZ_t1w_E2.fastq";
gunzip -c ${PROJ_DIR}"/../ajimenez_PCRseq_170503/AJ220317_12_1.fastq.gz" > ${PROJ_DIR}"/FastQ/U251_sgMSH6_hDDR_High_TMZ_t1w_E3.fastq";

# CONTROL
gunzip -c ${PROJ_DIR}"/../ajimenez_PCRseq_170503/AJ220317_7_1.fastq.gz" > ${PROJ_DIR}"/FastQ/U251_sgMSH6_hDDR_High_DMSO_t1w_E1.fastq";
gunzip -c ${PROJ_DIR}"/../ajimenez_PCRseq_170503/AJ220317_8_1.fastq.gz" > ${PROJ_DIR}"/FastQ/U251_sgMSH6_hDDR_High_DMSO_t1w_E2.fastq";
gunzip -c ${PROJ_DIR}"/../ajimenez_PCRseq_170503/AJ220317_9_1.fastq.gz" > ${PROJ_DIR}"/FastQ/U251_sgMSH6_hDDR_High_DMSO_t1w_E3.fastq";
