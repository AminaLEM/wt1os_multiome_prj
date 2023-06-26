# Set the path to the Cell Ranger output directory
dir_path="bams/ctrl1"

# Set the path to the reference transcriptome FASTA and index files
ref_fasta="Mus_musculus.GRCm38.cdna.all.fa"
ref_idx="Mus_musculus.GRCm38.cdna.all.fa.fai"

# Use BUStools to create a BUS file
# bustools index -i $ref_fasta.idx $ref_fasta
bustools correct -w $ref_idx -p ${dir_path}/outs/gex_possorted_bam.bam -o ${dir_path}/corrected_reads.bus --umi-dedup
bustools sort -t 4 -o ${dir_path}/sorted_reads.bus ${dir_path}/corrected_reads.bus
bustools count -o ${dir_path}/transcript_counts --genecounts -g $ref_fasta -e matrix.ec -t /path/to/refdata-gex-mm10-ensembl/transcripts.txt -b ${dir_path}/sorted_reads.bus