
bam=ctrl2
ref=mm10_genome/mm10


# export PATH=/data/lemsaraa/references/cellranger-arc-2.0.2:$PATH
# cellranger-arc count --id=$bam \
#                        --reference=$ref \
#                        --libraries=$bam.csv \
#                        --localcores=16 \
#                        --localmem=64
                       
                       
export PATH=/data/lemsaraa/references/cellranger-7.1.0:$PATH                       
cellranger count --id=$bam --transcriptome $ref --fastqs ../scRNAseq/Mai13_TXM04_Michael/ --sample A006850200_172967 --chemistry ARC-v1 --star_parameters="--outFilterMismatchNmax=0"

samtools view -b -h $bam/outs/gex_possorted_bam.bam chr2 > $bam.bam
samtools sort $bam.bam -o sorted_$bam.bam
samtools index sorted_$bam.bam


# grep -w "Wt1os" mm10-2020-A-build/gencode.vM23.primary_assembly.annotation.gtf.filtered > Wt1os.gtf
# bedtools getfasta -fi Mus_musculus.GRCm38.dna.primary_assembly.fa -bed Wt1os.bed > Wt1os.fa
# awk '{print $1, $2, $3, $4 - 105076538, $5 - 105076538, $6, $7, $8, $9}' Wt1os.gtf > Wt1os.mod.gtf

# awk -F "\t" '{$4=$4-105076537; $5=$5-105076538; print}'  OFS="\t" Wt1os.gtf > Wt1os.mod.gtf

