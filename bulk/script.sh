#!/bin/bash

SECONDS=0

# STEP 1: Run fastqc
# fastqc data/demo.fastq -o data/

# run fastp
# cat Sample_Names.tab | parallel -j 4 "fastp -i bastet2.ccg.uni-koeln.de/downloads/NGS_SF02_sfili_A006850289/A006850289_{}_L000_R1_001.fastq.gz -I bastet2.ccg.uni-koeln.de/downloads/NGS_SF02_sfili_A006850289/A006850289_{}_L000_R2_001.fastq.gz -o ./out_{}_R1.fastq.gz -O ./out_{}_R2.fastq.gz --html ./{}_fastp.html --json ./{}_fastp.json --report_title {} --adapter_fasta adapters.fa  2> {}_trimo.log"
# echo "fastp finished running!"

# # run trimmomatic to trim reads with poor quality
#  mkdir trimmedData
# cat Sample_Names.tab | parallel -j 4 "trimmomatic PE -threads 4 bastet2.ccg.uni-koeln.de/downloads/NGS_SF02_sfili_A006850289/A006850289_{}_L000_R1_001.fastq.gz bastet2.ccg.uni-koeln.de/downloads/NGS_SF02_sfili_A006850289/A006850289_{}_L000_R2_001.fastq.gz -baseout trimmedData/{} ILLUMINACLIP:adapters-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 -trimlog {}.log "
# echo "Trimmomatic finished running!"

# fastqc data/demo_trimmed.fastq -o data/

# mkdir ref
# cd ref
# wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
# wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
# wget https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
# cd ..

# STEP 2: Run STAR
# get the genome indices
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# cd ref
# gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
# gunzip Mus_musculus.GRCm39.109.gtf.gz
# cd ..
## set folder for the alignement based on modified and unmodified genome
fold_index=starIndex
fold_align=starAlign_rsem

# fold_index=starIndex_modified
# fold_align=starAlign_modified
 
 
# mkdir $fold_index
# mkdir $fold_align
# #create index
# STAR --runThreadN 32 --runMode genomeGenerate --genomeDir $fold_index --genomeFastaFiles ref/Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile ref/Mus_musculus.GRCm39.109.gtf --limitGenomeGenerateRAM=81727411808

# run alignment
ulimit -n 100000
cat Sample_Names.tab | parallel -j 2 "STAR --runThreadN 32 --genomeDir $fold_index  --readFilesIn out_{}_R1.fastq.gz out_{}_R2.fastq.gz --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix $fold_align/{} --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --readFilesCommand zcat"

# STEP 3: Run featureCounts - Quantification
# get gtf
# mkdir quants_exon
# featureCounts -T 8 -p -a ref/Mus_musculus.GRCm39.109.gtf -g 'exon_id' -o quants_exon/{}_featurecounts.txt $fold_align/*.bam
# echo "featureCounts finished running!"

# cat Sample_Names.tab | parallel -j 2 "samtools index $fold_align/{}Aligned.sortedByCoord.out.bam"

# cat Sample_Names.tab | parallel -j 2 "samtools view -b -h $fold_align/{}Aligned.sortedByCoord.out.bam 2 > starAlign/{}.chr2.bam"
# cat Sample_Names.tab | parallel -j 2 "samtools sort starAlign/{}.chr2.bam -o $fold_align/sorted_{}.chr2.bam"
# cat Sample_Names.tab | parallel -j 2 "samtools index $fold_align/sorted_{}.chr2.bam"


#### Mus musculus --------------------------------------------------------------

# mkdir salmon_out
# mkdir index_salmon
# mkdir index_salmon/Mus_musculus

# # # short
# salmon index \
#   -t ref/Mus_musculus.GRCm39.cdna.all.fa.gz \
#   -i index_salmon/Mus_musculus/short_index \
#   -k 23

# # long
# salmon index \
#   -t ref/Mus_musculus.GRCm38.cdna.all.fa.gz \
#   -i index_salmon/Mus_musculus/long_index \
#   -k 31


# cat Sample_Names.tab | parallel -j 2 "salmon quant -i index_salmon/Mus_musculus/short_index \
#     -l A \
#     -1 out_{}_R1.fastq.gz \
#     -2 out_{}_R2.fastq.gz \
#     -o salmon_out/{} \
#     --validateMappings \
#     --gcBias --seqBias \
#     --threads 4"
# mkdir rsem_out
# rsem-prepare-reference --gtf ref/Mus_musculus.GRCm39.109.gtf ref/Mus_musculus.GRCm39.dna.primary_assembly.fa rsem_index
cat Sample_Names.tab | parallel -j 2 "rsem-calculate-expression --no-bam-output --bam -p 8 $fold_align/{}Aligned.sortedByCoord.out.bam rsem_index/rsem_index rsem_out/{}"
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
