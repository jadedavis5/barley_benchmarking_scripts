#!/bin/bash -l

#SBATCH --job-name=comp
#SBATCH --account=fl3
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=10:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --export=ALL
#SBATCH --mem=30GB

module load singularity/4.1.0-slurm
module load blast/2.12.0--pl5262h3289130_0
module load interproscan/5.56-89.0

agat=/software/projects/fl3/jdavis/setonix/containers/agat_1.0.0--pl5321hdfd78af_0.sif
gffcompare=/software/projects/fl3/jdavis/setonix/containers/gffcompare_0.11.2.sif
isoquant=/software/projects/fl3/jdavis/.nextflow_singularity/depot.galaxyproject.org-singularity-isoquant%3A3.4.2--hdfd78af_0.img
gffread=/software/projects/fl3/jdavis/setonix/containers/gffread_0.11.7.sif

#Get BaRT files
wget https://ics.hutton.ac.uk/barleyrtd/downloads/BaRT2v18_transfix_pep.fasta.gz
gunzip BaRT2v18_transfix_pep.fasta.gz
bart_ref=/scratch/fl3/jdavis/FULL_TESTS/gffcompare/BaRT/BaRT2v18_transfix_pep.fasta

gff_id=STref
ref_annotation=../../ref_files/RGT_Planet.gff
gff=../../completed_gtf/outputAnnotation_STref.gff3
genome=../RGT_Planet_pseudomolecule_v1.fasta

#Find transcripts not in RGT gtf
singularity run $gffcompare gffcompare -R -r $annotation $gff
cat ../completed_gtf/gffcmp.outputAnnotation_${gff_id}.gff3.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq > ${gff_id}_novel.txt
awk '{print $0; print $0}' ${gff_id}_novel.txt > ${gff_id}_duplicated.txt
awk 'NR%2==1 {print "ID=" $0 ";"} NR%2==0 {print "Parent=" $0}' ${gff_id}_duplicated.txt > ${gff_id}_modified.txt
rm ${gff_id}_duplicated.txt
awk '{print "\\b" $0 "\\b"}' ${gff_id}_modified.txt | grep -E -f - $gff > ${gff_id}_transcripts.gff

#Convert to fa
singularity run $agat agat_sp_extract_sequences.pl --gff ${gff_id}_transcripts.gff -f $genome -o ${gff_id}-proteins.fa --merge

#Run blast
makeblastdb -in $bart_ref -title "BaRT-proteins" -dbtype prot

blastx -html -num_threads 256 -query ${gff_id}-proteins.fa -db $bart_ref \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp" > ${gff_id}-blastALL.tsv

#Identify transcripts which are also seen in BaRT
awk '{if ($3 > 95 && $13 > 95) print $0}' ${gff_id}-blastALL.tsv > ${gff_id}-blast-HC.tsv
cat ${gff_id}-blast-HC.tsv | cut -f1 | sort | uniq > ${gff_id}_BaRTnovel.txt

#Find all novel transcripts
grep -w -vFf ${gff_id}_BaRTnovel.txt ${gff_id}_novel.txt > ${gff_id}_NOVEL_FINAL.txt


#Generate gff file with only novel transcripts
awk '{print $0; print $0}' ${gff_id}_NOVEL_FINAL.txt > ${gff_id}_NOVEL_FINAL_duplicated.txt
awk 'NR%2==1 {print "ID=" $0 ";"} NR%2==0 {print "Parent=" $0}' ${gff_id}_NOVEL_FINAL_duplicated.txt > ${gff_id}_NOVEL_FINAL_modified.txt
awk '{print "\\b" $0 "\\b"}' ${gff_id}_NOVEL_FINAL_modified.txt | grep -E -f - $gff > ${gff_id}_NOVEL_FINAL_transcripts.gff


############ FIND HIGHLY EXPRESSED TRANSCRIPTS #########################
#Turn into GTF file for IsoQuant
singularity run $gffread gffread RGT_Planet_NBS.gff -T -o RGT_Planet_NBS.gtf


#Find transcripts with more than 50 reads mapping
singularity run $isoquant isoquant.py --reference $genome --genedb RGT_Planet.gtf \
--data_type nanopore -o ${gff_id}_ISOQUANT --no_model_construction -t 256 --matching_strategy precise \
--bam ../rna_011_ptt_f32_nb29_RGT_Planet_pseudomolecule_v1_aln_sorted.bam ../rna_011_ptt_f32_control_RGT_Planet_pseudomolecule_v1_aln_sorted.bam

awk '$2 > 50 && $3 > 50 {print $1}' ${gff_id}_ISOQUANT/OUT/OUT.transcript_grouped_counts.tsv > ${gff_id}_highlyExpressed_novel.txt

############ FUNCTIONAL CHARACTERISAITON #########################

#Convert to fa
singularity run $agat agat_sp_extract_sequences.pl --gff ${gff_id}_NOVEL_FINAL_transcripts.gff -f $genome -t cds -p --cfs --cis -o ${gff_id}-novel_proteins.fa

#Run interpro
mkdir input output temp
mv ${gff_id}-novel_proteins.fa input

curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/alt/interproscan-data-5.69-101.0.tar.gz
tar -pxzf interproscan-data-5.69-101.0.tar.gz


singularity exec \
    -B $PWD/interproscan-5.69-101.0/data:/opt/interproscan/data \
    -B $PWD/input:/input \
    -B $PWD/temp:/temp \
    -B $PWD/output:/output \
    interproscan_latest.sif \
    /opt/interproscan/interproscan.sh \
    --input /input/${gff_id}-novel_proteins.fa \ \
    --disable-precalc \
    -goterms -pa \
    --output-dir /output \
    --tempdir /temp

#Convert to space separated count format for Revigo
cat output/${gff_id}-novel_proteins.fa.tsv | cut -f 14 | sort | tr '|' '\n' \
| sort | uniq -c | sort -rn | sed 's/^[[:space:]]*\([0-9]*\)[[:space:]]*GO:\([0-9]*\)(.*)/GO:\2 \1/' > ${gff_id}_RevigoCount.tsv

