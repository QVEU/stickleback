#!/bin/bash -e
###$ -S /bin/bash
#$ -N threaded16_stickle
#$ -M Patrick.Dolan@nih.gov
#$ -m be
#$ -l h_vmem=15G
#$ -cwd
#$ -o OUT/
#$ -pe threaded 16

module purge
module load python

python /hpcdata/lvd_qve/QVEU_Code/StickleBack/stickleback.0.0.py /hpcdata/lvd_qve/Sequencing_Data/QVEU_Seq_0004_ONT-MinION_Spine-Library-QC/ev71top1-11_in/ev71top1-11_in/20220523_2115_MC-113212_FAT05465_e4da60b9/fastq_pass/merge.sam agcgggagaccggggtctctgagcg /hpcdata/lvd_qve/QVEU_Code/sequencing/template_fastas/puc19-ev71-twtainan1998_4643-bsmbi-and-bsai-free-deleted-1-annotations-1-7471.fasta 2000 3500
