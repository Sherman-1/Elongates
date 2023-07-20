#!/bin/bash

cov=$1

blastp -query output/${cov}/Nter_db.fasta -subject output/${cov}/five_prime_db.fasta -outfmt '6 qseqid sseqid evalue qstart qend sstart send qseq sseq length gapopen gaps' -out output/${cov}/nter_five_${cov}.tsv -evalue 100
blastp -query output/${cov}/Nter_db.fasta -subject output/${cov}/five_prime_db.fasta -out output/${cov}/nter_five_raw_${cov}.blast -evalue 100

blastp -query output/${cov}/Cter_db.fasta -subject output/${cov}/three_prime_db.fasta -outfmt '6 qseqid sseqid evalue qstart qend sstart send qseq sseq length gapopen gaps' -out output/${cov}/cter_three_${cov}.tsv -evalue 100
blastp -query output/${cov}/Cter_db.fasta -subject output/${cov}/three_prime_db.fasta -out output/${cov}/cter_three_raw_${cov}.blast -evalue 100