#!/bin/bash

cov=$1

blastp -query output/${cov}/Nter_db.fasta -subject output/${cov}/five_prime_db.fasta -outfmt '6 qseqid sseqid evalue qstart qend sstart send qseq sseq length gapopen gaps' -out output/${cov}/nter_five.tsv
blastp -query output/${cov}/Nter_db.fasta -subject output/${cov}/five_prime_db.fasta -out output/${cov}/nter_five_raw.blast

blastp -query output/${cov}/Cter_db.fasta -subject output/${cov}/three_prime_db.fasta -outfmt '6 qseqid sseqid evalue qstart qend sstart send qseq sseq length gapopen gaps' -out output/${cov}/cter_three.tsv
blastp -query output/${cov}/Cter_db.fasta -subject output/${cov}/three_prime_db.fasta -out output/${cov}/cter_three_raw.blast