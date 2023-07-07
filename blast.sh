#!/bin/bash



blastp -query output/Nter_db.fasta -subject output/five_prime_db.fasta -outfmt '6 qseqid sseqid evalue qstart qend sstart send qseq sseq length' -out output/nter_five.tsv

blastp -query output/Cter_db.fasta -subject output/three_prime_db.fasta -outfmt '6 qseqid sseqid evalue qstart qend sstart send qseq sseq length' -out output/cter_three.tsv