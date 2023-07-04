#!/bin/bash



blastp -query output/Nter_db.fasta -subject output/five_prime_db.fasta -outfmt 5 -out output/nter_five.xml 

blastp -query output/Cter_db.fasta -subject output/three_prime_db.fasta -outfmt 5 -out output/cter_three.xml