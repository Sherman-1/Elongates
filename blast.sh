#!/bin/bash



blastp -query output/Nter_db.fasta -subject output/five_prime_db.fasta -outfmt 5 -out output/nter_five.xml 