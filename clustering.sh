#!/bin/bash

cov=$1

cd work || exit

mkdir cov_${cov}

mkdir cov_${cov}/tmp 
mkdir cov_${cov}/clusters

echo "-------------------------"
echo -e "Working for coverage ${cov}"
echo -e "-------------------------\n\n"

cd ../input || exit

mmseqs createdb --dbtype 1 -v 1 Sarb_CDS.pep  Sbay_CDS.pep  Scer_NCBI_CDS.pep  Skud_CDS.pep  Smik_CDS.pep Spar_NCBI_CDS.pep ../work/cov_${cov}/tmp/DB 

cd ../work/cov_${cov}/tmp || exit

echo " Clustering . . . "
mmseqs cluster --cluster-mode 0 --min-seq-id 0.7 -c 0.5 --cov-mode 0 --remove-tmp-files 1 -v 1 DB clust .

echo " Parsing files . . . "
mmseqs createtsv -v 1 DB DB clust clust.tsv
mmseqs createseqfiledb -v 1 DB clust DB_clu_seq
mmseqs result2flat -v 1 DB DB DB_clu_seq clu_seq.fasta

cd .. || exit

python3 ../../flat2multi.py ${cov}

echo -e "\n\n\n-------------------------"
echo -e "Done for coverage ${cov}"
echo -e "-------------------------\n\n"