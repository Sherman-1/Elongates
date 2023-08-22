#!/bin/bash

cov=$1
species=$(yq '.Species_order.Scer_for_bash' env.yaml)

mkdir -p work || exit
mkdir -p output || exit
mkdir -p output/${cov} || exit
cd work || exit

mkdir -p ${cov}

mkdir -p ${cov}/tmp 
mkdir -p ${cov}/clusters


echo "-------------------------"
echo -e "Working for coverage ${cov}"
echo -e "-------------------------\n\n"

cd ../input || exit

mmseqs createdb --dbtype 1 -v 1 $(echo $species | sed 's/\([a-zA-Z_]*\)/\1_CDS_corr.pep/g') ../work/${cov}/tmp/DB

cd ../work/${cov}/tmp || exit

echo " Clustering . . . "
mmseqs cluster --cluster-mode 1 --min-seq-id 0.7 -c ${cov} --cov-mode 0 --remove-tmp-files 1 -v 1 DB clust .

echo " Parsing files . . . "
mmseqs createtsv -v 1 DB DB clust clust.tsv
mmseqs createseqfiledb -v 1 DB clust DB_clu_seq
mmseqs result2flat -v 1 DB DB DB_clu_seq clu_seq.fasta

cd .. || exit

python3 ../../generate_fastas.py ${cov}

echo -e "\n\n\n-------------------------"
echo -e "Done for coverage ${cov}"
echo -e "-------------------------\n\n"