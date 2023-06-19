#!/bin/bash

cov=$0

cd ./work || exit


mkdir cov_"${cov}" && cd cov_"${cov}" || exit


mkdir workdir && cd workdir || exit
echo "--------------------"
echo -e "Working for coverage ${cov}"
echo -e "--------------------\n"

mmseqs createdb --dbtype 1 -v 1 ../../input/Sarb_CDS.pep  ../../input/Sbay_CDS.pep  ../../input/Scer_NCBI_CDS.pep  ../../input/Skud_CDS.pep  ../../input/Smik_CDS.pep  ../../Spar_NCBI_CDS.pep DB 
mmseqs cluster --cluster-mode 0 --min-seq-id 0.7 -c "${cov}" --cov-mode 0 --remove-tmp-files 1 -v 1 DB clust tmp 


mmseqs createtsv -v 1 DB DB clust clust.tsv

mmseqs createseqfiledb -v 1 DB clust DB_clu_seq
mmseqs result2flat -v 1 DB DB DB_clu_seq clu_seq.fasta

echo -e "Clustering and flat generation done ... \n"

mv clu_seq.fasta ../
cd ..
echo -e "Generating fasta files for each cluster ... \n "
python3 /home/simon.herman/Bureau/Gits/ClusterGenes/flat2multi.py cov_"${cov}" "${cov}"

echo -e "Fasta generated for coverage cov_${cov} \n"


python3 /home/simon.herman/Bureau/Gits/ClusterGenes/parseMuscle.py "cov_${cov}" "${cov}"

cd ..
