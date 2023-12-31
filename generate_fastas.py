
import os, sys
from Bio import SeqIO 
import csv


'''
Use outputs of MMseqs2 clustering to generate fasta files for each clusters + size statistics

Args 
    directory automatically generated by main.sh
'''

cov = sys.argv[1]

reads = list(SeqIO.parse("tmp/clu_seq.fasta", format = "fasta"))

i = 0
buffer = []
sizes = []
dic = {}
for read in reads:

    if read.seq == "":

        cluster_size = len(buffer)

        if not os.path.isdir(f"clusters/size_{cluster_size}"):

            os.mkdir(f"clusters/size_{cluster_size}")
        
        handle = open(f"clusters/size_{cluster_size}/cluster_n{i}","w") 
        SeqIO.write(buffer, handle, format = "fasta")
        handle.close()

        buffer = []
        i = i+1
        
        if cluster_size in dic.keys():
            dic[cluster_size] = dic[cluster_size] + 1
        else:
            dic[cluster_size] = 1

    else:

        buffer.append(read)

# Last sequences still in buffer
cluster_size = len(buffer)

if not os.path.isdir(f"clusters/size_{cluster_size}"):

    os.mkdir(f"clusters/size_{cluster_size}")

handle = open(f"clusters/size_{cluster_size}/cluster_n{i}","w") 
SeqIO.write(buffer, handle, format = "fasta")
handle.close()

if cluster_size in dic.keys():
    dic[cluster_size] = dic[cluster_size] + 1
else:
    dic[cluster_size] = 1

# Write size statistics
with open(f'../../output/{cov}/yeasts_stats.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Cluster size', 'count'])  # writing headers
    for key, value in dic.items():
        writer.writerow([key, value])



