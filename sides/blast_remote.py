from Bio import SeqIO
import os

cov = 0.5

# Function to find smallest sequence in a fasta file
def find_smallest_seq(fasta_file, cluster) -> SeqIO.SeqRecord:
    smallest_seq = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if smallest_seq is None or len(record.seq) < len(smallest_seq.seq):
            smallest_seq = record

    smallest_seq.id = cluster
    return smallest_seq

record_list = list()

cluster_dir = f"work/{cov}/clusters"
for entry in os.scandir(cluster_dir):

        size = entry.name

        if size != "size_0" and size != "size_1":

            cluster_mem = set()
            size_dir = f"work/{cov}/clusters/{size}"
            for cluster_name in os.listdir(size_dir): # get every cluster for a given cluster size

                print(cluster_name)
                record_list.append(find_smallest_seq(f"{size_dir}/{cluster_name}", cluster_name))
                

                
SeqIO.write(record_list, f"output/representants_clusters.faa", "fasta")                  
                
                