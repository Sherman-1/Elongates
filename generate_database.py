import polars as pl
import gff3_parser
import yaml
from utils import multifasta_to_dict
from utils.chimeric_sequences import * 
import yaml 
from collections import defaultdict
from Bio import SeqIO
import subprocess

THRESHOLD = 15

## Load data
species = yaml.safe_load(open('env.yaml'))["Species"] # Species will be passed as argument in the future

gff_dict = dict()
genome_dict = dict()

for specie in species:

    gff_dict[specie] = pl.from_pandas(gff3_parser.parse_gff3(f"input/{specie}.gff", parse_attributes = True))
    genome_dict[specie] = multifasta_to_dict(f"input/{specie}.fna", genome = True)

elongates = pl.read_csv("output/0.5_elongates.csv", has_header = True)

infosList = get_infos_for_chimeres(threshold=THRESHOLD, elongates = elongates)

five_prime_db = list()
three_prime_db = list()

for cluster in infosList.keys():
    print(cluster)
    
    for index, infos in enumerate(infosList[cluster]):

        print(infos)
        
        
        tmp = get_extended_UTRs(infos, gff_dict, genome_dict, cluster)

        if tmp == None:
            continue
        
        five_prime_db.extend(tmp["5utr"].values())
        three_prime_db.extend(tmp["3utr"].values())

SeqIO.write(five_prime_db, f"output/five_prime_db.fasta", "fasta")
SeqIO.write(three_prime_db, f"output/three_prime_db.fasta", "fasta")