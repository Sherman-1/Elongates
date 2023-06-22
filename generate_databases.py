import polars as pl
import gff3_parser
import yaml
from utils import multifasta_to_dict
from utils.chimeric_sequences import * 
import yaml 
from collections import defaultdict
from Bio import SeqIO
import subprocess

DB_PATH = "/home/simon.herman/Bureau/NC_DB"
THRESHOLD = 10
## Load data


species = yaml.safe_load(open('env.yaml'))["Species"] # Species will be passed as argument in the future

gff_dict = defaultdict(pl.DataFrame)
genome_dict = defaultdict(dict)

for specie in species:

    gff_dict[specie] = pl.from_pandas(gff3_parser.parse_gff3(f"/home/simon.herman/Bureau/SCER_NCBI_NEXT/{specie}.gff", parse_attributes = True))
    genome_dict[specie] = multifasta_to_dict(f"/home/simon.herman/Bureau/SCER_NCBI_NEXT/{specie}.fna", genome = True)

elongates = pl.read_csv("output/0.5_elongates.csv", has_header = True)
elongate_ids = set(
    elongates.filter(
    (pl.col('max_length_Nter') > THRESHOLD) | (pl.col("max_length_Cter") > THRESHOLD)
    )["num_cluster"].to_list()
    
)

filtered_data = elongate_data.filter(
    pl.col("num_cluster").is_in(elongate_ids)
)

infos = get_infos_for_chimeres(clusters = elongate_ids, filtered_data = filtered_data)

five_prime_db = list()
three_prime_db = list()
chimeres_db = list()


for key in infos.keys():
    
    for i, liste in enumerate(infos[key]):
        
        tmp = get_extended_UTRs(liste, gff_dict, genome_dict)

        if tmp == None:
            continue
        
        five_prime_db.extend(tmp["5utr"].values())
        three_prime_db.extend(tmp["3utr"].values())
        chimeres_db.extend(create_chimeric_sequences(tmp, liste, cds_dict))


SeqIO.write(five_prime_db, f"{DB_PATH}/five_prime_db.fasta", "fasta")
SeqIO.write(three_prime_db, f"{DB_PATH}/three_prime_db.fasta", "fasta")
SeqIO.write(chimeres_db, f"{DB_PATH}/chimeres_db.fasta", "fasta")