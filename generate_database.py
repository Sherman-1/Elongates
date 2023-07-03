import polars as pl
import gff3_parser
from Bio import SeqIO
import yaml
from utils import multifasta_to_dict
 

## Load data
species = yaml.safe_load(open('env.yaml'))["Species"] # Species will be passed as argument in the future
elongates = pl.read_csv("output/0.5_elongates.csv", has_header = True)

gff_dict = dict()
genome_dict = dict()

for specie in species:

    gff_dict[specie] = pl.from_pandas(gff3_parser.parse_gff3(f"input/{specie}.gff", parse_attributes = True, verbose = False))
    genome_dict[specie] = multifasta_to_dict(f"input/{specie}.fna", genome = True)


from utils.chimeric_sequences import * 
THRESHOLD = 15


infosDict = get_infos_for_UTRs(threshold=THRESHOLD, elongates = elongates)

five_prime_db = list()
three_prime_db = list()


for cluster_id in infosDict.keys():

    for sequence in infosDict[cluster_id]:

        UTR_dict = get_extended_UTRs(sequence, gff_dict, genome_dict, cluster_id)

        if UTR_dict == None:

            continue

        five_prime_db.extend(UTR_dict["5utr"].values())
        three_prime_db.extend(UTR_dict["3utr"].values())


elongates_db = list()





        
SeqIO.write(five_prime_db, f"output/five_prime_db.fasta", "fasta")
SeqIO.write(three_prime_db, f"output/three_prime_db.fasta", "fasta")


