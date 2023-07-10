import polars as pl
import gff3_parser
from Bio import SeqIO
import yaml
from utils import multifasta_to_dict
from utils.chimeric_sequences import * 
THRESHOLD = 15

def prepare_db(cov):

    ## Load data
    species = yaml.safe_load(open('env.yaml'))["Species_order"]["Scer"] # Species will be passed as argument in the future
    elongates = pl.read_csv(f"output/{cov}_elongates.csv", has_header = True)

    gff_dict = dict()
    genome_dict = dict()

    for specie in species:

        gff_dict[specie] = pl.from_pandas(gff3_parser.parse_gff3(f"input/{specie}.gff", parse_attributes = True, verbose = False))
        genome_dict[specie] = multifasta_to_dict(f"input/{specie}.fna", genome = True)

    infosDict = get_infos_for_UTRs(threshold=THRESHOLD, elongates = elongates)

    five_prime_db = list()
    three_prime_db = list()


    for cluster_id in infosDict.keys():

        for cds_infos in infosDict[cluster_id]:

            UTR_dict = get_extended_UTRs(cds_infos, gff_dict, genome_dict, cluster_id)

            if UTR_dict == None:

                continue

            five_prime_db.extend(UTR_dict["5utr"].values() if UTR_dict["5utr"] != None else [])
            three_prime_db.extend(UTR_dict["3utr"].values() if UTR_dict["3utr"] != None else [])

    SeqIO.write(five_prime_db, f"output/five_prime_db.fasta", "fasta")
    SeqIO.write(three_prime_db, f"output/three_prime_db.fasta", "fasta")


    # Get elongated sequences
    Nter_db = list()
    Cter_db = list()

    for row in elongates.filter(

        (pl.col("Nter_elongate_length") >= THRESHOLD) 

                                ).iter_rows(named=True):
        
        Nter_db.append(SeqRecord(seq = Seq(row["Nter_elongate"].replace("-","")), id = f"{row['seq_id']}-{row['cluster_name']}", description = ""))


    for row in elongates.filter(

        (pl.col("Cter_elongate_length") >= THRESHOLD) 

                                ).iter_rows(named=True):
        
        Cter_db.append(SeqRecord(seq = Seq(row["Cter_elongate"].replace("-","")), id = f"{row['seq_id']}-{row['cluster_name']}", description = ""))
        

    SeqIO.write(Nter_db, f"output/Nter_db.fasta", "fasta")
    SeqIO.write(Cter_db, f"output/Cter_db.fasta", "fasta")


        



