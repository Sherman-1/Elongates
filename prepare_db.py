import polars as pl
import pandas as pd
import gff3_parser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
import argparse


from utils import multifasta_to_dict
from utils.handle_UTRs import * 
THRESHOLD = 15

def prepare_db(cov):

    ## Load data
    species = yaml.safe_load(open('env.yaml'))["Species_order"]["Scer"] # Species will be passed as argument in the future
    elongates = pl.from_pandas(pd.read_csv(f"output/{cov}/{cov}_elongates.csv", low_memory = False))

    gff_dict = dict()
    genome_dict = dict()

    for specie in species:

        gff_dict[specie] = pl.from_pandas(gff3_parser.parse_gff3(f"input/{specie}.gff", parse_attributes = True, verbose = False))
        genome_dict[specie] = multifasta_to_dict(f"input/{specie}.fna", genome = True)

    infosDict = get_infos_for_UTRs(threshold=THRESHOLD, elongates = elongates)

    trans_five_prime_db = list()
    trans_three_prime_db = list()
    raw_five_prime_db = list()
    raw_three_prime_db = list()
    untranslated_nter = list()
    untranslated_cter = list()
    
    


    for cluster_id in infosDict.keys():

        for cds_infos in infosDict[cluster_id]: # list of dicts

            translated_UTRs, raw_UTRs, untranslated_elongates = get_extended_UTRs(cds_infos, gff_dict, genome_dict, cluster_id, cov)

            if translated_UTRs["5utr"]and raw_UTRs["5utr"]:

                trans_five_prime_db.extend(translated_UTRs["5utr"].values())
                raw_five_prime_db.append(raw_UTRs["5utr"])

            if translated_UTRs["3utr"] and raw_UTRs["3utr"]:

                trans_three_prime_db.extend(translated_UTRs["3utr"].values())
                raw_three_prime_db.append(raw_UTRs["3utr"])

            if untranslated_elongates["Nter"]:

                untranslated_nter.append(untranslated_elongates["Nter"])

            if untranslated_elongates["Cter"]:

                untranslated_cter.append(untranslated_elongates["Cter"])
    
   
    SeqIO.write(trans_five_prime_db, f"output/{cov}/five_prime_db.fasta", "fasta")
    SeqIO.write(trans_three_prime_db, f"output/{cov}/three_prime_db.fasta", "fasta")
    SeqIO.write(raw_five_prime_db, f"output/{cov}/five_prime_db_untrans.fasta", "fasta")
    SeqIO.write(raw_three_prime_db, f"output/{cov}/three_prime_db_untrans.fasta", "fasta")
    SeqIO.write(untranslated_nter, f"output/{cov}/Nter_untrans.fasta", "fasta")
    SeqIO.write(untranslated_cter, f"output/{cov}/Cter_untrans.fasta", "fasta")
    

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
        

    SeqIO.write(Nter_db, f"output/{cov}/Nter_db.fasta", "fasta")
    SeqIO.write(Cter_db, f"output/{cov}/Cter_db.fasta", "fasta")


        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--cov', type=str, help='Coverage')
    args = parser.parse_args()

    prepare_db(args.cov)

