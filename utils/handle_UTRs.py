"""
chimeric_sequences.py

Author: HERMAN Simon
Date: 18/04/23
Version: 0.01

These functions facilitate the generation of chimeric sequences by combining the coding DNA sequence (CDS)
with different translated versions of 5' UTR and 3' UTR sequences. 

Functions included in this script are:

    translate_frames(dna_sequence, specie, seq_id, length, utr): Translates DNA sequence into protein for
    3 different reading frames and returns a dictionary of translated sequences for each frame.

    get_infos_for_chimeres(cluster, filtered_data): Retrieves information for chimeric sequences from the
    filtered dataset.

    get_extended_CDS(tuple_, gff_dict, genome_dict, sequence=None): Retrieves the extended coding sequence (CDS)
    of a specified gene from a given species.

    create_chimeric_sequences(utr5_versions, utr3_versions, specie, seq_id, length, cds_dict): Constructs chimeric
    sequences by combining the CDS with different versions of 5' UTR and 3' UTR sequences.

Please refer to the function documentation for additional details and usage examples.

Note: Make sure to properly handle exceptions and provide appropriate error messages to ensure
smooth functioning of the script.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import polars as pl
import csv
from utils.files import write_dicts_to_csv

__all__ = ["translate_frames", "get_infos_for_UTRs", "get_extended_UTRs"]

def translate_frames(dna_sequence, specie, seq_id, length, utr, cluster) -> dict():

    """
    Translates a given DNA sequence into protein sequences for the three different reading frames. Frame 1 is
    complemented for a multiplicity of 3, and the translated sequences are returned as a dictionary of SeqRecord
    objects.

    Parameters:
    --------
    dna_sequence (str or Bio.Seq): The DNA sequence to be translated. It can be a string or a Bio.Seq object.
    specie (str): The name of the species the DNA sequence belongs to.
    seq_id (str): The identifier of the input DNA sequence.
    length (int): The length of the untranslated region (UTR) in the input DNA sequence.
    utr (str): The untranslated region (UTR) type, either 'five_prime' or 'three_prime'.

    Returns:
    --------
    result (dict): A dict containing the translated protein sequences for each reading frame as SeqRecord objects. 
    """

    result = dict()
    dna_seq = dna_sequence if isinstance(dna_sequence, Seq) else Seq(dna_sequence)


    # Translate DNA sequence into protein for 3 different reading frames
    # We divide length by 3 because the length of the UTR is given in nucleotides

    result["frame_0"] = SeqRecord(seq = Seq(dna_seq[:length - length % 3].translate()), 
                                id = f"{seq_id}-{cluster}-{utr}-f0-{int(length/3)}", 
                                description = "") # last modif : tried to remove description = specie

    result["frame_1"] = SeqRecord(seq = Seq(dna_seq[1:length - (length - 1) % 3].translate()), 
                                id = f"{seq_id}-{cluster}-{utr}-f1-{int(length/3)}", 
                                description = "")

    result["frame_2"] = SeqRecord(seq = Seq(dna_seq[2:length - (length - 2) % 3].translate()), 
                                id = f"{seq_id}-{cluster}-{utr}-f2-{int(length/3)}", 
                                description = "")

    return result

def get_infos_for_UTRs(threshold, elongates):

    """
    Retrieves information for UTRs extraction from the filtered dataset. Input clusters correspond 
    to clusters where max Nter or Cter elongation is greater than fixed threshold defined at the beginning
    of the pipeline. 

    Parameters:
    --------
    threshold (int): The threshold value from which we consider a cluster to be elongated.
    elongates (pl.DataFrame): A polars DataFrame containing the extracted data for all clusters, regardless of their elongation status.

    Returns:
    --------
    out (dict): A dict of dicts

    """ 

    # Filter clusters to get only the ones above threshold ( 15 by default )
    elongates_filtered = elongates.filter(
        (pl.col('max_Nter') >= threshold) | (pl.col("max_Cter") >= threshold)
        )[["cluster_name", "species", "seq_id", "max_Nter", 
           "max_Cter","Nter_elongate","Cter_elongate",
           "is_max_Nter","is_max_Cter","Nter_nb_aa","Cter_nb_aa",]]

    out = {}
    
    # A bit dirty but it works for now
    def add_to_dict(row):
        
        cluster_name = row[0]
        info = {
            "species": row[1],
            "seq_id": row[2],
            "five_length": row[3]*3, # We multiply by 3 to get the length in nucleotides
            "three_length": row[4]*3, # Longest elongate is not used for UTR extraction
            "is_max_Nter": row[7],
            "is_max_Cter": row[8],
            "Nter_nb_aa": row[9],
            "Cter_nb_aa": row[10] 
        }
        
        if cluster_name in out:
            out[cluster_name].append(info)
        else:
            out[cluster_name] = [info]
            
        return tuple(row)
            
            
    elongates_filtered.apply(add_to_dict)

    return out

def get_extended_UTRs(cds_infos, gff_dict, genome_dict, cluster, cov):

    """
    Compute the translations of UTR regions for a given CDS. Each UTR region ( 5' or 3' ) is translated
    into protein for the three different reading frames. The translated sequences are returned as a dictionary

    Parameters:
    --------
    cds_infos (dict): dict containing the specie, sequence id and length of the CDS for which the UTRs are to be
                    translated.
    gff_dict (dict): Dictionary containing the gff files for each specie.
    genome_dict (dict): Dictionary containing the FASTA files for each specie.
    
    Returns:
    --------
    result_dict (dict): Dict of dict containing the translated UTRs for each reading frame. 
    
    Structure of result_dict:
    -------- 
    {"five_prime": {"frame_1": SeqRecord, "frame_2": SeqRecord, "frame_3": SeqRecord},
    "three_prime": {"frame_1": SeqRecord, "frame_2": SeqRecord, "frame_3": SeqRecord}}
    """

    # Initilize
    specie = cds_infos["species"]
    seq_id = cds_infos["seq_id"]
    is_max_Nter = cds_infos["is_max_Nter"]
    is_max_Cter = cds_infos["is_max_Cter"]

    if is_max_Nter:
        FIVE_LENGTH = cds_infos["five_length"] 
    else:
        FIVE_LENGTH = cds_infos["five_length"]*2

    if is_max_Cter:
        THREE_LENGTH = cds_infos["three_length"]
    else:
        THREE_LENGTH = cds_infos["three_length"]*2
    
    coordinates = [] # End and start coordinates of each CDS features
    raw_UTRs = {} # Dictionary to store the raw UTRs
    translated_UTRs = {} # Dictionary to store the translated UTRs
    untranslated_elongates = {}

    if specie == "unknown":

        return None

    # Filter gff file to get only CDS for the given specie and gene
    
    gff = gff_dict[specie].filter(

        (pl.col("Type") == "CDS") & ((pl.col("Name") == seq_id) | (pl.col("Parent") == seq_id))

    )[["Start","End","Strand","Seqid"]] # Keep only necessary columns
    
    # Store datas necessary to compute the elongate sequence

    strand = gff[0]["Strand"].to_list()[0] # + or -
    strand_id = gff[0]["Seqid"].to_list()[0] # chromosome or scaffold id

    for row in gff.iter_rows(named=True): # Named = True to iter with column names

        coordinates.append(sorted((int(row['Start'])-1, int(row['End'])-1))) # -1 for python indexing
        
    coordinates = sorted(coordinates, key=lambda x: x[0]) # Sort coordinates by start position

    # coordinates[0] = (start, end) of the first exon
    # coordinates[-1] = (start, end) of the last exon

    if strand == "+":

        if not is_max_Nter:

            start_5 = coordinates[0][0]-FIVE_LENGTH if coordinates[0][0]-FIVE_LENGTH >= 0 else 0 # Get the start position of the 5' UTR
            
            five_prime = genome_dict[specie][strand_id]["seq"][

                start_5:coordinates[0][0]

            ] # Get the 5'UTR corresponding to the non-coding region of the elongate in the cluster

        else:

        
            five_prime = genome_dict[specie][strand_id]["seq"][

                coordinates[0][0]:coordinates[0][0]+FIVE_LENGTH

            ] # Get the 5' sequence corresponding to the actual, detected, elongate 

        
        if not is_max_Cter:

            end_3 = coordinates[-1][1]+1+THREE_LENGTH if coordinates[-1][1]+1+THREE_LENGTH <= genome_dict[specie][strand_id]["len"] else genome_dict[specie][strand_id]["len"] # Get the end position of the 3' UTR
            # +1 for -1,1 because GFF points to the last nucleotide of the stop codon

            three_prime = genome_dict[specie][strand_id]["seq"][
                coordinates[-1][1]+1:end_3
                ] # Get the 3' sequence

        else:

            three_prime = genome_dict[specie][strand_id]["seq"][
                coordinates[-1][1]-THREE_LENGTH:coordinates[-1][1]+1
            ]

        with open(f"output/{cov}/{specie}_to_remove.gff", "a") as f:
                
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([

                strand_id, "to_remove", "gene", coordinates[0][0] - len(five_prime) + 1, coordinates[-1][1] + len(three_prime) + 1, ".", "+", ".", f"ID={seq_id};Species={specie};Cluster={cluster}"

            ])

    
    # Reverse complement if the strand is negative, don't forget to reverse the coordinates
    if strand == "-":

        if not is_max_Nter:

            end_5 = coordinates[-1][1]+1+FIVE_LENGTH if coordinates[-1][1]+1+FIVE_LENGTH <= genome_dict[specie][strand_id]["len"] else genome_dict[specie][strand_id]["len"] # Get the start position of the 5' UTR

            five_prime = genome_dict[specie][strand_id]["seq"][
                coordinates[-1][1]+1:end_5
            ].reverse_complement() # Get the 5' sequence
        
        else:

            five_prime = genome_dict[specie][strand_id]["seq"][
                coordinates[-1][1]-FIVE_LENGTH:coordinates[-1][1]
            ].reverse_complement()

        if not is_max_Cter:
        
            start_3 = coordinates[0][0]-THREE_LENGTH if coordinates[0][0]-THREE_LENGTH >= 0 else 0 # Get the end position of the 3' UTR

            three_prime = genome_dict[specie][strand_id]["seq"][
                start_3:coordinates[0][0]
            ].reverse_complement() # Get the 3' sequence
            
        else:

            three_prime = genome_dict[specie][strand_id]["seq"][
                coordinates[0][0]:coordinates[0][0]+THREE_LENGTH+1
            ].reverse_complement()
        

        with open(f"output/{cov}/{specie}_to_remove.gff", "a", newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([

                strand_id, "to_remove", "gene", coordinates[0][0]-len(three_prime)+1, coordinates[-1][1] + len(five_prime) +1, ".", "-", ".", f"ID={seq_id};Species={specie};Cluster={cluster}"

            ])
        

    if len(five_prime) != FIVE_LENGTH or len(three_prime) != THREE_LENGTH:

            write_dicts_to_csv([{"seq_id" : seq_id, "cluster" : cluster, "contig_length " : genome_dict[specie][strand_id]["len"], 
                                "FIVE_LENGTH" : FIVE_LENGTH, "THREE_LENGTH" : THREE_LENGTH,
                                "start" : coordinates[0],
                                "end" : coordinates[-1]}], 
                                f"output/{cov}/truncated_utr.csv", mode  = "a")


    if not is_max_Nter:    
        raw_UTRs["5utr"] = SeqRecord(seq = Seq(five_prime), id = f"{seq_id}-{cluster}", description = "")
        translated_UTRs["5utr"] = translate_frames(five_prime, specie = specie, 
                                            seq_id = seq_id, length = len(five_prime), 
                                            utr = "5utr", cluster = cluster)
        untranslated_elongates["Nter"] = None
        
    else:
        
        raw_UTRs["5utr"] = None
        translated_UTRs["5utr"] = None

        if five_prime != "":
            untranslated_elongates["Nter"] = SeqRecord(seq = Seq(five_prime), id = f"{seq_id}-{cluster}", description = "")
        else:
            untranslated_elongates["Nter"] = None

    if not is_max_Cter:
        raw_UTRs["3utr"] = SeqRecord(seq = Seq(three_prime), id = f"{seq_id}-{cluster}", description = "")
        translated_UTRs["3utr"] = translate_frames(three_prime, specie = specie, 
                                           seq_id = seq_id, length= len(three_prime), 
                                           utr = "3utr", cluster = cluster)
        untranslated_elongates["Cter"] = None
        
    else:

        raw_UTRs["3utr"] = None
        translated_UTRs["3utr"] = None
        if three_prime != "":
            untranslated_elongates["Cter"] = SeqRecord(seq = Seq(three_prime), id = f"{seq_id}-{cluster}", description = "")
        else:
            untranslated_elongates["Cter"] = None
    
        
    
    return translated_UTRs, raw_UTRs, untranslated_elongates



# Deprecated
def create_chimeric_sequences(chimeric_utr_dict, cds_infos, cds_dict):
    
    """
    Constructs chimeric sequences by combining the CDS with different versions of 5' UTR and 3' UTR sequences.

    Parameters:
    utr5_versions (dict of Bio.SeqRecord): Translated 5' UTR sequences for each frame.
    utr3_versions (dict of Bio.SeqRecord): Translated 3' UTR sequences for each frame.
    cds_dict (dict): Dictionary containing the CDS sequences for each specie.
    specie (str): Specie name.
    seq_id (str): Sequence id.
    length (int): Length of the elongation of 3' and 5' UTR.

    Returns:
    list of Bio.SeqRecord: A list containing nine chimeric sequences created by combining the CDS with each version of 5' UTR and 3' UTR.
    """

    chimeric_sequences = []
    specie = cds_infos[0]
    seq_id = cds_infos[1]
    length = cds_infos[2]/3
    cds = cds_dict[specie][seq_id]
    utr5_versions = chimeric_utr_dict["5utr"]
    utr3_versions = chimeric_utr_dict["3utr"]

    for i, (frame5,utr5) in enumerate(utr5_versions.items(), start = 1):
        for i, (frame3,utr3) in enumerate(utr3_versions.items(), start = 1):

            chimera = SeqRecord(seq = Seq(utr5.seq + cds + utr3.seq), 
                                id = f"{specie}_{seq_id}_frames:{frame5}_{frame3}_length:{length}", 
                                description = "")
            
            chimeric_sequences.append(chimera)

    return chimeric_sequences

