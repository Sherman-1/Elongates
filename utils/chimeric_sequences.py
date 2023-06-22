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
from collections import defaultdict


def translate_frames(dna_sequence,specie, seq_id,length,utr):

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

    remainder = len(dna_seq) % 3
    if remainder != 0:
        dna_seq += "N" * (3 - remainder)

    # Translate DNA sequence into protein for 3 different reading frames
    # We divide length by 3 because the length of the UTR is given in nucleotides

    result["frame_0"] = SeqRecord(seq = Seq(dna_seq.translate()), 
                                  id = f"{seq_id}_{specie}_{utr}_frame:0_length:{length/3}", 
                                  description = specie)
    
    result["frame_1"] = SeqRecord(seq = Seq(dna_seq[1:].translate()), 
                                  id = f"{seq_id}_{specie}_{utr}_frame:1_length:{length/3}", 
                                  description = specie)
    
    result["frame_2"] = SeqRecord(seq = Seq(dna_seq[2:].translate()), 
                                  id = f"{seq_id}_{specie}_{utr}_frame:2_length:{length/3}", 
                                  description = specie)

    return result


def get_infos_for_chimeres(clusters, filtered_data):

    """
    Retrieves information for chimeric sequences from the filtered dataset. Input clusters correspond 
    to clusters where max Nter or Cter elongation is greater than fixed threshold defined at the beginning
    of the pipeline. 

    Parameters:
    --------
    clusters (set): The clusters ids to filter the dataset by. 

    Returns:
    --------
    out (dict): A dict of lists of tuples, where each tuple contains the specie, sequence id and length of the CDS.

    Structure of out:
    --------
    {cluster_id: [(specie, sequence_id, length), (specie, sequence_id, length), ...], ...}
    """

    out = defaultdict(list)
    
    for cluster in clusters:

        infos = filtered_data.filter(pl.col("num_cluster") == cluster)[["specie","sequence_id","max_length_Nter","max_length_Cter"]]

        length = max(
            int(infos[1]["max_length_Nter"].item()),
            int(infos[1]["max_length_Cter"].item())
        )*3 # 3 nucleotides per amino acid

        for row in infos.iter_rows(named=True):
            out[cluster].append((row["specie"], row["sequence_id"], length))

    return out

def get_extended_UTRs(cds_infos, gff_dict, genome_dict):

    """
    Compute the translations of UTR regions for a given CDS. Each UTR region ( 5' or 3' ) is translated
    into protein for the three different reading frames. The translated sequences are returned as a dictionary

    Parameters:
    --------
    cds_infos (tuple): Tuple containing the specie, sequence id and length of the CDS for which the UTRs are to be
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
    specie = cds_infos[0]
    seq_id = cds_infos[1]
    length = cds_infos[2]
    coordinates = [] # End and start coordinates of each CDS features
    result_dict = {} # Dictionary to store the results

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

        
        coordinates.append((int(row['Start'])-1, int(row['End'])-1))
        
    coordinates = sorted(coordinates, key=lambda x: x[0]) # Sort coordinates by start position

    # coordinates[0] = (start, end) of the first exon
    # coordinates[-1] = (start, end) of the last exon
    
    start_5 = coordinates[0][0]-length*2 if coordinates[0][0]-length*2 >= 0 else 0 # Get the start position of the 5' UTR
    end_3 = coordinates[-1][1]+1+length*2 if coordinates[-1][1]+1+length*2 <= genome_dict[specie][strand_id]["len"] else genome_dict[specie][strand_id]["len"] # Get the end position of the 3' UTR
    
    five_prime = genome_dict[specie][strand_id]["seq"][

        start_5:coordinates[0][0]

        ] # Get the 5' sequence
    
        # +1 for -1,1 because GFF points to the last nucleotide of the stop codon

    three_prime = genome_dict[specie][strand_id]["seq"][
        coordinates[-1][1]+1:end_3
        ] # Get the 3' sequence


    # Reverse complement if the strand is negative, don't forget to reverse the coordinates
    if strand == "-":

        five_prime, three_prime = three_prime.reverse_complement(), five_prime.reverse_complement()
    
    # Translate the sequences for each frame
    result_dict["5utr"] = translate_frames(five_prime, specie = specie, seq_id = seq_id, length = length*2, utr = "5utr")
    result_dict["3utr"] = translate_frames(three_prime, specie = specie, seq_id = seq_id, length= length*2, utr = "3utr")

    return result_dict

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

