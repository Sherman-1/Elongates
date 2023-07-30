import polars as pl

cov = 0.5
seuil = 4

elongates = pl.read_csv(f"output/{cov}/{cov}_elongates.csv", infer_schema_length =10000)

########
# Get the clusters that are elongated in Scer but not in Sbay
# Knowing well that such clusters include clusters without Sbay, or only with Scer or any other scenario 
########

Nter_scer_conditions = ((pl.col("species") == "Scer_NCBI") & (abs(pl.col("max_Nter") - pl.col("Nter_nb_aa")) < seuil) & (pl.col("max_Nter") >= 10))
Nter_sbay_conditions = ((pl.col("species") == "Sbay") & (pl.col("Nter_nb_aa") > seuil))

Nter_scer_clusters = set(elongates.filter(Nter_scer_conditions)["cluster_name"].to_list())
Nter_sbay_clusters = set(elongates.filter(Nter_sbay_conditions)["cluster_name"].to_list())

Nter_clusters = elongates.filter(pl.col("cluster_name").is_in(Nter_scer_clusters- Nter_sbay_clusters))



Cter_scer_conditions = ((pl.col("species") == "Scer_NCBI") & (abs(pl.col("max_Cter") - pl.col("Cter_nb_aa")) < seuil) & (pl.col("max_Cter") >= 10))
Cter_sbay_conditions = ((pl.col("species") == "Sbay") & (pl.col("Cter_nb_aa") > seuil))

Cter_scer_clusters = set(elongates.filter(Cter_scer_conditions)["cluster_name"].to_list())
Cter_sbay_clusters = set(elongates.filter(Cter_sbay_conditions)["cluster_name"].to_list())

Cter_clusters = elongates.filter(pl.col("cluster_name").is_in(Cter_scer_clusters- Cter_sbay_clusters))


########
# From those clusters, we want to keep only the ones that have Scer AND Sbay
########


# Nter
scer_df_nter = set(Nter_clusters.filter(pl.col("species") == "Scer_NCBI")["cluster_name"].to_list())
sbay_df_nter = set(Nter_clusters.filter(pl.col("species") == "Sbay")["cluster_name"].to_list())
common = scer_df_nter.intersection(sbay_df_nter)

full_filtered_Nter = Nter_clusters.filter(pl.col("cluster_name").is_in(common)).sort("cluster_name")

# Cter
scer_df_cter = set(Cter_clusters.filter(pl.col("species") == "Scer_NCBI")["cluster_name"].to_list())
sbay_df_cter = set(Cter_clusters.filter(pl.col("species") == "Sbay")["cluster_name"].to_list())
common = scer_df_cter.intersection(sbay_df_cter)

full_filtered_Cter = Cter_clusters.filter(pl.col("cluster_name").is_in(common)).sort("cluster_name")


# Question : are there some clusters that correspond to both Nter and Cter elongation criteria ?

# Answer : 

# full_filtered_Nter.filter(pl.col("cluster_name").is_in(full_filtered_Cter["cluster_name"].unique().to_list()))


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils.handle_UTRs import translate_frames


def custom_target_elongate(cluster, scer_length, seq_id, specie, elongate_length, side, genome_dict, gff_dict): 


    coordinates = []

    if side != "Nter" and side != "Cter":

        raise ValueError("Side must be either Nter or Cter")
    
    gff = gff_dict[specie].filter(

        (pl.col("Type") == "CDS") & ((pl.col("Name") == seq_id) | (pl.col("Parent") == seq_id))

    )[["Start","End","Strand","Seqid"]] # Keep only necessary columns
    
    # Store datas necessary to compute the elongate sequence

    strand = gff[0]["Strand"].to_list()[0] # + or -
    strand_id = gff[0]["Seqid"].to_list()[0] # chromosome or scaffold id

    for row in gff.iter_rows(named=True): # Named = True to iter with column names

        coordinates.append(sorted((int(row['Start'])-1, int(row['End'])-1))) # -1 for python indexing
        
    coordinates = sorted(coordinates, key=lambda x: x[0]) # Sort coordinates by start position


    custom_elongate_length = scer_length - elongate_length 

    if strand == "+":

        if side == "Nter":

            start_5 = coordinates[0][0]-custom_elongate_length if coordinates[0][0]-custom_elongate_length >= 0 else 0 # Get the start position of the 5' UTR
            
            elongate = genome_dict[specie][strand_id]["seq"][

                start_5:coordinates[0][0]

            ] 

        elif side == "Cter": # Useless check but it's for the sake of clarity
        
            end_3 = coordinates[-1][1]+1+custom_elongate_length if coordinates[-1][1]+1+custom_elongate_length <= genome_dict[specie][strand_id]["len"] else genome_dict[specie][strand_id]["len"] # Get the end position of the 3' UTR
            # +1 for -1,1 because GFF points to the last nucleotide of the stop codon

            elongate = genome_dict[specie][strand_id]["seq"][
                coordinates[-1][1]+1:end_3
                ] # Get the 3' sequence

    
    # Reverse complement if the strand is negative, don't forget to reverse the coordinates
    if strand == "-":

        if side == "Nter":

            end_5 = coordinates[-1][1]+1+custom_elongate_length if coordinates[-1][1]+1+custom_elongate_length <= genome_dict[specie][strand_id]["len"] else genome_dict[specie][strand_id]["len"] # Get the start position of the 5' UTR

            elongate = genome_dict[specie][strand_id]["seq"][
                coordinates[-1][1]+1:end_5
            ].reverse_complement() # Get the 5' sequence
        

        if side == "Cter":
        
            start_3 = coordinates[0][0]-custom_elongate_length if coordinates[0][0]-custom_elongate_length >= 0 else 0 # Get the end position of the 3' UTR

            elongate = genome_dict[specie][strand_id]["seq"][
                start_3:coordinates[0][0]
            ].reverse_complement() # Get the 3' sequence


    nucleotic_seq = SeqRecord(seq = Seq(elongate), id = f"{seq_id}-{cluster}", description = "")

    # def translate_frames(dna_sequence, specie, seq_id, length, utr, cluster)

    frames_dict = translate_frames(dna_sequence = nucleotic_seq, specie = specie, seq_id = seq_id, length = custom_elongate_length)

    return 0


dataframes = {
    "Nter": full_filtered_Nter,
    "Cter": full_filtered_Cter,
}

for side in ["Nter", "Cter"]:

    df = dataframes[side]

    #os.mkdir(f"{current_path}/local_align_files/{side}")

    for cluster, sequences in df.groupby("cluster_name"):

        #os.mkdir(f"{current_path}/local_align_files/{side}/{cluster}")

        scer_length = sequences.filter(pl.col("species") == "Scer_NCBI")[f"{side}_nb_aa"].max() # Maybe several Scer sequences in the cluster, we take the longest elongate

        input_dict = dict()

        for sequence in sequences.iter_rows(named = True): 

            if scer_length - sequence[f"{side}_nb_aa"] >= 10:

                dict_ = custom_target_elongate(cluster, scer_length, sequence["seq_id"], sequence["species"], sequence[f"{side}_nb_aa"], side, genome_dict, gff_dict)

                input_dict[sequence["seq_id"]] = sequence["Nter_nb_aa"] 
                #os.mkdir(f"{current_path}/local_align_files/{side}/{cluster}/{sequence['seq_id']}")
                #os.mkdir(f"{current_path}/local_align_files/{side}/{cluster}/{sequence['seq_id']}/nucleotide")
                #os.mkdir(f"{current_path}/local_align_files/{side}/{cluster}/{sequence['seq_id']}/protein")