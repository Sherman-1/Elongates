import polars as pl 
import re
import yaml 
from utils.process import get_specie

# Load yaml file
with open('/home/simon.herman/Bureau/Gits/Elongates/env.yaml', 'r') as f:
    yaml_data = yaml.safe_load(f)
    species_dict = {k: v for d in yaml_data['Regex'].values() for k, v in d.items()}
    species_order = yaml_data['Species_order']['Scer']
    re_dict = yaml_data['Regex']['Scer']


def is_monophyletic(binary_vector):
    try:
        # find the first occurrence of '2'
        two_index = binary_vector.index(2)

        # find the first and last occurrence of '1' after '2'
        ones_after_two = binary_vector[two_index:]
        start_index = ones_after_two.index(1) + two_index
        end_index = len(ones_after_two) - 1 - ones_after_two[::-1].index(1) + two_index

        # check for '0' or '1' before '2' and '0' in between 'start_index' and 'end_index'
        if 1 in binary_vector[:two_index] or 0 in binary_vector[start_index:end_index+1]:
            return False
        else:
            return True

    except ValueError:
        return None


def parse_blast_dataframe(row : tuple) -> tuple:

    query_pattern = r'(.*?)-(cluster_n\d+)'
    subject_pattern = r'(.*?)-(cluster_n\d+)-.*-(f\d+)-.*'

    query_id = row[0]
    subject_id = row[1]
    evalue = row[2]
    qstart = row[3]
    qend = row[4]
    sstart = row[5]
    send = row[6]
    qseq = row[7]
    sseq = row[8]
    length = row[9]



    query_matches = re.match(query_pattern, query_id)
    query_seq_id = query_matches.group(1) if query_matches else None
    query_cluster = query_matches.group(2) if query_matches else None

    subject_matches = re.match(subject_pattern, subject_id)
    subject_seq_id = subject_matches.group(1) if subject_matches else None
    subject_cluster = subject_matches.group(2) if subject_matches else None
    subject_relative_frame = subject_matches.group(3) if subject_matches else None

    is_same_cluster = int(query_cluster == subject_cluster)
    q_specie = get_specie(re_dict, query_seq_id)
    s_specie = get_specie(re_dict, subject_seq_id)

    tuple_ = tuple([query_seq_id, subject_seq_id, evalue, qstart, qend, sstart, send, 
                    qseq, sseq, length, query_cluster, subject_cluster, subject_relative_frame, 
                    is_same_cluster, q_specie, s_specie])
    return tuple_


columns =  ["qseqid", "sseqid", "evalue", "qstart", "qend", "sstart", "send", "qseq", "sseq", "length"]
df_nter  = pl.read_csv("output/nter_five.tsv", separator="\t", has_header = False)
df_cter  = pl.read_csv("output/cter_three.tsv", separator="\t", has_header = False)
df_nter.columns = columns
df_cter.columns = columns

df_nter = df_nter.with_columns(
    pl.lit("NA").alias('query_cluster'),
    pl.lit("NA").alias('subject_cluster'),
    pl.lit("NA").alias('relative_frame'),
    pl.lit("NA").alias('same_cluster'),
    pl.lit("NA").alias('q_specie'),
    pl.lit("NA").alias('s_specie')
)

df_cter = df_cter.with_columns(
    pl.lit("NA").alias('query_cluster'),
    pl.lit("NA").alias('subject_cluster'),
    pl.lit("NA").alias('relative_frame'),
    pl.lit("NA").alias('same_cluster'),
    pl.lit("NA").alias('q_specie'),
    pl.lit("NA").alias('s_specie')
)

columns = df_nter.columns
df_nter = df_nter.apply(parse_blast_dataframe)
df_nter.columns = columns

columns = df_cter.columns
df_cter = df_cter.apply(parse_blast_dataframe)
df_cter.columns = columns              

nter_data = pl.read_csv("output/0.5_elongates.csv").select("cluster_size","seq_id",
                                                           "Nter_gaps","Nter_gap_openings","Nter_nb_aa",
                                                           "Nter_elongate_length","Nter_ratio",
                                                           "is_max_Nter","is_min_Nter",
                                                           "Nter_event_ID","Nter_events",
                                                           "Meth_after_Nter")

cter_data = pl.read_csv("output/0.5_elongates.csv").select("cluster_size","seq_id",
                                                            "Cter_gaps","Cter_gap_openings","Cter_nb_aa",
                                                            "Cter_elongate_length","Cter_ratio",
                                                            "is_max_Cter","is_min_Cter",
                                                            "Cter_event_ID","Cter_events")

df_nter = df_nter.join(nter_data, left_on="qseqid", right_on="seq_id").unique()
df_cter = df_cter.join(cter_data, left_on="qseqid", right_on="seq_id").unique()



tmp = df_nter.filter(

    (pl.col("same_cluster") == 1) &
    (pl.col("Nter_gap_openings") == 0) & 
    (pl.col("Nter_elongate_length") < 50) &
    (pl.col("evalue") < 1e-2)

).sort("qseqid")


tmp_id = list()
tmp_vector = list()
tmp_bool = list()
for query, matches in tmp.groupby("qseqid"):

    q_specie = matches["q_specie"].unique().to_list()[0]
    species = matches["s_specie"].unique().to_list()
    species = [1 if item in species else 2 if item == str(q_specie) else 0 for item in species_order]
    tmp_id.append(query)
    tmp_vector.append(str(species))
    tmp_bool.append(is_monophyletic(species))


test = tmp.join(pl.DataFrame({
    'query': tmp_id,
    'vector': tmp_vector,
    'is_monophyletic': tmp_bool
}), left_on="qseqid", right_on="query")
