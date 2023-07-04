import re
import polars as pl
from Bio.Blast import NCBIXML

# Define regex patterns for the query and subject
query_pattern = r'(.*?)-(cluster_n\d+)'
subject_pattern = r'(.*?)-(cluster_n\d+)-.*-(f\d+)-.*'

# Parse the BLAST XML output
blast_records = NCBIXML.parse(open("output/nter_five.xml","r"))

# Create an empty list to hold the parsed data
parsed_data = []

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            # Extract the seq_id, cluster from the query
            query_matches = re.match(query_pattern, blast_record.query)
            query_seq_id = query_matches.group(1) if query_matches else None
            query_cluster = query_matches.group(2) if query_matches else None
            
            # Extract the seq_id, cluster, and relative frame from the subject
            subject_matches = re.match(subject_pattern, alignment.hit_def)
            subject_seq_id = subject_matches.group(1) if subject_matches else None
            subject_cluster = subject_matches.group(2) if subject_matches else None
            subject_relative_frame = subject_matches.group(3) if subject_matches else None
            
            # Check if the query and the subject are from the same cluster
            same_cluster = int(query_cluster == subject_cluster)
            
            # Append the parsed data to the list
            parsed_data.append([query_seq_id, query_cluster, subject_seq_id, subject_relative_frame, same_cluster, hsp.expect])

# Convert the parsed data into a polars DataFrame
df = pl.DataFrame(parsed_data, schema=['Query', 'Cluster', 'Subject', 'Frame', 'Intra', 'E-value'])

# Sort the data by E-value in ascending order
df = df.sort(['E-value'], descending=False)

# Group the data by 'Query Seq ID' and 'Same Cluster', and select the first row of each group
df = df.groupby(['Query', 'Intra']).apply(lambda df: df.head(1))

df = df.sort(["Cluster"])


elongates = pl.read_csv("output/0.5_elongates.csv", has_header = True)

columnsToKeep = ['Query', 'Cluster', 'cluster_size', 'Subject', 'Frame', 'Intra', 'E-value', 'Nter_elongate', 'Nter_gaps','Nter_gap_openings', 'Nter_ratio', 'Meth_after_Nter' ]
df = df.join(elongates, left_on="Query", right_on="seq_id")[columnsToKeep]

df.write_csv("output/nter_five.csv")
