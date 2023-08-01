from Bio import SeqIO
import yaml
import gff3_parser
from utils.files import multifasta_to_dict
import polars as pl

with open('/home/simon.herman/Bureau/Gits/Elongates/env.yaml', 'r') as f:
    yaml_data = yaml.safe_load(f)
    species = yaml_data['Species_order']['Scer']

gff_dict = dict()
genome_dict = dict()

for specie in species:

    gff_dict[specie] = pl.from_pandas(gff3_parser.parse_gff3(f"input/{specie}.gff", parse_attributes = True, verbose = False))
    genome_dict[specie] = multifasta_to_dict(f"input/{specie}.fna", genome = True)

def is_good(protein_seq):

    bool_ = (not ('.' in protein_seq) or ('*' in protein_seq) or ('Z' in protein_seq)) & (protein_seq[0] == "M")
    return bool_


no_internal_stops = dict()
for specie in species:
    good_sequences = set()
    i = 0
    with open(f'input/{specie}_CDS.pep', 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if is_good(record.seq):

                good_sequences.add(record.id)

        no_internal_stops[specie] = good_sequences
    

bad = dict()
# Only those species are poorly annotated 
for specie in ["Skud","Smik","Sbay"]:

    no_stop_end = set()
    for sequence in gff_dict[specie].filter(pl.col("Type") == "CDS").iter_rows(named = True):

        cds_name = sequence["Name"]
        sense = sequence["Strand"]
        stop = int(sequence["End"])
        contig = sequence["Seqid"]

        if sense == "+":

            if genome_dict[specie][contig]["seq"][stop-3:stop] != "TGA" or genome_dict[specie][contig]["seq"][stop-3:stop] != "TAA" or genome_dict[specie][contig]["seq"][stop-3:stop] != "TAG":

                no_stop_end.add(cds_name)

        elif sense == "-":

            if genome_dict[specie][contig]["seq"][stop:stop+3][::-1] != "TGA" or genome_dict[specie][contig]["seq"][stop:stop+3][::-1] != "TAA" or genome_dict[specie][contig]["seq"][stop:stop+3][::-1] != "TAG":
                no_stop_end.add(cds_name)


    bad[specie] = good_sequences


for specie in ["Skud","Smik","Sbay"]:


    no_internal_stops[specie] == no_internal_stops[specie] - bad[specie]


    with open(f'input/{specie}_CDS.pep', 'r') as fasta_file:

        _ = set(no_internal_stops[specie])
        SeqIO.write([record for record in SeqIO.parse(fasta_file, 'fasta') if record.id in _ ], f"input/{specie}_CDS_corr.pep", "fasta")

