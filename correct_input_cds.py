from Bio import SeqIO
import yaml

with open('/home/simon.herman/Bureau/Gits/Elongates/env.yaml', 'r') as f:
    yaml_data = yaml.safe_load(f)
    species = yaml_data['Species_order']['Scer']

def is_good(protein_seq):

    boolo = not ('.' in protein_seq) or ('*' in protein_seq) or ('Z' in protein_seq)
    return boolo

with open("bad_sequences.txt", "w") as f:

    
    for specie in species:
        good_sequences = list()
        i = 0
        with open(f'input/{specie}_CDS.pep', 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if is_good(record.seq):

                    good_sequences.append(record)
      
        SeqIO.write(good_sequences, f'input/{specie}_CDS_corr.pep', 'fasta')