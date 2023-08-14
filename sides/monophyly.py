import re
from Bio import SeqIO
import yaml
import csv
import os

base_dir = '/home/sherman/Bureau/Gits/Elongates/work/0.5/clusters'

with open('/home/sherman/Bureau/Gits/Elongates/env.yaml', 'r') as f:
    yaml_data = yaml.safe_load(f)
    species_dict = {k: v for d in yaml_data['Regex'].values() for k, v in d.items()}
    species_order = yaml_data['Species_order']['Scer']

# Function to check which species are present in a multifasta file
def species_in_cluster(fasta_file, species_dict, species_order):
    species_in_file = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        for species in species_order:
            regex = species_dict.get(species, None)
            if regex and re.search(regex, record.description):
                species_in_file.append(species)
    return species_in_file

def is_monophyletic(binary_vector):
    try:
        # find the first and last occurrence of '1'
        start_index = binary_vector.index(1)
        end_index = len(binary_vector) - 1 - binary_vector[::-1].index(1)

        # check for '0' in between 'start_index' and 'end_index'
        if 0 in binary_vector[start_index:end_index+1]:
            return False
        else:
            return True
        
    except ValueError:
        return None

# Open output file
with open('output/species.csv', 'w') as f:
    writer = csv.writer(f)

    # Write header
    writer.writerow(['Cluster'] + species_order + ['Monophyletic'])

    # Go through each cluster size directory
    for subdir in os.listdir(base_dir):
        subdir_path = os.path.join(base_dir, subdir)
        if os.path.isdir(subdir_path):
            # Go through each multifasta file in the directory
            for file in os.listdir(subdir_path):
                file_path = os.path.join(subdir_path, file)
                if os.path.isfile(file_path) and not file_path.startswith('temp.fasta'):  # adjust the condition if needed
                    # Check which species are present in the multifasta file
                    species_in_file = species_in_cluster(file_path, species_dict, species_order)
                    # Write binary vector to the output file
                    binary_vector = [1 if species in species_in_file else 0 for species in species_order]
                    mono = is_monophyletic(binary_vector)
                    writer.writerow([file] + binary_vector + [mono])