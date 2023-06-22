"""
process_files.py

Author: HERMAN Simon
Date: 18/04/23
Version: 0.01

"""


__all__ = ["read_multifasta", "multifasta_to_dict", "extract_ids", "str_to_raw", "remove_temp_file", "write_dicts_to_csv"]

import os
import csv
from Bio.SeqIO import parse
from Bio.Seq import Seq
from .classes import _Seq


def remove_temp_file(cov, size, file = "temp.fasta",skip_check=False):

    """
    Remove the temporary FASTA file for a given coverage and cluster size.

    This function checks if the specified temporary FASTA file exists, and if it does,
    attempts to remove it. If the file is not found during the removal process, an
    appropriate error message is printed. The existence check can be skipped if needed.

    Parameters:
    cov (str): The coverage value used as a part of the file path.
    size (int or str): The cluster size used as a part of the file path.
    skip_check (bool, optional): Whether to skip checking for file existence. Defaults to False.
    """


    temp_file_path = f"work/{cov}/clusters/{size}/{file}"
    
    if not skip_check:
        if not os.path.exists(temp_file_path):
            return

    try:
        os.remove(temp_file_path)
    except FileNotFoundError:
        print(f"File error : {cov}, {size}, {temp_file_path}")

def read_multifasta(filename) -> list:


    """
    Reads a multi-FASTA file and returns a list of sequences.
    
    Args:
    filename (str): The path to the multi-FASTA file to be read.

    Returns:
    list: A list of sequence objects, where each object contains a sequence identifier and the corresponding sequence.
    """

    seq_list = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        seq_id = None
        sequence = ''
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    seq_list.append(_Seq(seq_id,sequence))
                    sequence = ""
                seq_id = line[1:]

            else:

                sequence = sequence + line

        if seq_id is not None:

            seq_list.append(_Seq(seq_id, sequence))

    return seq_list

def multifasta_to_dict(path, genome = False):


    """
    Reads a FASTA file and returns a dictionary where each record in the file
    is a key-value pair with the record identifier as the key and the record sequence (in uppercase) as the value.

    Args:
    path (str): The path to the FASTA file to be read.

    Returns:
    dict: A dictionary where the keys are the record identifiers and the values
          are the corresponding record sequences as Seq objects in uppercase letters.
    """
    
    records = parse(path, "fasta")

    if genome:
        dico = {}
        for record in records:
            dico[record.id] = {}
            dico[record.id]["seq"] = Seq(str(record.seq).upper())
            dico[record.id]["len"] = len(record.seq)
            
        return dico
    
    else:

        return {record.id: Seq(str(record.seq).upper()) for record in records}

def extract_ids(path):

    """
    Reads FASTA file and return a list of existing ids 

    Args:
    path (str) : Path to fasta file

    Returns:
    list: List of existing ids inside the fasta provided
    """

    ids = []
    with open(path,"r") as handle:

        for line in handle.readlines():
            if line.startswith("#"):
                continue
            elif line.startswith(">"):
                ids.append(line.strip().replace(">",""))

    return ids


def str_to_raw(s):

    """
    Converts a string to a raw string, escaping special characters.
    Taken from https://stackoverflow.com/questions/21605526/how-to-create-raw-string-from-string-variable-in-python 
    https://stackoverflow.com/users/1099876/ndpu
    """
    raw_map = {8:r'\b', 7:r'\a', 12:r'\f', 10:r'\n', 13:r'\r', 9:r'\t', 11:r'\v'}
    return r''.join(i if ord(i) > 32 else raw_map.get(ord(i), i) for i in s)


def extend_list_with_dict_values(input_dict, input_list):
    """
    Extend the input_list with all values from the input_dict.

    This function creates a new list that combines the input_list and all values from the input_dict.
    The input_list remains unchanged.

    Args:
        input_dict (dict): A dictionary containing lists as values.
        input_list (list): A list to be extended with the values from the input_dict.

    Returns:
        list: A new list containing the elements of the input_list and all the values from the input_dict.
    """
    result_list = input_list.copy()
    for value in input_dict.values():
        result_list.extend(value)
    return result_list

def write_dicts_to_csv(dicts, filename):
    """
    Appends a list of dictionaries to a CSV file, creating the file and writing a header row if necessary.

    Arguments:
    - dicts: List of dictionaries. All dictionaries must have the same keys.
    - filename: Output CSV filename.
    """
    # Check if file exists and is non-empty
    file_exists = os.path.isfile(filename) and os.path.getsize(filename) > 0

    with open(filename, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=dicts[0].keys())

        if not file_exists:
            # File didn't exist or was empty, so write the header
            writer.writeheader()

        for d in dicts:
            writer.writerow(d)
