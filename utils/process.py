
import warnings
import re 
import numpy as np
import itertools

def get_specie(re_dict, id):

    """
    Function reads dictionnary of key:value specie:RNA name re.pattern. Associate id to one of 
    the species found in the dictionnary. 
    If more than one specie is found, a warning is raised. If no specie is found, returns "unknown".
    For know, dict is defined in env.yaml file.

    Parameters:
    ----------
    re_dict (dict): A dictionary containing species as keys and their corresponding regular expression patterns as values.
    id (str): An identifier that represents the specie in the input data.

    Returns:
    ----------
    specie (str): The specie corresponding to the given identifier. If more than one specie is found, a warning is raised.
                If no specie is found, the function returns "unknown".
    
    Example:
    ----------

    rna name could be : RNA-XM-12345678.1
    rna pattern than will be RNA-XM-([0-9]+).1
    If XM corresponds to Saccharomyces cerevisiae, the dictionnary will be :
    re_dict = {"Saccharomyces cerevisiae" : "RNA-XM-([0-9]+).1"}
    """

    i = 0
    for key in re_dict.keys():

        if re_dict[key] in id:
            specie = key
            i += 1

    if i > 1:
        warnings.warn("More than one specie found for %s" % id)
    
    elif i == 0:
        specie = "unknown"
        
    return specie


def is_methionine_after_nter(sequence: str, nter_length: int) -> int:

    """
    This function checks if the first amino acid following an N-terminal elongation is a Methionine.
    If the Methionine is found just after the elongation, it suggests that the elongation could be due 
    to a misplaced start codon in the genome annotation (GFF file), rather than a biological event.
    
    Parameters:
    sequence (str): The amino acid sequence of the protein.
    nter_length (int): The length of the N-terminal elongation.
    
    Returns:
    bool: True if the first amino acid following the elongation is Methionine, False otherwise.
    """
    try:
        next_aa_after_elongation = sequence[nter_length]
    except IndexError:
        return False  # nter_length is longer than the sequence itself
    return int(next_aa_after_elongation.upper() == 'M')

def count_dashes(string,reverse = False):

    """
    Count the number of dashes "-" on N-ter and C-ter ( 3' 5' ) 

    Args:
    --------
        string (str): The input string 
        reverse (bool, optional): If True, count dashes from the end of the string. Defaults to False.

    Returns:
    --------
        int: The number of dashes at the N-ter or C-ter of the input string.

    Why do we need this function ?:
    --------

        In MUSCLE alignment format, dashes "-" are used to represent gaps in the alignment.
        In our case, we want to count the number of dashes at the N-ter or C-ter of the sequences to 
        determine the length of the N-ter or C-ter elongation for the other sequences. 
        
        YYYYYYYYY-----
        YYYYYYYYYYYYYY

        In this example, our goal is to retrieve the Cter elongation of the second sequence. To compute it, we first
        need to count the number of dashes at the N-ter of the OTHER sequence. 

    """

    count = 0
    i = 0
    if reverse == True : string = string[::-1]

    while i < len(string) and string[i] == "-":
        count += 1
        i += 1
    return count

import networkx as nx

def subgraph_counts(lengths_list,cluster_id, cluster_size, localisation, threshold = 5):

    """
    Create a complete graph with the nodes as the numbers in the input list. Edges are drawn between nodes only if 
    the absolute difference between the two nodes is less than threshold value.
    In this case, the nodes are the lengths of every Nter OR Cter elongates. The weight of an edge is the absolute
    difference between the lengths of the two nodes. 
    The function removes edges with weight greater than threshold and returns the number of subgraphs in the remaining graph.

    Because 15 amino acids is the minimum length for a peptide to form an helix, it's set as the default threshold value.
    Consideration only valable for our study scope.

    Args:
        num_list (list): A list of numbers to be used as nodes in the graph.
        threshold (int): The maximum allowed weight for an edge in the graph.

    Returns:
        int: The number of subgraphs in the remaining graph after removing edges that have a weight greater than the
        threshold value.

    Raises:
        NetworkXError: If the graph contains self-loops.

    Example:
        >>> lst = [1,4,40,41,60]
        >>> subgraph_counts(lst, threshold = 10)
        The number of subgraphs is: 3
    """
    G = nx.Graph()
    toplist = list()
    cluster_size = re.match(pattern="size_(\d+)", string=cluster_size).group(1)


    # Add nodes to the graph with custom names and values
    for name, value in enumerate(lengths_list):
        G.add_node(name, value=value)


    # Add edges between nodes based on the absolute difference of their values
    threshold = 10
    for node1, node2 in itertools.combinations(G.nodes, 2):

        value1 = G.nodes[node1]['value']
        value2 = G.nodes[node2]['value']
        
        if abs(value1 - value2) < threshold:
            G.add_edge(node1, node2)

    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    nb_subgraphs = nx.number_connected_components(G)-1
    subgraph_id = 1
    for sub in S:

        lengths = [data["value"] for name,data in sub.nodes(data=True)]
        subdict = {
            'cluster_id': cluster_id,
            'cluster_size': cluster_size,
            'subgraph_id': f"subgraph_{subgraph_id}",
            'localisation': localisation,
            'nb_subgraphs': nb_subgraphs,
            'sequences_in_event': len(lengths),
            'std_dev_lengths': np.std(lengths),
            'mean_lengths': np.mean(lengths),
        }
        toplist.append(subdict)
        subgraph_id += 1

    return toplist, nb_subgraphs

def get_elongates(sequence, max_length, upstream = False):

    """
    Extract a subsequence from a given sequence starting from extremities up to a specified maximum length, considering
    either upstream or downstream direction. The function also returns the number of gaps, gap openings,
    the number of amino acids and the total length of the subsequence.

    Specified maximum length and length of subsequence are not necessarily the same since the function starts to
    extract the subsequence only after encountering the first non-gap character.

    Args:
    --------
        sequence (str): The input sequence as a string, where "-" represents a gap.
        max_length (int): The maximum length of the subsequence to be extracted.
        upstream (bool, optional): If True, extract the subsequence from the start of the sequence (upstream),
            otherwise extract from the end of the sequence (downstream). Defaults to False.

    Returns:
    --------
        dict: A dict containing the following elements:
            - elongated_subsequence (str): The extracted elongated subsequence as a string.
            - gaps (int): The total number of gaps in the elongated subsequence.
            - gap_openings (int): The number of gap openings in the elongated subsequence.
            - nb_aa (int): The number of amino acids in the elongated subsequence.
            - len(elongate_seq) (int, optional): The length of the elongated subsequence. Only
                returned if upstream is True.
    """

    i = 0
    elongate_seq = []
    gaps = 0
    gap_openings = 0
    nb_aa = 0

    if upstream == False:
        while sequence[len(sequence)-i-1] == "-" and i < max_length:
            i+=1

        while i < max_length:

            if sequence[len(sequence)-i-1] == "-":
                gaps += 1
                if sequence[len(sequence)-i] != "-":
                    gap_openings += 1
            else: 
                nb_aa += 1

            elongate_seq.append(sequence[len(sequence)-i-1])
            i += 1

        elongate_seq = ''.join(elongate_seq)[::-1]
        return {
            'Cter_elongate': elongate_seq,
            'Cter_gaps': gaps,
            'Cter_gap_openings': gap_openings,
            'Cter_nb_aa': nb_aa,
            'Cter_elongate_length': len(elongate_seq.replace("-", ""))
                }
    
    elif upstream == True:
        
        while sequence[i] == "-" and i < max_length:
            i +=1

        while i < max_length:

            if sequence[i] == "-":
                gaps += 1
                if sequence[i-1] != "-":
                    gap_openings += 1

            else:
                nb_aa += 1

            elongate_seq.append(sequence[i])
            i += 1

        elongate_seq = ''.join(elongate_seq)   
        return {
            'Nter_elongate': elongate_seq,
            'Nter_gaps': gaps,
            'Nter_gap_openings': gap_openings,
            'Nter_nb_aa': nb_aa,
            'Nter_elongate_length': len(elongate_seq.replace("-", ""))
                }

def compute_ratios(elongate_length, seq_length):

    """
    Compute the Nter/Cter elongation length to sequence length.

    Parameters:
    -----------
    elongate : int
        Length of the Nter/Cter elongation.
    seq_length : int
        Length of the sequence.

    Returns:
    --------
    float
        Ratio of dashes to sequence length. Returns 0 if the number of dashes is 0.
    """
    return elongate_length / seq_length if elongate_length > 0 else 0


def process_record(record, cluster_name, cluster_size, regex_dict):

    """
    Process a single sequence record and return the corresponding sublist.

    Parameters:
    -----------
    record : SeqRecord
        A single sequence record from a MUSCLE MSA.
    cluster_name : str
        The name of the cluster.
    cluster_size : str
        The size of the cluster, formatted as "size_X" where X is the size.
    regex_dict : dict
        A dictionary containing regex patterns for extracting species information from the sequence IDs.

    Returns:
    --------
    record_info : dict
        A dict containing information for a single sequence in the cluster.
    
    """

    id = record.id
    tmp_seq = str(record.seq).replace("\n", "")
    Nter_dashes = count_dashes(tmp_seq)
    Cter_dashes = count_dashes(tmp_seq, reverse=True)
    
    cluster_size = re.match(pattern="size_(\d+)", string=cluster_size).group(1)

    record_info = {

        'cluster_size': cluster_size,
        'cluster_name': str(cluster_name),
        'seq_id': id,
        'species': get_specie(regex_dict, id),
        'sequence': tmp_seq,
        'sequence_length': len(tmp_seq),
        'Nter_dashes': Nter_dashes,
        'Cter_dashes': Cter_dashes

    }


    # Nter_dashes and Cter_dashes are returned separately to avoid computing them twice 
    return record_info

def update_record_info(record_infos, max_Nter, max_Cter):

    """
    Update the information of a record with the elongation lengths.

    Parameters:
    record_infos (dict): dictionary containing the record's information.
    max_Nter (int): maximum elongation length at the N-terminus.
    max_Cter (int): maximum elongation length at the C-terminus.

    Returns:
    record_infos (dict): Updated dictionary of the record's information.
    """

    Nter_infos = get_elongates(record_infos["sequence"], max_length=max_Nter, upstream=True)
    Cter_infos = get_elongates(record_infos["sequence"], max_length=max_Cter)

    record_infos["max_Nter"] = max_Nter
    record_infos["max_Cter"] = max_Cter
    record_infos.update(Nter_infos)
    record_infos.update(Cter_infos)
    record_infos["Nter_ratio"] = compute_ratios(record_infos["Nter_elongate_length"], record_infos["sequence_length"])
    record_infos["Cter_ratio"] = compute_ratios(record_infos["Cter_elongate_length"], record_infos["sequence_length"])
    record_infos["is_max_Nter"] = 1 if record_infos["Nter_elongate_length"] == max_Nter else 0
    record_infos["is_max_Cter"] = 1 if record_infos["Cter_elongate_length"] == max_Cter else 0
    record_infos["is_min_Nter"] = 1 if record_infos["Nter_elongate_length"] == 0 else 0
    record_infos["is_min_Cter"] = 1 if record_infos["Cter_elongate_length"] == 0 else 0
    return record_infos

def compute_max_Nter_Cter(records, cluster_name, cluster_size, regex_dict):

    """
    Compute the maximum elongation length for a set of records.

    Parameters:
    records (list): list of records.
    cluster_name (str): name of the cluster.
    cluster_size (int): size of the cluster.
    regex_dict (dict): dictionary of regex patterns.

    Returns:
    max_Nter (int): maximum elongation length at the N-terminus.
    max_Cter (int): maximum elongation length at the C-terminus.
    infos_list (list): list of dictionaries containing record information.
    """

    max_Nter = 0
    max_Cter = 0
    infos_list = []

    for record in records:
        record_infos = process_record(record, cluster_name, cluster_size, regex_dict)
        infos_list.append(record_infos)
        max_Nter = max(max_Nter, record_infos['Nter_dashes'])
        max_Cter = max(max_Cter, record_infos['Cter_dashes'])
        
    return max_Nter, max_Cter, infos_list

def get_elongation_events(Nter_lengths, Cter_lengths, cluster_name, cluster_size):

    """
    Get the elongation events for a set of records.

    Parameters:
    Nter_lengths (list): list of elongation lengths at the N-terminus.
    Cter_lengths (list): list of elongation lengths at the C-terminus.
    cluster_name (str): name of the cluster.
    cluster_size (int): size of the cluster.

    Returns:
    Nter_events_infos (list): list of dictionaries containing N-terminus elongation event information.
    Cter_events_infos (list): list of dictionaries containing C-terminus elongation event information.
    """

    Nter_events_infos, Number_nter_events = subgraph_counts(
        Nter_lengths,
        cluster_id = cluster_name, 
        cluster_size = cluster_size, 
        localisation = "Nter", 
        threshold = 5)
    
    Cter_events_infos, Number_cter_events = subgraph_counts(
        Cter_lengths,
        cluster_id = cluster_name, 
        cluster_size = cluster_size, 
        localisation = "Cter", 
        threshold = 5)

    return Nter_events_infos, Cter_events_infos, Number_nter_events, Number_cter_events


def process_multiple_records(records, cluster_name, cluster_size, regex_dict):

    """
    Process multiple records by computing max elongations, updating records with elongations, and gathering events.

    Parameters:
    records (list): list of records.
    cluster_name (str): name of the cluster.
    cluster_size (int): size of the cluster.
    regex_dict (dict): dictionary of regex patterns.

    Returns:
    infos_list (list): list of dictionaries containing updated record information.
    event_lists (list): list of dictionaries containing elongation event information.
    """
    max_Nter, max_Cter, elongates_infos_list = compute_max_Nter_Cter(records, cluster_name, cluster_size, regex_dict)
    Nter_lengths = []
    Cter_lengths = []
    for record_infos in elongates_infos_list:
        record_infos = update_record_info(record_infos, max_Nter, max_Cter)
        Nter_lengths.append(record_infos["Nter_elongate_length"])
        Cter_lengths.append(record_infos["Cter_elongate_length"])

    Nter_events_infos, Cter_events_infos, Number_nter_events, Number_cter_events = get_elongation_events(Nter_lengths, Cter_lengths, cluster_name, cluster_size)
    
    for record_infos in elongates_infos_list:
        record_infos["Nter_events"] = Number_nter_events
        record_infos["Cter_events"] = Number_cter_events
        record_infos["Meth_after_Nter"] = is_methionine_after_nter(sequence = record_infos["sequence"],
                                                                   nter_length = record_infos["Nter_elongate_length"])

    event_lists = []
    if Nter_events_infos is not None:
        event_lists.extend(Nter_events_infos)
    if Cter_events_infos is not None:
        event_lists.extend(Cter_events_infos)
    
    return elongates_infos_list, event_lists
