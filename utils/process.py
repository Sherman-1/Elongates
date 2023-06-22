
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

def subgraph_counts(lengths_list,cluster_id, cluster_size, localisation, threshold = 15):

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
        sublist = [cluster_id,
                    cluster_size, 
                    f"subgraph_{subgraph_id}", # Id of the event
                    localisation, # Nter or Cter
                    nb_subgraphs, # Total Nb of elongation events in the cluster
                    len(lengths), # Nb of sequences in this particular elongation event 
                    np.std(lengths), # Standard deviation of the lengths of the elongates in this particular elongation event
                    np.mean(lengths), # Mean of the lengths of the elongates in this particular elongation event
        ]
        toplist.append(sublist)
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
        tuple: A tuple containing the following elements:
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

        return (''.join(elongate_seq)[::-1], gaps, gap_openings, nb_aa, len(elongate_seq))
    
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
            
        return (''.join(elongate_seq), gaps, gap_openings, nb_aa, len(elongate_seq))

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
    sublist : list
        A list containing information for a single sequence in the cluster.
    Nter_dashes : int
        Number of dashes at the N-terminus of the sequence.
    Cter_dashes : int
        Number of dashes at the C-terminus of the sequence.
    """

    id = record.id
    tmp_seq = str(record.seq).replace("\n", "")
    Nter_dashes = count_dashes(tmp_seq)
    Cter_dashes = count_dashes(tmp_seq, reverse=True)
    
    cluster_size = re.match(pattern="size_(\d+)", string=cluster_size).group(1)

    sublist = [
        cluster_size,
        str(cluster_name),
        id,
        get_specie(regex_dict, id),
        tmp_seq,
        len(tmp_seq),
        Nter_dashes,
        Cter_dashes
    ]

    # Nter_dashes and Cter_dashes are returned separately to avoid computing them twice 
    return sublist, Nter_dashes, Cter_dashes



def process_multiple_records(records, cluster_name, cluster_size, regex_dict):

    elongates_list = []
    infos_lists = []
    Nter_lengths = []
    Cter_lengths = []
    max_Nter = 0
    max_Cter = 0

    # Process each record using process_record()
    for record in records:

        # sublist = [cluster size, cluster_name, seq id, specie, CDS, CDS length, Nter_dashes, Cter_dashes]
        sublist, Nter_dashes, Cter_dashes = process_record(record, cluster_name, cluster_size, regex_dict)

        
        elongates_list.append(sublist)

        # For each sequence in the cluster, we keep track of the maximum N-terminus and C-terminus number 
        # of dashes to compute the elongated sequences later. 
        # For two sequences aligned by MUSCLE, dashes represent the gaps, the " lack " of amino acids 
        # compared to the other sequence. 
        max_Nter = max(max_Nter, Nter_dashes)
        max_Cter = max(max_Cter, Cter_dashes)
    
    
    # Once every record from records has been seen, we can fix the maximum N-terminus 
    # and C-terminus length difference between sequences in the cluster
    # We use these infos to compute the elongated sequences for every sequence in the cluster
    for list in elongates_list:

        # infos = [elongated sequence, gaps, gap_openings, nb_aa, length of elongated sequence]
        Nter_infos = get_elongates(list[4], max_length=max_Nter, upstream=True)
        Cter_infos = get_elongates(list[4], max_length=max_Cter)

        list.append(max_Nter)
        list.append(max_Cter)
        list.extend(Nter_infos)
        list.extend(Cter_infos)
        list.append(compute_ratios(Nter_infos[-1], list[5]))
        list.append(compute_ratios(Cter_infos[-1], list[5]))

        # Checks if the elongated sequence is the longest one for the N-terminus and C-terminus
        if Nter_infos[-1] == max_Nter:
            list.append(1)
        else:
            list.append(0)
        if Cter_infos[-1] == max_Cter:
            list.append(1)
        else:
            list.append(0)

        
        Nter_lengths.append(Nter_infos[-1])
        Cter_lengths.append(Cter_infos[-1])

    # Compute the number of elongation events per cluster for each side
    # See subgraph_counts() function for more details

    # Also compute infos about the elongation events
    # elongations_infos is a list of lists, each sublist containing the following information :
    # [ cluster id, cluster_size, event_id, nb of event in the cluster, nb of sequences in a particular elongation event, sd of elongate length for current event ]
    Nter_events_infos, Nter_nb_elongation_events = subgraph_counts(Nter_lengths, cluster_id = cluster_name, cluster_size = cluster_size, localisation = "Nter", threshold = 10)
    Cter_events_infos, Cter_nb_elongation_events = subgraph_counts(Cter_lengths, cluster_id = cluster_name, cluster_size = cluster_size, localisation = "Cter", threshold = 10)

    if Nter_events_infos != None:
        infos_lists.extend(Nter_events_infos)
    if Cter_events_infos != None:
        infos_lists.extend(Cter_events_infos)
    
    for list in elongates_list:

        list.extend((Nter_nb_elongation_events, Cter_nb_elongation_events)) # need to turn elongations_infos into an OD to avoid hard coding the index

    # At the end of this function, every sublist in infos_list contains the following information :

    # [size, cluster_name, seq_id, specie, sequence, length, Nter_dashes, Cter_dashes, 
    # max_Nter, max_Cter, Elongate length, gaps, gap_openings, nb_aa, length_elongate, Nter_subgraphs, Cter_subgraphs, 
    # Nter_ratio, Cter_ratio,] to be updated

    return elongates_list, infos_lists

