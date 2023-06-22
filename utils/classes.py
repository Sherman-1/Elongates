

__all__ = ["_Seq", "Cluster"]

class _Seq:

    """
    Own temporary Seq class, BioPython can't handle MUSCLE fasta format output for some reason
    """
    
    def __init__(self, seq_id, sequence):
        self.id = seq_id
        self.seq = sequence


class Cluster:

    def __init__(self, cluster_id, size):
        self.id = cluster_id
        self.size = size
        self.seq_list = None

    def add_seq(self, seq):
        self.seq_list.append(seq)