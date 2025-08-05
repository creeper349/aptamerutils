import zlib

def _compress_size(s:str):
    return len(zlib.compress(s.encode("utf-8")))

def _kmer_set(seq:str, kmer:int) -> set:
    return set([seq[i:i + kmer] for i in range(len(seq) - kmer + 1)])

def ncd(seq1:str, seq2:str, **kwargs):
    c1 = _compress_size(seq1)
    c2 = _compress_size(seq2)
    cm = _compress_size(seq1 + seq2)
    return (cm - min(c1, c2)) / max(c1, c2)

def editing_distance(seq1:str, seq2:str, **kwargs):
    score = 0
    if len(seq1) > len(seq2):
        m = seq1
        seq1 = seq2
        seq2 = m
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score += 1
    return (len(seq1) - score) / len(seq1)

def jaccard(seq1:str, seq2:str, kmer = 5, **kwargs):
    k1 = _kmer_set(seq1, kmer)
    k2 = _kmer_set(seq2, kmer)
    intersection = len(k1 & k2)
    union = len(k1 | k2)
    return 1 - intersection / union if union > 0 else 1.0