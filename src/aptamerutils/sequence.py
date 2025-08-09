import RNA
import re
import numpy as np
from .find import Find
from typing import Literal
from ._dist import *

class Sequence:

    def __init__(self, seq:str):
        self.seq = seq.upper()
        
    def __str__(self):
        return self.seq
    
    def __repr__(self):
        return f"Sequence(\"{self.seq}\")"
    
    def __len__(self):
        return len(self.seq)
    
    def __getitem__(self, key):
        result = self.seq[key]
        return Sequence(result) if isinstance(result, str) else result
    
    def __contains__(self, item):
        if isinstance(item, Sequence):
            item = str(item)
        return item in self.seq
    
    def __iter__(self):
        return iter(self.seq)
    
    def __eq__(self, other):
        if not isinstance(other, Sequence):
            other = Sequence(other)
        return self.seq == other.seq
    
    def __hash__(self):
        return hash(self.seq)
    
    def __add__(self, other):
        if isinstance(other, Sequence):
            return Sequence(self.seq + other.seq)
        return Sequence(self.seq + str(other))
    
    def __radd__(self, other):
        if isinstance(other, Sequence):
            return Sequence(other.seq + self.seq)
        return Sequence(str(other) + self.seq)
    
    def __mul__(self, n):
        return Sequence(self.seq * n)
    
    def __rmul__(self, n):
        return Sequence(n * self.seq)
    
    def __getattr__(self, name):
        seq = object.__getattribute__(self, 'seq')
        if hasattr(seq, name):
            attr = getattr(self.seq, name)
            if callable(attr):
                def wrapper(*args, **kwargs):
                    result = attr(*args, **kwargs)
                    return Sequence(result) if isinstance(result, str) else result
                return wrapper
            return attr
        raise AttributeError(f"'Sequence' object has no attribute '{name}'")
    
    def isDNA(self):
        return ("T" in self.seq)
    
    def isRNA(self):
        return ("U" in self.seq)
    
    def toDNA(self):
        return Sequence(self.seq.replace("U", "T"))
    
    def toRNA(self):
        return Sequence(self.seq.replace("T", "U"))
    
    def getComplement(self):
        newStrand = ""
        for char in self.seq:
            if self.isDNA():
                newChar = "A" if char == "T" else "T" if char == "A" else "G" if char == "C" else "C" if char == "G" else ""
            elif self.isRNA():
                newChar = "A" if char == "U" else "U" if char == "A" else "G" if char == "C" else "C" if char == "G" else ""
            newStrand += newChar
            
        return Sequence(newStrand)
    
    def reverse(self):
        return Sequence(self.seq[::-1])
    
    def reverseComplement(self):
        return self.getComplement().reverse()
    
    def getStemLoopMap(self):
        pattern = RNA.fold(self.seq)[0]
        return pattern
    
    def getLoopPos(self):
        pattern = self.getStemLoopMap()
        pos = []
        for i in range(len(self)):
            if pattern[i] == ".":
                pos.append(i)
        return pos
    
    def getkmerSet(self, kmer:int) -> set:
        return set([self[i:i + kmer] for i in range(len(self) - kmer + 1)])
    
    def getkmerCount(self, kmer:int) -> dict:
        counts = {}
        for i in range(len(self) - kmer + 1):
            part = self[i: i + kmer]
            counts[part] = counts.get(part, 0) + 1
        return counts
    
    def getTmEstimation(self):
        if len(self.seq) < 14:
            a = self.seq.count('A')
            t = self.seq.count('T')
            g = self.seq.count('G')
            c = self.seq.count('C')
            return 2 * (a + t) + 4 * (g + c)
        else:
            gc_frac = self.GC_Content()
            return 64.9 + 41 * (gc_frac - 0.16)
        
    def getMW(self):
        weight_table = {
            'A': 347.224,
            'T': 338.210,
            'G': 363.223,
            'C': 323.198,
            'U': 324.182    # Data from ChemSpider
        }
        return sum(weight_table[base] for base in self.seq)\
                    - 18.015 * (len(self) - 1) - int(self.isDNA()) * 16 * len(self)
    
    def GC_Content(self):
        g = self.seq.count("G")
        c = self.seq.count("C")
        return (g + c) / len(self) if self else 0
    
    def findall(self, sub: str):
        positions = []
        i = 0
        while i <= len(self) - len(sub):
            if self[i:i+len(sub)] == sub:
                positions.append(list(range(i, i + len(sub))))
                i += len(sub)
            else:
                i += 1
        return positions
    
    def findFeaturePosList(self, kmerSet:set[str]):
        positionlist = []
        for kmer in kmerSet:
            positions = self.findall(kmer)
            for position in positions:
                positionlist += position
        return np.unique(positionlist)
    
    def matchCustomSeqCombination(self, pattern:Find):
        return pattern.match(self.seq)
    
    def trimTwoEnds(self, start, end, fixed_length = None, fixed_length_tol = 0):
        pattern = r"%s([ATCGU]+)%s" % (re.escape(str(start)), re.escape(str(end)))
        match = re.search(pattern, self.seq)
        if match:
            target = match.group(1)
            if (not fixed_length) or (fixed_length - fixed_length_tol <= len(target) <= fixed_length + fixed_length_tol):
                return Sequence(target)
        return None
    
    def trimTwoEndsWithLength(self, start, fixed_length):
        pattern = r"%s([ATCGU]{%d})" % (re.escape(str(start)), fixed_length)
        match = re.search(pattern, self.seq)
        if match:
            return Sequence(match.group(1))
        else:
            return None
        
    def trimWithFuzzyPattern(self, pattern):
        match = pattern.search(self.seq)
        if match:
            return Sequence(match.group(2))
        else:
            return None
    
    def distWith(self, other, dist: Literal["ed", "ncd", "kmer"] = "kmer", **kwargs):
        dist_matrics = editing_distance if dist == "ed" else ncd if dist == "ncd" else jaccard
        return dist_matrics(self.seq, str(other), **kwargs)
    
if __name__ == "__main__":
    a = Sequence("atcg")
    b = Sequence("GGG")
    seq = Sequence("GACGACCATGGAGGTATACAATCTGTCGTC")

    print(a.upper(), a.isdecimal())
    print(a + b)
    print("TC" in a)
    print(a[1:3])
    print(len(a))
    print(a.replace("T", "U"))
    print(a == "ATCG")
    print(3 * b, b * 3)
    print(seq.isDNA(), seq.reverseComplement())
    print(seq.getStemLoopMap(), seq.trimTwoEndsWithLength("GACGAC", 12))
    print(seq.matchCustomSeqCombination(Find(Sequence("GACGAC")) & Find("GTCGTC") & Find(b)))
    print(Sequence("ATGC").getMW())
    print(Sequence("ATGCAGC").distWith("ATGCATG", "kmer", kmer = 4))