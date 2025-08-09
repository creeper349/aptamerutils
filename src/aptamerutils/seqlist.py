import os
import re
import json
import regex
import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from .sequence import Sequence
from ._distmap import DistMatrix
from ._textrenderer import *
from .find import Find
from typing import cast, Literal

class SeqList():
    def __init__(self):
        self.seqs = {}
        self.order = []
        
    def __len__(self):
        return len(self.order)
    
    def __getitem__(self, key):
        if isinstance(key, int):
            seq = self.order[key]
            labels = self.seqs[seq]
            return seq, labels
        else:
            key = Sequence(key)
            return key, self.seqs[key]
    
    def __setitem__(self, index, seq):
        if not isinstance(seq, Sequence):
            seq = Sequence(seq)
        self.order.insert(index, seq)
        if (seq not in self.seqs):
            self.seqs[seq] = {"Count": 1, "label": None}
        else:
            self.seqs[seq]["Count"] += 1
        
    def __delitem__(self, index):
        seq = self.order[index]
        del self.order[index]
        if self.seqs[seq]["Count"] == 1:
            del self.seqs[seq]
        else:
            self.seqs[seq]["Count"] -= 1
        
    def __iter__(self):
        return iter(self.order)
    
    def __contains__(self, item):
        if not isinstance(item, Sequence):
            item = Sequence(item)
        return item in self.seqs
    
    def __repr__(self):
        return str(self.seqs)
    
    def __str__(self):
        return str(self.seqs)
        
    def append(self, seq, label = None):
        if not isinstance(seq, Sequence):
            seq = Sequence(seq)
        self.order.append(seq)
        if (seq not in self.seqs):
            self.seqs[seq] = {"Count": 1, "label": label}
        else:
            self.seqs[seq]["Count"] += 1
        return self
    
    def addTwoEnds(self, header = "", end = ""):
        lst = list(self.seqs.items())
        newseqlist = SeqList()
        for item in lst:
            newseqlist.seqs[Sequence(header + str(item[0]) + end)] = item[1]
            newseqlist.order += [Sequence(header + str(item[0]) + end)] * item[1]["Count"]
        return newseqlist
            
    def fromList(self, lst:list[Sequence], label = None):
        for seq in lst:
            self.append(seq, label = label)
        return self
    
    def fromLinesInFile(self, path:str, fromline:int = 0):
        lines = []
        with open(path, "r") as f:
            lines_ = f.readlines()
            for i in range(len(lines_)):
                if not i < fromline:
                    lines.append(lines_[i].strip())
        return self.fromList(lines, os.path.basename(path))
    
    def fromfastq(self, path:str):
        with open(path, "r") as f:
            lines_ = f.readlines()
            lines = []
            for i in range(len(lines_)):
                if lines_[i].strip() == "+":
                    lines.append(lines_[i - 1].strip())
        return self.fromList(lines, os.path.basename(path))
    
    def fromfastqFolder(self, dir:str):
        all_files = os.listdir(dir)
        fastq_file = [f for f in all_files if re.fullmatch(r".*\.fastq", f, re.IGNORECASE)]
        for file in fastq_file:
            self.fromfastq(os.path.join(dir, file))
        return self
    
    def fromJson(self, path:str):
        with open(path, "r") as f:
            dict_ = json.load(f)
        lst, seqs = [], {}
        for item in list(dict_.items()):
            lst += [item[0]] * item[1]["Count"]
            seqs[Sequence(item[0])] = {"Count":item[1]["Count"], "label":item[1]["label"]}
        self.seqs = seqs
        self.order = lst
        return self
    
    def clear(self):
        self.seqs = {}
        self.order = []
        return self
    
    def saveJson(self, path:str = "save.json"):
        seqs = {}
        for item in list(self.seqs.items()):
            seqs[str(item[0])] = item[1]
        json.dump(seqs, open(path, "w"), indent = 4)
        
    def saveLines(self, path:str = "save.txt"):
        seqs = []
        for item in list(self.seqs.items()):
            seqs.append(f"Sequence:{str(item[0])}  Count:{item[1]['Count']}   Label:{item[1]['label']}")
        with open(path, "w") as f:
            for seq in seqs:
                f.write(seq + "\n")
                
    def trimTwoEnds(self, start, end, fixed_length = None, fixed_length_tol = 0):
        trimmed = {}
        sequences = []
        for seq in self.order:
            seq = cast(Sequence, seq)
            tseq = seq.trimTwoEnds(start, end, fixed_length, fixed_length_tol)
            sequences.append(tseq)
            if not tseq is None:
                if not tseq in trimmed:
                    trimmed[tseq] = {"Count": 1, "label": self.seqs[seq]["label"]}
                else:
                    trimmed[tseq]["Count"] += 1
        self.seqs = trimmed
        self.order = sequences
        return self
    
    def trimTwoEndsWithLength(self, start, fixed_length):
        trimmed = {}
        sequences = []
        for seq in self.order:
            seq = cast(Sequence, seq)
            tseq = seq.trimTwoEndsWithLength(start, fixed_length)
            if not tseq is None:
                sequences.append(tseq)
                if not tseq in trimmed:
                    trimmed[tseq] = {"Count": 1, "label": self.seqs[seq]["label"]}
                else:
                    trimmed[tseq]["Count"] += 1
        self.seqs = trimmed
        self.order = sequences
        return self
    
    def trimWithFuzzyPattern(self, start:str, end:str, fixed_length:int, fixed_length_tol = 0, pre_ed_tol = 1, end_ed_tol = 1):
        header_pattern = f"(?b:({start}){{s<={pre_ed_tol}}})"
        end_pattern = f"(?b:({end}){{s<={end_ed_tol}}})"
        middle_pattern = f"([ACGT]{{{fixed_length - fixed_length_tol},{fixed_length + fixed_length_tol}}})"
        pattern = regex.compile(header_pattern + middle_pattern + end_pattern)
        trimmed = {}
        sequences = []
        for seq in self.order:
            seq = cast(Sequence, seq)
            tseq = seq.trimWithFuzzyPattern(pattern)
            if not tseq is None:
                sequences.append(tseq)
                if not tseq in trimmed:
                    trimmed[tseq] = {"Count": 1, "label": self.seqs[seq]["label"]}
                else:
                    trimmed[tseq]["Count"] += 1
        self.seqs = trimmed
        self.order = sequences
        return self
    
    def generateDistMap(self, dist_func:Literal["ed", "ncd", "kmer"] = "kmer", parallel = False, **kwargs):
        if parallel:
            return DistMatrix.fromSequencesParallel(list(self.seqs.keys()), dist_func, **kwargs)
        else:
            return DistMatrix.fromSequences(list(self.seqs.keys()), dist_func, **kwargs)
        
    def getClustersLabeled(self, clusters):
        keys = list(self.seqs.keys())
        for i in range(len(keys)):
            flielabel = self.seqs[keys[i]]["label"]
            self.seqs[keys[i]]["label"] = {"cluster": int(clusters[i]), "filelabel": flielabel}
        return self
    
    def getFeaturePositionSet(self, kmerSets:list):
        positions = []
        keys = list(self.seqs.keys())
        for i in range(len(keys)):
            positions.append(list(keys[i].findFeaturePosList(kmerSets[i])))
        return positions
    
    def getUniqueClusters(self):
        uniqueclusters = []
        for value in self.seqs.values():
            if isinstance(value["label"], dict):
                cluster = value["label"]["cluster"]
                if not cluster in uniqueclusters:
                    uniqueclusters.append(cluster)
        return uniqueclusters
    
    def getUniqueSeqs(self):
        return list(self.seqs.keys())
    
    def getkmerCount(self, kmer:int) -> dict:
        count = {}
        items = list(self.seqs.items())
        for i in range(len(items)):
            kmers = cast(Sequence, items[i][0]).getkmerCount(kmer)
            for part in kmers.keys():
                count[part] = count.get(part, 0) + items[i][1]["Count"] * kmers.get(part, 0)
        return count
    
    def getClusterTotalCount(self, cluster_id):
        count = 0
        for item in list(self.seqs.items()):
            if isinstance(item[1]["label"], dict):
                cluster = item[1]["label"]["cluster"]
                if cluster == cluster_id:
                    count += item[1]["Count"]
            else:
                raise RuntimeError("Sequence clusters is not defined; Use Seqlist.getClustersLabeled first")
        return count
    
    def getClusterSeqs(self, cluster_id):
        seqs = SeqList()
        for item in list(self.seqs.items()):
            if isinstance(item[1]["label"], dict):
                cluster = item[1]["label"]["cluster"]
                if cluster == cluster_id:
                    seqs.seqs[item[0]] = item[1]
                    seqs.order += [item[0]] * item[1]["Count"]
            else:
                raise RuntimeError("Sequence clusters is not defined; Use Seqlist.getClustersLabeled first")
        return seqs
    
    def getClusterFeature(self, cluster_id, featureminfrac = 0.8, kmer = 5):
        from ._dist import _kmer_set
        from ._textrenderer import union_fraction
        
        sets = []
        clusterseqs = self.getClusterSeqs(cluster_id)
        for seq in clusterseqs:
            sets.append(_kmer_set(str(seq), kmer))
        return union_fraction(sets, featureminfrac)
    
    def sortbyCount(self, topk:int = None):
        items = list(self.seqs.items())
        counts = [item[1]["Count"] for item in items]
        argsort = np.argsort(counts)[::-1]
        newdict, newlist = {}, []
        if topk is None:
            for arg in argsort:
                newdict[items[arg][0]] = items[arg][1]
                newlist += [items[arg][0]] * items[arg][1]["Count"]
            self.seqs, self.order = newdict, newlist
            return self
        else:
            seqlist = SeqList()
            for arg in argsort[:topk]:
                newdict[items[arg][0]] = items[arg][1]
                newlist += [items[arg][0]] * items[arg][1]["Count"]
            seqlist.seqs, seqlist.order = newdict, newlist
            return seqlist
        
    def sortbyCluster(self):
        newdict, newlist = {}, []
        clusters = self.getUniqueClusters()
        items = list(self.seqs.items())
        for cluster in clusters:
            for i in range(len(items)):
                if isinstance(items[i][1]["label"], dict):
                    if items[i][1]["label"]["cluster"] == cluster:
                        newdict[items[i][0]] = items[i][1]
                        newlist += [items[i][0]] * items[i][1]["Count"]
                else:
                    raise RuntimeError("Sequence clusters is not defined; Use Seqlist.getClustersLabeled first")
        self.seqs, self.order = newdict, newlist
        return self
    
    def drawText(self, filename:str = "save.pdf", displaykmerfeature:int = None, 
                 showLoop = False, header = "", end = "", fontsize = 5, 
                 featureminfrac = 0.8, style:Literal["A", "B"] = "B", 
                 figuresize = (8, 12), **kwargs):
        fig, ax = plt.subplots(figsize = figuresize)
        ax = getDisplayClusterFigure(ax, self, displaykmerfeature, showLoop, header, end, fontsize, featureminfrac, style = style, **kwargs)
        plt.savefig(filename)
        return ax
    
    def findCustomSeqCombination(self, pattern: Find):
        keys = list(self.seqs.keys())
        match = SeqList()
        for key in keys:
            if key.matchCustomSeqCombination(pattern):
                match.seqs[key] = self.seqs[key]
                match.order += [key] * self.seqs[key]["Count"]
        return match
    
if __name__ == "__main__":
    s = SeqList().fromfastqFolder(".\Tobph7").trimTwoEndsWithLength("GACGAC", 30)
    clusters = s.generateDistMap().createPosMap().getCluster(eps = 0.6)
    s.getClustersLabeled(clusters)
    s = s.sortbyCount(topk = 50)
    print(s)
    print(s.sortbyCluster())
    print(s.findCustomSeqCombination(Find("ACCTGTA") & Find("GGGGTCC")))
    s.drawText("./proj/save.pdf", displaykmerfeature = 5, showLoop = True, featureminfrac = 0.8)