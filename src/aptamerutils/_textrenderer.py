import numpy as np
import random
import colorsys
from collections import Counter
from typing import List, Set
from ._dist import _kmer_set

def union_topk(sets:List[Set[str]], topk:int = 3):
    elem_counter = Counter()
    for s in sets:
        for elem in s:
            elem_counter[elem] += 1
            
    return elem_counter.most_common(topk)

def union_fraction(sets:List[Set[str]], m:float = 0.8):
    assert (m > 0 and m < 1)
    threshold = m * len(sets)
    element_counter = Counter()
    for s in sets:
        for elem in s:
            element_counter[elem] += 1
    return [elem for elem, count in element_counter.items() if count > threshold]

def _colorGenerator(n):
    colors = []
    for i in range(n):
        h = i / n
        l = 0.5
        s = 0.9
        rgb = colorsys.hls_to_rgb(h, l, s)
        rgb_255 = tuple(int(255 * x) for x in rgb)
        colors.append('#{:02X}{:02X}{:02X}'.format(*rgb_255))
    random.shuffle(colors)
    return colors

def _drawSequence(ax, seq, position:list, drawposition = (0,0), color = None, showLoop = False, space = 0.04, fontsize = 14):
    color = np.random.rand(3) if color == None else color
    pattern = seq.getStemLoopMap() if showLoop else "x" * len(seq)
    for i in range(len(seq)):
        if i in position:
            ax.text(drawposition[0] + i * space, drawposition[1], seq[i], fontsize = fontsize, color = color, fontweight = "bold" if pattern[i] == "." else "normal")
        else:
            ax.text(drawposition[0] + i * space, drawposition[1], seq[i], fontsize = fontsize, color = "black", fontweight = "bold" if pattern[i] == "." else "normal")
        i += 1
    return ax

def _drawSequenceSet(ax, seqs, position:list, preText:list = None, endText:list = None, colors:list = None, showLoop = False, spaceH = 0.04, spaceV = 0.06, fontsize = 14):
    ax.set_axis_off()
    ax.invert_yaxis()
    posX, posY = 0, 0
    maxPosX = 0
    seqs_ = list(seqs.seqs.keys())
    if (not preText is None):
        length = max([len(text) for text in preText])
    for i in range(len(seqs_)):
        if (not preText is None) and (len(preText) > i):
            ax.text(posX, posY, preText[i], fontsize = fontsize, color = "black")
            posX += spaceH * length
        color = colors[i] if ((not colors is None) and (i < len(colors))) else None
        _drawSequence(ax, seqs_[i], position[i], drawposition = (posX, posY), color = color, showLoop = showLoop, space = spaceH, fontsize = fontsize)
        posX += (len(seqs_[i]) + 1) * spaceH
        if (not endText is None) and (len(endText) > i):
            ax.text(posX, posY, endText[i], fontsize = fontsize, color = "black")
        posX += spaceH * len(endText[i])
        posY += spaceV
        if maxPosX < posX: maxPosX = posX
        posX = 0
    ax.set_xlim(-0.1, maxPosX)
    ax.set_ylim(posY + 0.1, -0.1)
    return ax

def getDisplayClusterFigure(ax, seqlist, displaykmerunion:int = None, showLoop = False, header = "", end = "", fontsize = 5, featureminfrac = 0.8):
    from .seqlist import SeqList
    from .sequence import Sequence
    
    preText, endText, clusterSets, positions, colors = [], [], [], [], []
    seqs = SeqList()
    unique_cluster = seqlist.getUniqueClusters()
    colors_ = _colorGenerator(len(unique_cluster))
    for i in range(len(unique_cluster)):
        currentClusterSet = []
        clusterCount = 0
        for item in seqlist.seqs.items():
            if isinstance(item[1]["label"], dict):
                cluster = item[1]["label"]["cluster"]
                if cluster == unique_cluster[i]:
                    seqs.append(Sequence(header + str(item[0]) + end))
                    preText.append(f"Family{i}")
                    endText.append(f"Count:{item[1]['Count']}")
                    currentClusterSet.append(_kmer_set(str(item[0]), displaykmerunion))
                    clusterCount += 1
            else:
                raise RuntimeError("Sequence clusters is not defined; Use Seqlist.getClustersLabeled first")
        for j in range(clusterCount):
            clusterSets.append(union_fraction(currentClusterSet, featureminfrac))
            colors.append(colors_[i])
    positions = seqs.getFeaturePositionSet(clusterSets)
    ax = _drawSequenceSet(ax, seqs, positions, preText, endText, colors = colors, fontsize = fontsize, showLoop = showLoop)
    return ax