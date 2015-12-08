from __future__ import division
import json
import numpy as np
import sys
import os
import random

def main(length):
    HindIII = 'AAGCTT'
    BamHI = 'GGATCC'
    total = 1
    
    s = ''.join([random.choice('ATCG') for j in range(length)])
 
    seq = cleanup(s)
 
    Gcount = 0
    pos = seq.find('G',0)
    while not pos == -1:
        Gcount = Gcount + 1
        pos = seq.find('G', pos+1)
 
    gc = 2*Gcount/length
 
    if gc < 0.49 and gc > 0.45:
        sequen = seq
        return sequen
            
    else:
        sequen = main(length)
        return sequen


def checkScaffoldCorrectness(s):
    SacII = 'CCGCGG'
    SalI = 'GTCGAC'
    SacI = 'GAGCTC'
#     BtsCI = 'GGATG'

    forbiddenSeqs = [SacII, SalI, SacI]
    
    for forbiddenSeq in forbiddenSeqs:
        if s.count(forbiddenSeq) > 0:
            return False
        #end if
    #end for
    
    return True
#end def

def cleanup(s):
    aaaaa = 'AAAAA'
    ccccc = 'CCCCC'
    ggggg = 'GGGGG'
    ttttt = 'TTTTT'
    SacII = 'CCGCGG'
    SalI = 'GTCGAC'
    SacI = 'GAGCTC'
    BtsCI = 'GGATG'
    BtsCIComplement = 'CCTAC'

 
    done = False
    while not done:
        if s.count(aaaaa) > 0:
            pos = s.find(aaaaa, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(aaaaa))])
            s = s[:pos] + new + s[pos + len(aaaaa):]
            # print "replace polyA", pos
            continue
        elif s.count(ccccc) > 0:
            pos = s.find(ccccc, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(ccccc))])
            s = s[:pos] + new + s[pos + len(ccccc):]
            # print "replace polyC", pos
            continue
        elif s.count(ggggg) > 0:
            pos = s.find(ggggg, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(ggggg))])
            s = s[:pos] + new + s[pos + len(ggggg):]
            # print "replace polyG", pos
            continue
        elif s.count(ttttt) > 0:
            pos = s.find(ttttt, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(ttttt))])
            s = s[:pos] + new + s[pos + len(ttttt):]
            # print "replace polyT", pos
            continue
        elif s.count(SacII) > 0:
            pos = s.find(SacII, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(SacII))])
            s = s[:pos] + new + s[pos + len(SacII):]
            # print "replace BbvCIB", pos
            continue
        elif s.count(SalI) > 0:
            pos = s.find(SalI, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(SalI))])
            s = s[:pos] + new + s[pos + len(SalI):]
            # print "replace BbvCIB", pos
            continue
        elif s.count(SacI) > 0:
            pos = s.find(SacI, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(SacI))])
            s = s[:pos] + new + s[pos + len(SacI):]
            # print "replace BbvCIB", pos
            continue
        elif s.count(BtsCI) > 0:
            pos = s.find(BtsCI, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(BtsCI))])
            s = s[:pos] + new + s[pos + len(BtsCI):]
            # print "replace BbvCIB", pos
            continue
        elif s.count(BtsCIComplement) > 0:
            pos = s.find(BtsCIComplement, 0)
            new = ''.join([random.choice('ATCG') for j in range(len(BtsCI))])
            s = s[:pos] + new + s[pos + len(BtsCI):]
            # print "replace BbvCIB", pos
            continue
        break
    return s 

