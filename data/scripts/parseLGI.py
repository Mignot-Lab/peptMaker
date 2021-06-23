import csv
import os
import re
from collections import defaultdict
from Bio import SeqIO

lgiFasta = Seq=[str(i.seq) for i in SeqIO.parse('data/LGI1_Human.fasta', 'fasta')][0]
convDict={'A':"N6-acK",
'G':"Gly",
'M':"K-Methyl",
'O':"Gly",
'U':"Ub"}
phos = {'S':"pSer",#S
'T':"pThr",#T
'Y':"pTyr"}
#O	O-(mucin type)glycosylation
#U	ubiquitylation
ptmDict = defaultdict(lambda:defaultdict(set))
with open('data/LGI1-PTM.csv') as inFile:
    reader = csv.reader(inFile)
    next(reader)
    for row in reader:
        print(row)
        core=row[3]
        pos=row[9:]
        n =0
        for i, j in zip(list(core), pos):
            n += 1
            if j:
                if i in phos:
                    j = phos.get(i)
                elif j in convDict:
                    j = convDict.get(j)
                else:
                    j = ''
            else:
                j = ''
            ptmDict[core][n].add(j)

outFile=open('data/LGI_PTMS.txt', 'w')
writer = csv.writer(outFile, delimiter='\t')
writer.writerow(['POS', 'AA', 'PTM'])
for seq, val in ptmDict.items():
    if seq in lgiFasta:
        getIndex  = lgiFasta.find(seq)
        end=getIndex+len(seq)
        kmer = lgiFasta[getIndex:end]
        #print(kmer, seq)
        assert kmer == seq
        for pos, item in val.items():
            for ptm in item:
                if ptm:
                    print(ptm, kmer, pos, getIndex)
                    actPos = getIndex+pos
                    makeRow = [actPos, seq[pos-1], ptm]
                    print(makeRow)
                    writer.writerow(makeRow)
outFile.close()

