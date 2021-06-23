import csv
import re
from collections import defaultdict
from Bio import SeqIO

#fastaIn = 'data/tauIsoforms.fasta'
#recordDict = SeqIO.to_dict(SeqIO.parse(fastaIn, "fasta"))
recordDict = SeqIO.to_dict(SeqIO.parse('data/tauIsoforms.fasta', "fasta"))
tau = str(recordDict.get('sp|P10636-8|TAU_HUMAN').seq)
kmerMatch = defaultdict(int)
RRepeat = 'PGGG'
NRepeat = ''

pepLen = 3
for i in range(0, len(tau), 1):
    kmer = tau[i:i+pepLen]
    print(kmer)
    if len(kmer) == pepLen:
        matches = re.findall(kmer, tau)
        if len(matches) > 1:
            kmerMatch[kmer]=len(matches)

## create a conv dict
convNames = {'sp|P10636-5|TAU_HUMAN':'2N3R', 
'sp|P10636-4|TAU_HUMAN':'1N3R_Tau-B', 
'sp|P10636-8|TAU_HUMAN':'2N4R_Tau-F', 
'sp|P10636-2|TAU_HUMAN':'0N3R_Tau-Fetal', 
'sp|P10636-7|TAU_HUMAN':'1N4R_Tau-E', 
'sp|P10636-6|TAU_HUMAN':'0N4R_Tau-D', 
'sp|P10636-9|TAU_HUMAN':'Tau-G', 
'sp|P10636|TAU_HUMAN':'Canonical', 
'sp|P10636-3|TAU_HUMAN':'Tau-A'}

## create a conv dict
# convNames = {'sp|P10636-5|TAU_HUMAN/1-410':'2N3R', 
# 'sp|P10636-4|TAU_HUMAN/1-381':'1N3R_Tau-B', 
# 'sp|P10636-8|TAU_HUMAN/1-441':'2N4R_Tau-F', 
# 'sp|P10636-2|TAU_HUMAN/1-352':'0N3R_Tau-Fetal', 
# 'sp|P10636-7|TAU_HUMAN/1-412':'1N4R_Tau-E', 
# 'sp|P10636-6|TAU_HUMAN/1-383':'0N4R_Tau-D', 
# 'sp|P10636-9|TAU_HUMAN/1-776':'Tau-G', 
# 'sp|P10636|TAU_HUMAN/1-758':'Canonical', 
# 'sp|P10636-3|TAU_HUMAN/1-316':'Tau-A'}
##check the aligned fasta
#fastaIn = 'data/tauAligned.fa'
#recordDict = SeqIO.to_dict(SeqIO.parse(fastaIn, "fasta"))

recordDict = SeqIO.to_dict(SeqIO.parse('data/tauIsoforms.fasta', "fasta"))
orginalFasta = SeqIO.to_dict(SeqIO.parse('data/tauIsoforms.fasta', "fasta"))
## get the 2N4R
tau_=orginalFasta.get('sp|P10636-8|TAU_HUMAN')#.seq
tauSeq= ' '+str(tau_.seq)
positionDict = defaultdict(list)
for id, rec in recordDict.items():
    #recId = rec.description
    makeName = convNames.get(id)
    for n, aa in enumerate(rec.seq):
        n += 1
        if aa != "-":
            print('{} AA {} POS FROM {}'.format(aa, n, id))
            positionDict['{}_{}'.format(aa, n)].append(makeName)


## parse the weissling ptms
convDict={'K-Acetyl':"Acetyl",
 'K-GlyGly':'diGLY',
 'K-Methyl':'K-Methyl',
 'K-Phospho':'pLys',
 'R-Methyl':'R-Methyl',
 'S-Phospho':'pSer',
 'T-Phospho':'pThr',
 'Y-Phospho':'pTyr'}
ptmFile = 'data/ptmList_Wessling_et_al.csv'
wTauPtm = set()
outFile=open('data/PTMList_Weissling_isoforms.txt', 'w')# as outFile:
writer=csv.writer(outFile, delimiter='\t')
header=['PTMID', 'POS', 'AA', 'PTM', 'AD', 'CT', '2N4R_Unaligned']
isoNames = [i for i in convNames.values()]
header.extend(isoNames)
writer.writerow(header)
with open(ptmFile) as csvFile:
    reader = csv.reader(csvFile)
    next(reader)
    for row in reader:
        if row[0]:
            ptmId, pos, ptm, ad, ct = row
            wTauPtm.add(ptm)
            if 'or' in pos:
                #print(row)
                for p in pos.split('or'):
                    posParsed = int(re.findall('[0-9]+', p)[0])
                    aaParsed = re.findall('[A-Z]+', p)[0]
                    key='{}-{}'.format(aaParsed,ptm)
                    convKey = convDict.get(key)
                    poskey = '{}_{}'.format(aaParsed,posParsed)
                    makeRow = [ptmId, posParsed, aaParsed, convKey, ad, ct]
                    ## check pos
                    isoVec = ['' for _ in range(len(isoNames))]
                    if poskey in positionDict:
                        isoList=positionDict.get(poskey)
                        for iso in isoList:
                            getInd=isoNames.index(iso)
                            isoVec[getInd] = '*'
                    if tauSeq[posParsed] == aaParsed:
                        makeRow.append('*')
                    else:
                        makeRow.append(' ')
                    makeRow.extend(isoVec)
                    writer.writerow(makeRow)
                        
            else:
                posParsed = int(re.findall('[0-9]+', pos)[0])
                aaParsed = re.findall('[A-Z]+', pos)[0]
                key='{}-{}'.format(aaParsed,ptm)
                convKey = convDict.get(key)
                poskey = '{}_{}'.format(aaParsed,posParsed)
                ## check pos
                makeRow = [ptmId, posParsed, aaParsed, convKey, ad, ct]
                ## check pos
                isoVec = ['' for _ in range(len(isoNames))]
                if poskey in positionDict:
                    isoList=positionDict.get(poskey)
                    for iso in isoList:
                        getInd=isoNames.index(iso)
                        isoVec[getInd] = '*'
                if tauSeq[posParsed] == aaParsed:
                    makeRow.append('*')
                else:
                    makeRow.append(' ')
                makeRow.extend(isoVec)
                writer.writerow(makeRow)
outFile.close()