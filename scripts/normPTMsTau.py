import csv
import os
import re
from collections import defaultdict

from Bio import SeqIO

## Daniel PTM list
S3P = [[46,"S","pSer","?"],[61,"S","pSer","?"],[214,"S","pSer","SGK1"],[396,"S","pSer","PHF-tau"],[502,"S","pSer","?"],[508,"S","pSer","?"],[512,"S","pSer","?"],[515,"S","pSer","PDPK1 and TTBK1"],[515,"S","pSer","PDPK1 and TTBK1"],[516,"S","pSer","PDPK1 and TTBK1"],[519,"S","pSer","CK1 PDPK1 and TTBK1"],[531,"S","pSer","PKA"],[552,"S","pSer","PDPK1"],[554,"S","pSer","PHK"],[579,"S","pSer","MARK1, MARK2, MARK3, MARK4, BRSK1, BRSK2 and PHK"],[602,"S","pSer","PHK"],[610,"S","pSer","?"],[622,"S","pSer","PHK"],[641,"S","pSer","?"],[669,"S","pSer","PHK"],[673,"S","pSer","?"],[713,"S","pSer","CK1 and PDPK1"],[721,"S","pSer","CK1 and PDPK1"],[726,"S","pSer","?"],[733,"S","pSer","CaMK2 and TTBK1"],[739,"S","pSer","PDPK1 and TTBK1"]]
Kac = [[480,"K","N6-acK","alternate"],[542,"K","N6-acK","?"],[576,"K","N6-acK","alternate"],[615,"K","N6-acK","alternate"],[628,"K","N6-acK","alternate"],[634,"K","N6-acK","alternate"],[638,"K","N6-acK","alternate"],[648,"K","N6-acK","alternate"],[660,"K","N6-acK","alternate"],[664,"K","N6-acK","alternate"],[686,"K","N6-acK","alternate"],[702,"K","N6-acK","alternate"]]
S3P.extend(Kac)
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


recordDict = SeqIO.to_dict(SeqIO.parse('data/tauIsoforms.fasta', "fasta"))
fastaDict= defaultdict(str)
for k, v in recordDict.items():
    k_ = k.replace('|', '_')+'.fasta'
    fastaDict[convNames.get(k)] = k_
    fName = os.path.join('data','{}'.format(k_))
    print(fName)
    with open(fName, "w") as output_handle:
        SeqIO.write(v, output_handle, "fasta")

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

ddConv={'pSer':'pSer', 'N6-acK':'Acetyl'}
for dList in [S3P, Kac]:
    for row in dList:
        print(ptm)
        posParsed, aaParsed, ptm = row[:3]
        poskey = '{}_{}'.format(aaParsed,posParsed)
        if ptm in ddConv:
            ptmParsed = ddConv.get(ptm)
            ptmId = '{}{}_{}'.format(aaParsed, posParsed, ptmParsed)
        makeRow = [ptmId, posParsed, aaParsed, ptmParsed, '', '']
        isoVec = ['' for _ in range(len(isoNames))]
        if poskey in positionDict:
            isoList=positionDict.get(poskey)
            for iso in isoList:
                getInd=isoNames.index(iso)
                isoVec[getInd] = '*'
        if posParsed <= len(tauSeq):
            if tauSeq[posParsed] == aaParsed:
                makeRow.append('*')
        else:
            makeRow.append(' ')
        makeRow.extend(isoVec)
        writer.writerow(makeRow)
outFile.close()

## did some manual editing of these positions T_69, K_87, S_56
# 6 2N4R_Unaligned
# 7 2N3R
# 8 1N3R_Tau-B
# 9 2N4R_Tau-F
# 10 0N3R_Tau-Fetal
# 11 1N4R_Tau-E
# 12 0N4R_Tau-D
# 13 Tau-G
# 14 Canonical
# 15 Tau-A
ptmFile = 'data/PTMList_Weissling_isoforms_June21_2021.txt'
ptmDict = defaultdict(list)
with open(ptmFile, 'r') as tabFile:
    reader = csv.reader(tabFile, delimiter='\t')
    head = []
    for n, row in enumerate(reader):
        if n == 0:
            head.extend(row)
        else:
            makeKey = row[2]+'_'+row[1]
            for n, item in enumerate(row[6:]):
                n += 6
                if item == '*':
                    key=head[n]
                    ptmDict[key].append(row[1:4])

##write out the ptm lists and commandds
comFile=open('scripts/tauptmCommands.sh', 'w')
length = 15
span = 11
for k, v in ptmDict.items():
    makeFile = os.path.join('data/', k+'.txt')
    if not os.path.exists(makeFile):
        fFile = fastaDict.get(k)
        oFile = k+'_{}mer_{}overlap_PTM'.format(length, span)
        comFile.write('python3 scripts/parseProt.py -fasta data/{} -ptm {} -mer {} -overlap {} -out outputs/{}\n'.format(fFile, makeFile, length, span, oFile))
        outFile = open(makeFile, 'w')
        writer = csv.writer(outFile, delimiter='\t')
        print('TRIGGERING NEW FILE {}'.format(makeFile))
        writer.writerow(['POS', 'AA', 'PTM'])
        for i in v:
            writer.writerow(i)
        outFile.close()
    else:
        outFile = open(makeFile, 'a')
        writer = csv.writer(outFile, delimiter='\t')
        for i in v:
            writer.writerow(i)
        outFile.close()
    outFile.close()
comFile.close()
    



