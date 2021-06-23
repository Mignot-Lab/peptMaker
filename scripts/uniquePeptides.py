import csv
import os
import re
from collections import defaultdict, Counter

# fileList = [
# 'outputs/1N4R_Tau-E_15mer_11overlap_PTM_PEPLIST.csv',
# 'outputs/0N4R_Tau-D_15mer_11overlap_PTM_PEPLIST.csv',
# 'outputs/2N3R_15mer_11overlap_PTM_PEPLIST.csv',
# 'outputs/Tau-A_15mer_11overlap_PTM_PEPLIST.csv',
# 'outputs/1N3R_Tau-B_15mer_11overlap_PTM_PEPLIST.csv',
# 'outputs/2N4R_Tau-F_15mer_11overlap_PTM_PEPLIST.csv',
# 'outputs/0N3R_Tau-Fetal_15mer_11overlap_PTM_PEPLIST.csv',
# 'outputs/Tau-G_15mer_11overlap_PTM_PEPLIST.csv',
# 'outputs/Canonical_15mer_11overlap_PTM_PEPLIST.csv']
fileList = [
'outputs/1N4R_Tau-E_11mer_10overlap_PTM_PEPLIST.csv',
'outputs/0N4R_Tau-D_11mer_10overlap_PTM_PEPLIST.csv',
'outputs/2N3R_11mer_10overlap_PTM_PEPLIST.csv',
'outputs/Tau-A_11mer_10overlap_PTM_PEPLIST.csv',
'outputs/1N3R_Tau-B_11mer_10overlap_PTM_PEPLIST.csv',
'outputs/2N4R_Tau-F_11mer_10overlap_PTM_PEPLIST.csv',
'outputs/0N3R_Tau-Fetal_11mer_10overlap_PTM_PEPLIST.csv',
'outputs/Tau-G_11mer_10overlap_PTM_PEPLIST.csv',
'outputs/Canonical_11mer_10overlap_PTM_PEPLIST.csv']
peptideDict = defaultdict(set)
peptideCounter = Counter()
for file_ in fileList:
    if os.path.exists(file_):
        fName = os.path.basename(file_)
        fNameClean=fName.split('_')[0]
        print(fNameClean)
        fileHandle = open(file_, 'rt')
        reader = csv.reader(fileHandle)
        next(reader)
        for row in reader:
            print(row[0])
            peptideCounter[fNameClean] += 1
            peptideDict[row[0]].add(fNameClean)

##make n*m matrix
nMatrixCounter = Counter()
isoList = set()
for k, v in peptideDict.items():
    if len(v) == 1:
        iso = v.pop()
        isoList.add(iso)
        nMatrixCounter[(iso, iso)] += 1
    else:
        iso1 = v.pop()
        for iso2 in v:
            if iso1 != iso2:
                isoList.add(iso1)
                isoList.add(iso2)
                nMatrixCounter[(iso1, iso2)] += 1

outFile = open('outputs/sharingMatrix11mer.csv', 'w')
outFile.write('ISO,'+','.join([i for i in isoList]) + '\n')
for iso1 in isoList:
    outFile.write(iso1+',')
    row = []
    for iso2 in isoList:
        makeKey=(iso1, iso2)
        #print(nMatrixCounter.get(makeKey), makeKey)
        if makeKey in nMatrixCounter:
            getCount = nMatrixCounter.get(makeKey)
        else:
            getCount = 0
        row.append(str(getCount))
    print(row)
    outFile.write(','.join(row)+'\n')
outFile.close()
