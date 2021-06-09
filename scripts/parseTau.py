import argparse
import csv
import os
import re
from collections import defaultdict

from Bio import SeqIO


### parse the fasta files
def parseFasta(fasta):
    if os.path.exists(fasta):
        Seq=[str(i.seq) for i in SeqIO.parse(fasta, 'fasta')][0]
    ## create a tau dict pos id is key amino acid is string value
    print('FOUND 1 FASTA WITH {} AA'.format(len(Seq)))
    posDict = defaultdict()
    for n, aa in enumerate(Seq):
        posDict[(n+1, aa)] = set()
    return posDict, Seq


## make ptm dict
def parsePTMs(ptm, posDict):
    with open(ptm, 'r') as ptmFile:
        reader = csv.reader(ptmFile, delimiter='\t')
        next(reader)
        for row in reader:
            n, aa, p = row[:3]
            #print(len(row))
            key=(int(n), aa)
            print(key)
            if key in posDict:
                posDict[key].add(p)
    return posDict


## main call  
def makePeptides(posDict, mer, overlap, out, Seq):
    ##intitalize params
    combTrac = 0
    pepLen = mer
    span = pepLen-overlap
    outFile = open(out, 'w')
    header = ['{}Mer'.format(pepLen), 'Start', 'End']
    header.extend(['P{}'.format(i) for i in range(1, pepLen+1)])
    header.extend(['NCombinations'])
    outFile.write(','.join(header)+'\n')
    for i in range(0, len(Seq), span):
        startPos = i+1
        endPos = startPos+(pepLen-1)
        kmer = Seq[i:i+pepLen]
        if len(kmer) == pepLen:
            #break
            print(startPos, 'to ', endPos, ' ', Seq[i:i+pepLen], i)
            ptmList =[]
            combPTM = []
            for j in range(startPos, endPos+1): ##change pointer 
                key = (int(j), Seq[j-1])
                if key in posDict:
                    ptmOut=[i for i in posDict.get(key)]
                    combPTM.extend(ptmOut)
                    ptms='|'.join(ptmOut)
                    ptmList.append('{} {}'.format(str(j), ptms))
                else:
                    ptmList.append('{}'.format(str('')))
            print(ptmList) 
            print(combPTM)
            outFile.write('{},{},{},'.format(kmer, startPos, endPos))
            if combPTM: ##calculate combinations needed 
                combins = len(combPTM)
                totalN=combins*combins
                totalN += 1
                combTrac += totalN
                outFile.write(','.join(ptmList)+',{}\n'.format(totalN))
            else:
                outFile.write(','.join(ptmList)+',{}\n'.format(1))
    outFile.close()
    print('OUT WRITTEN {} WITH TOTAL COMB {}'.format(out, combTrac))      


def main():
    parser = argparse.ArgumentParser('A script to make user defined peptide mers with an option for overlap, additionally tag ptms')
    parser.add_argument('-fasta', required=True, help='FASTA file of protein sequence to be processed')
    parser.add_argument('-ptm', required=True, help='A tab delimited PTM list should have pos aminoacid and ptm as first three columns')
    parser.add_argument('-mer', required=True, help='total length of the kmer', type=int)
    parser.add_argument('-overlap', required=True, help='overlap needed', type=int)
    parser.add_argument('-out', required=True, help='file to write output')
    args = parser.parse_args()
    fasta = args.fasta
    ptm = args.ptm
    mer = args.mer
    overlap = args.overlap
    out = args.out
    posDict, Seq = parseFasta('data/tauMain.fasta')
    posDict=parsePTMs(posDict=posDict, ptm='data/PTMList_DANIEL.txt')
    makePeptides(posDict, mer, overlap, out, Seq)

if __name__ == '__main__':main()