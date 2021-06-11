import argparse
import csv
import os
import re
from collections import defaultdict
from itertools import combinations

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
    outFile = open(out+'.csv', 'w')
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
                    if ptms:
                        ptmList.append('{}'.format(ptms))
                    else:
                        ptmList.append('{}'.format(str('')))
                else:
                    ptmList.append('{}'.format(str('')))
            print(ptmList) 
            print(combPTM)
            outFile.write('{},{},{},'.format(kmer, startPos, endPos))
            if combPTM: ##calculate combinations needed 
                combins = len(combPTM)
                totalN=(combins*combins)
                combTrac += totalN
                outFile.write(','.join(ptmList)+',{}\n'.format(totalN))
            else:
                outFile.write(','.join(ptmList)+',{}\n'.format(1))
                #combTrac += 1
    outFile.close()
    print('OUT WRITTEN {} WITH TOTAL COMB {}'.format(out, combTrac))      


##create combination vector
#kmer = list('SLPTPPTREPK')
def combPeptides(posList, kmer):
    #posList = ['', '', '', '', '', 'Gly|N6-acK', 'N6-acK', '', '', '', '']
    #kmer = ['P', 'G', 'G', 'G', 'N', 'K', 'K', 'I', 'E', 'T', 'H'] 
    print(posList, kmer)
    posDict=[]
    for n, pos in enumerate(posList):
        if pos:
            for ptm in pos.split('|'):
                posDict.append((ptm, n))
    nComb = len(posDict)+1
    posIndex = []
    outKmer = set()
    for i in range(1, nComb):
        for sub in combinations(posDict, i): 
            #print(sub)
            kmerCp = kmer.copy()
            posIndex.append(sub)
            for item in sub:
                ptm, pos = item
                #print(ptm, pos)
                makePTM=kmerCp[pos]
                if not re.search("\\[", makePTM): #if a ptm has already been assigned at the same position do nothing
                    motif='{}[{}]'.format(makePTM, ptm)
                    kmerCp[pos] = motif
            outKmer.add(''.join(kmerCp))
    print('BUILT {} COMBINATIONS'.format(len(outKmer)))
    return outKmer


def makePeptidelist(out):
    pepOut=open('{}_PEPLIST.csv'.format(out), 'w')
    writer = csv.writer(pepOut)
    writer.writerow(['PEP', 'START', 'END', 'BASE'])
    with open(out+'.csv') as pepFile:
        reader = csv.reader(pepFile)
        next(reader)
        for row in reader:
            pep, st, end= row[:3]
            writer.writerow([pep, st, end, '*'])
            ptmCheck = '\t'.join(row[3:14])
            if ptmCheck:
                ptmList = row[3:14]
                print('MAKING PTM COMBINATIONS {} FOR {}'.format(ptmCheck, pep))### implement bool check for ptms here ?
                kmer = list(pep)
                outKmer = combPeptides(posList=ptmList, kmer=kmer)
                for pep in outKmer:
                    writer.writerow([pep, st, end, ''])
    print('FINISHED WRITE TO {}'.format(pepOut))
    pepOut.close()

#makePeptidelist(out='outputs/peptideList.csv')

#"-ptm", "data/PTMList_DANIEL.txt"
def main():
    parser = argparse.ArgumentParser('A script to make user defined peptide mers with an option for overlap, additionally tag ptms')
    parser.add_argument('-fasta', required=True, help='FASTA file of protein sequence to be processed')
    parser.add_argument('-ptm', required=False, help='A tab delimited PTM list should have pos aminoacid and ptm as first three columns')
    parser.add_argument('-mer', required=True, help='total length of the kmer', type=int)
    parser.add_argument('-overlap', required=True, help='overlap needed', type=int)
    parser.add_argument('-out', required=True, help='file to write output')
    args = parser.parse_args()
    print(args)
    fasta = args.fasta
    mer = args.mer
    overlap = args.overlap
    out = args.out
    if args.ptm is not None:
        ptm = args.ptm
        posDict, Seq = parseFasta(fasta)
        posDict=parsePTMs(posDict=posDict, ptm=ptm)
        makePeptides(posDict, mer, overlap, out, Seq)
        makePeptidelist(out)
    else:
        posDict, Seq = parseFasta(fasta)
        makePeptides(posDict, mer, overlap, out, Seq)

if __name__ == '__main__':main()
