import csv
import re

from Bio import SeqIO

fastaIn = 'data/tauIsoforms.fasta'
recordDict = SeqIO.to_dict(SeqIO.parse(fastaIn, "fasta"))
tauIso = recordDict.get('sp|P10636|TAU_HUMAN')
tauSeq = str(tauIso.seq)
tau = "-MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQEPESGKVVQEGFLREPGPPGLSHQLMSGMPGAPLLPEGPREATRQPSGTGPEDTEGGRHAPELLKHQLLGDLHQEGPPLKGAGGKERPGSKEEVDEDRDVDESSPQDSPPSKASPAQDGRPPQTAAREATSIPGFPAEGAIPLPVDFLSKVSTEIPASEPDGPSVGRAKGQDAPLEFTFHVEITPNVQKEQAHSEEHLGRAAFPGAPGEGPEARGPSLGEDTKEADLPEPSEKQPAAAPRGKPVSRVPQLKARMVSKSKDGTGSDDKKAKTSTRSSAKTLKNRPCLSPKHPTPGSSDPLIQPSSPAVCPEPPSSPKYVSSVTSRTGSSGAKEMKLKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
S_Phospo = [[46,"S","pSer","?"],[61,"S","pSer","?"],[214,"S","pSer","SGK1"],[396,"S","pSer","PHF-tau"],[502,"S","pSer","?"],[508,"S","pSer","?"],[512,"S","pSer","?"],[515,"S","pSer","PDPK1 and TTBK1"],[515,"S","pSer","PDPK1 and TTBK1"],[516,"S","pSer","PDPK1 and TTBK1"],[519,"S","pSer","CK1 PDPK1 and TTBK1"],[531,"S","pSer","PKA"],[552,"S","pSer","PDPK1"],[554,"S","pSer","PHK"],[579,"S","pSer","MARK1, MARK2, MARK3, MARK4, BRSK1, BRSK2 and PHK"],[602,"S","pSer","PHK"],[610,"S","pSer","?"],[622,"S","pSer","PHK"],[641,"S","pSer","?"],[669,"S","pSer","PHK"],[673,"S","pSer","?"],[713,"S","pSer","CK1 and PDPK1"],[721,"S","pSer","CK1 and PDPK1"],[726,"S","pSer","?"],[733,"S","pSer","CaMK2 and TTBK1"],[739,"S","pSer","PDPK1 and TTBK1"]]
K_Acetyl = [[480,"K","N6-acK","alternate"],[542,"K","N6-acK","?"],[576,"K","N6-acK","alternate"],[615,"K","N6-acK","alternate"],[628,"K","N6-acK","alternate"],[634,"K","N6-acK","alternate"],[638,"K","N6-acK","alternate"],[648,"K","N6-acK","alternate"],[660,"K","N6-acK","alternate"],[664,"K","N6-acK","alternate"],[686,"K","N6-acK","alternate"],[702,"K","N6-acK","alternate"]]
N_glyc_K = [[87,"K","Gly",""],[383,"K","Gly",""],[467,"K","Gly",""],[480,"K","Gly",""],[491,"K","Gly",""],[542,"K","Gly",""],[576,"K","Gly",""],[597,"K","Gly",""],[598,"K","Gly",""],[664,"K","Gly",""],[686,"K","Gly",""]]
N_glyc_D = [[87,"D","Gly",""],[383,"D","Gly",""],[467,"D","Gly",""],[480,"D","Gly",""],[491,"D","Gly",""],[542,"D","Gly",""],[576,"D","Gly",""],[597,"D","Gly",""],[598,"D","Gly",""],[664,"D","Gly",""],[686,"D","Gly",""]]
N_glyc_isoD = [[87,"isoD","Gly",""],[383,"isoD","Gly",""],[467,"isoD","Gly",""],[480,"isoD","Gly",""],[491,"isoD","Gly",""],[542,"isoD","Gly",""],[576,"isoD","Gly",""],[597,"isoD","Gly",""],[598,"isoD","Gly",""],[664,"isoD","Gly",""],[686,"isoD","Gly",""]]
K_total = [[i,"K","N6-acK",""] for i in range(0, len(tau)) if tau[i] == "K"]
S_total = [[i,"S","pSer",""] for i in range(0, len(tau)) if tau[i] == "S"]



DDPTMS=[S_Phospo, K_Acetyl, N_glyc_K, N_glyc_D, N_glyc_isoD, K_total, S_total]
with open('data/PTMList_DANIEL.txt', 'w') as outFile:
    writer=csv.writer(outFile, delimiter='\t')
    writer.writerow(['POS','AA','PTM','OTHER', 'SOURCE', 'FREQ'])
    for item in DDPTMS:
        for ptm in item:
            ptm.extend(['DD', '-'])
            writer.writerow(ptm)

## read the ptms & ## harmonize the 2 lists
convDict={'K-Acetyl':"N6-acK",
 'K-GlyGly':'Gly',
 'K-Methyl':'K-Methyl',
 'K-Phospho':'pLys',
 'R-Methyl':'R-Methyl',
 'S-Phospho':'pSer',
 'T-Phospho':'pThr',
 'Y-Phospho':'pTyr'}
ptmFile = 'data/ptmList_Wessling_et_al.csv'
wTauPtm = set()
outFile=open('data/PTMList_DANIEL.txt', 'a')# as outFile:
writer=csv.writer(outFile, delimiter='\t')
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
                    ## check pos
                    if tauSeq[posParsed-1] == aaParsed:
                        print(row, key)
                        makeRow = [posParsed, aaParsed, convKey, '-', 'Wessling_et_al', '{}-{}'.format(ad, ct)]
                        writer.writerow(makeRow)
            else:
                posParsed = int(re.findall('[0-9]+', pos)[0])
                aaParsed = re.findall('[A-Z]+', pos)[0]
                convKey = convDict.get(key)
                if tauSeq[posParsed-1] == aaParsed:
                    print(row, key)
                    makeRow = [posParsed, aaParsed, convKey, '-', 'Wessling_et_al', '{}-{}'.format(ad, ct)]
                    writer.writerow(makeRow)
outFile.close()