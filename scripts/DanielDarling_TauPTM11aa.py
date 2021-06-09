tau = "-MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQEPESGKVVQEGFLREPGPPGLSHQLMSGMPGAPLLPEGPREATRQPSGTGPEDTEGGRHAPELLKHQLLGDLHQEGPPLKGAGGKERPGSKEEVDEDRDVDESSPQDSPPSKASPAQDGRPPQTAAREATSIPGFPAEGAIPLPVDFLSKVSTEIPASEPDGPSVGRAKGQDAPLEFTFHVEITPNVQKEQAHSEEHLGRAAFPGAPGEGPEARGPSLGEDTKEADLPEPSEKQPAAAPRGKPVSRVPQLKARMVSKSKDGTGSDDKKAKTSTRSSAKTLKNRPCLSPKHPTPGSSDPLIQPSSPAVCPEPPSSPKYVSSVTSRTGSSGAKEMKLKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
# Sites based on uniprot and select sources
S3P = [[46,"S","pSer","?"],[61,"S","pSer","?"],[214,"S","pSer","SGK1"],[396,"S","pSer","PHF-tau"],[502,"S","pSer","?"],[508,"S","pSer","?"],[512,"S","pSer","?"],[515,"S","pSer","PDPK1 and TTBK1"],[515,"S","pSer","PDPK1 and TTBK1"],[516,"S","pSer","PDPK1 and TTBK1"],[519,"S","pSer","CK1 PDPK1 and TTBK1"],[531,"S","pSer","PKA"],[552,"S","pSer","PDPK1"],[554,"S","pSer","PHK"],[579,"S","pSer","MARK1, MARK2, MARK3, MARK4, BRSK1, BRSK2 and PHK"],[602,"S","pSer","PHK"],[610,"S","pSer","?"],[622,"S","pSer","PHK"],[641,"S","pSer","?"],[669,"S","pSer","PHK"],[673,"S","pSer","?"],[713,"S","pSer","CK1 and PDPK1"],[721,"S","pSer","CK1 and PDPK1"],[726,"S","pSer","?"],[733,"S","pSer","CaMK2 and TTBK1"],[739,"S","pSer","PDPK1 and TTBK1"]]
Kac = [[480,"K","N6-acK","alternate"],[542,"K","N6-acK","?"],[576,"K","N6-acK","alternate"],[615,"K","N6-acK","alternate"],[628,"K","N6-acK","alternate"],[634,"K","N6-acK","alternate"],[638,"K","N6-acK","alternate"],[648,"K","N6-acK","alternate"],[660,"K","N6-acK","alternate"],[664,"K","N6-acK","alternate"],[686,"K","N6-acK","alternate"],[702,"K","N6-acK","alternate"]]
N_glyc_K = [[87,"K","Gly",""],[383,"K","Gly",""],[467,"K","Gly",""],[480,"K","Gly",""],[491,"K","Gly",""],[542,"K","Gly",""],[576,"K","Gly",""],[597,"K","Gly",""],[598,"K","Gly",""],[664,"K","Gly",""],[686,"K","Gly",""]]
N_glyc_D = [[87,"D","Gly",""],[383,"D","Gly",""],[467,"D","Gly",""],[480,"D","Gly",""],[491,"D","Gly",""],[542,"D","Gly",""],[576,"D","Gly",""],[597,"D","Gly",""],[598,"D","Gly",""],[664,"D","Gly",""],[686,"D","Gly",""]]
N_glyc_isoD = [[87,"isoD","Gly",""],[383,"isoD","Gly",""],[467,"isoD","Gly",""],[480,"isoD","Gly",""],[491,"isoD","Gly",""],[542,"isoD","Gly",""],[576,"isoD","Gly",""],[597,"isoD","Gly",""],[598,"isoD","Gly",""],[664,"isoD","Gly",""],[686,"isoD","Gly",""]]

# All possible acetylation and phosphorylation sites at S and K
K_total = [[i,"K","N6-acK",""] for i in range(len(tau)-9) if tau[i] == "K"]
S_total = [[i,"S","pSer",""] for i in range(len(tau)-9) if tau[i] == "S"]

# function outputs a list of both the modified version, and unmodified version of the sequence, including position,
# modification type and modifier if available

def PTM_11aa_seqgen(seq,ptmlist):
    outputseqlist = []
    for i in ptmlist:
        # each of the 9 possible core positions are appended to the output sequence list for each ptm in the
        # inputed list
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 1", "11aa_seq: " + seq[(int(i[0])-1):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+10],"11aa_seq_modified: "+seq[(int(i[0])-1):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+10], "Modifier: " + i[3]])
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 2", "11aa_seq: " + seq[(int(i[0])-2):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+9],"11aa_seq_modified: "+seq[(int(i[0])-2):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+9], "Modifier: " + i[3]])
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 3", "11aa_seq: " + seq[(int(i[0])-3):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+8],"11aa_seq_modified: "+seq[(int(i[0])-3):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+8], "Modifier: " + i[3]])
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 4", "11aa_seq: " + seq[(int(i[0])-4):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+7],"11aa_seq_modified: "+seq[(int(i[0])-4):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+7], "Modifier: " + i[3]])
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 5", "11aa_seq: " + seq[(int(i[0])-5):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+6],"11aa_seq_modified: "+seq[(int(i[0])-5):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+6], "Modifier: " + i[3]])
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 6", "11aa_seq: " + seq[(int(i[0])-6):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+5],"11aa_seq_modified: "+seq[(int(i[0])-6):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+5], "Modifier: " + i[3]])
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 7", "11aa_seq: " + seq[(int(i[0])-7):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+4],"11aa_seq_modified: "+seq[(int(i[0])-7):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+4], "Modifier: " + i[3]])
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 8", "11aa_seq: " + seq[(int(i[0])-8):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+3],"11aa_seq_modified: "+seq[(int(i[0])-8):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+3], "Modifier: " + i[3]])
        outputseqlist.append(["PTM_pos: " + str(i[0]), "Modification_type: " + i[2], "Core_Placement: 9", "11aa_seq: " + seq[(int(i[0])-9):int(i[0])]+i[1]+seq[(int(i[0])+1):int(i[0])+2],"11aa_seq_modified: "+seq[(int(i[0])-9):int(i[0])]+i[1]+"("+i[2]+")"+seq[(int(i[0])+1):int(i[0])+2], "Modifier: " + i[3]])
    return outputseqlist




sortedptms_known = sorted(N_glyc_D+N_glyc_K+N_glyc_isoD+S3P+Kac, key= lambda x:x[0])


knownseqoutput = [[i[3][10:],str(i[0][9:])+"_"+i[1][19:]+"_cp"+i[2][16:]] for i in PTM_11aa_seqgen(tau,sortedptms_known)]
knownseqoutput += [[i[4][19:],str(i[0][9:])+"_modified_"+i[1][19:]+"_cp"+i[2][16:]] for i in PTM_11aa_seqgen(tau,sortedptms_known)]
setseqoutput = []
for i in knownseqoutput:
    test = True
    for j in setseqoutput:
        if j[0] == i[0]:
            test = False
    if test == False:
        continue
    else:
        setseqoutput.append(i)


print(len(setseqoutput))

f = open("PTM_sequences_taunames.txt", "w")
for i in setseqoutput:
    if setseqoutput.index(i) == len(setseqoutput)-1:
        f.write(i[1])
    else:
        f.write(i[1]+", ")


f.close()
print()
