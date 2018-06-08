# This program calculates a maximum tankyrase-targeting score (TTS) for each protein.

# 8 × 20 PSSM generated in Guettler et al.
fP = [-0.5000, 0.0546, 0.1025, 0.1852, -0.5000, -0.5000, -0.5000, 0.0468]
fG = [-0.5000, 0.0505, 0.0256, 0.6792, -0.5000, 1.0000, 0.0318, 0.0318]
fA = [-0.5000, 0.0426, 0.1132, 0.0896, -0.5000, -0.5000, 0.0721, 0.0674]
fV = [-0.5000, 0.0297, 0.0365, -0.5000, 0.0690, -0.5000, 0.0305, 0.0374]
fL = [-0.5000, 0.0718, 0.0259, -0.5000, -0.5000, -0.5000, 0.0215, 0.0255]
fI = [-0.5000, 0.0307, 0.0197, -0.5000, 0.0380, -0.5000, 0.0400, 0.0291]
fM = [-0.5000, 0.0606, 0.0456, -0.5000, -0.5000, -0.5000, 0.0522, 0.0337]
fC = [-0.5000, 0.0384, 0.0394, 0.0460, 0.0361, -0.5000, 0.1082, 0.0481]
fS = [-0.5000, 0.0635, 0.0688, -0.5000, -0.5000, -0.5000, 0.0597, 0.0481]
fT = [-0.5000, 0.0738, 0.0385, -0.5000, -0.5000, -0.5000, 0.0549, 0.0301]
fR = [1.0000, 0.0222, 0.0560, -0.5000, -0.5000, -0.5000, 0.0301, 0.0102]
fK = [-0.5000, 0.0195, 0.0431, -0.5000, -0.5000, -0.5000, 0.0308, 0.0090]
fH = [-0.5000, 0.0337, 0.0265, -0.5000, -0.5000, -0.5000, 0.0484, 0.0312]
fD = [-0.5000, 0.0738, 0.1108, -0.5000, 0.6246, -0.5000, 0.0816, 0.1872]
fE = [-0.5000, 0.1137, 0.1045, -0.5000, 0.1249, -0.5000, 0.1178, 0.2106]
fN = [-0.5000, 0.0255, 0.0300, -0.5000, -0.5000, -0.5000, 0.0530, 0.0455]
fQ = [-0.5000, 0.0479, 0.0584, -0.5000, 0.0592, -0.5000, 0.0763, 0.0259]
fW = [-0.5000, 0.0455, 0.0335, -0.5000, -0.5000, -0.5000, 0.0295, 0.0183]
fF = [-0.5000, 0.0440, -0.5000, -0.5000, -0.5000, -0.5000, 0.0341, 0.0272]
fY = [-0.5000, 0.0581, 0.0216, -0.5000, 0.0482, -0.5000, 0.0274, 0.0366]
sAminoAcidDic = {"P": fP, "G": fG, "A": fA, "V": fV, "L": fL, "I": fI, "M": fM, "C": fC, "S": fS, "T": fT, "R": fR, "K": fK, "H": fH, "D": fD, "E": fE, "N": fN, "Q": fQ, "W": fW, "F": fF, "Y": fY}

# CalTTS calculates a maximum TTS for each protein.
def CalTTS(sProteinSeq):

    fMaxTTS = -0.770
    nMaxPos = 0

    if len(sProteinSeq) <8:

        print("This protein is smaller than 8 aa.")

    else:

        for nPos in range(len(sProteinSeq)-7):

            fTTS = 0.000

            if ("U" not in sProteinSeq[nPos:nPos+8]) and ("X" not in sProteinSeq[nPos:nPos+8]) and ("B" not in sProteinSeq[nPos:nPos+8]):

                for nNum in range(8):

                    fTTS += sAminoAcidDic[sProteinSeq[nPos + nNum]][nNum]

                if fTTS > fMaxTTS:

                    fMaxTTS = max(fMaxTTS, fTTS)

                    nMaxPos = nPos

        if fMaxTTS == -0.770:

            print("fMaxTTS is -0.770.")

        else:

            return "%.3f"%(fMaxTTS/3.86), nMaxPos + 1, sProteinSeq[nMaxPos:nMaxPos+8]

def main():

    # MsProteins.fasta should contain all protein sequences of your organism of interest in the Uniprot database.
    InFile = open("MsProteins.fasta", "r")
    
    sMsProteinEntryList = []
    sMsProteinNameDic = {}
    sMsProteinSeqDic = {}
    nProteinNum = 0

    for sReadLine in InFile.readlines():

        if sReadLine[0] == ">":

            sMsProteinEntryList.append(sReadLine[4:10])
            sString = sReadLine[11:]
            sProteinName = ""

            for cChar in sString:

                if cChar == " ":

                    break

                else:

                    sProteinName += cChar

            sMsProteinNameDic[sMsProteinEntryList[nProteinNum]] = sProteinName

            nProteinNum += 1

        else:

            try:

                sMsProteinSeqDic[sMsProteinEntryList[nProteinNum - 1]] += sReadLine.replace("\n", "")

            except:

                sMsProteinSeqDic[sMsProteinEntryList[nProteinNum - 1]] = sReadLine.replace("\n", "")

    InFile.close()

    # IdentifiedProteins.txt should contain Uniprot Accession Numbers of your identified tankyrase-binding proteins.
    InFile = open("IdentifiedProteins.txt", "r")

    for sReadLine in InFile.readlines():

        sEntry = sReadLine[0:6]

        fMaxTTS, fMaxPos, sMaxSeq = CalTTS(sMsProteinSeqDic[sEntry])

        print(sEntry, sMsProteinNameDic[sEntry], fMaxTTS, fMaxPos, sMaxSeq)

        for nPos in range(8):

            print(sAminoAcidDic[sMaxSeq[nPos]][nPos])

        if "U" in sMsProteinSeqDic[sEntry]:

            print("U is in this protein seq.\n")

        if "X" in sMsProteinSeqDic[sEntry]:

            print("X is in this protein seq.\n")

        if "B" in sMsProteinSeqDic[sEntry]:

            print("B is in this protein seq.\n")

main()
