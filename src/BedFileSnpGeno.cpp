#include "BedFileSnpGeno.h"

BedFileSnpGeno::BedFileSnpGeno(string bFile, AncestrySnps *aSnps, BimFileAncestrySnps *bSnps, FamFileSamples *fSmps)
{
    bedFile = bFile;
    ancSnps = aSnps;
    bimSnps = bSnps;
    famSmps = fSmps;

    unsigned long base = 1;
    for (int i = 0; i < 64; i++) {
        baseNums[i] = base;
        base = base << 1;
    }

    numAncSnps = ancSnps->GetNumAncestrySnps();
    numBimSnps = bimSnps->GetNumBimSnps();
    numBimAncSnps = bimSnps->GetNumBimAncestrySnps();
    numSamples = famSmps->GetNumFamSamples();

    ancSnpSmpGenos = {};
    ancSnpSnpIds = {};

    vtxExpGd0 = new SampleGenoDist(&aSnps->vtxPopExpGds[0], &aSnps->vtxPopExpGds[1],
    &aSnps->vtxPopExpGds[2], &aSnps->vtxPopExpGds[0]);
    vtxExpGd0->TransformAllDists();
    vtxExpGd0->CalculateBaryCenters();
}

BedFileSnpGeno::~BedFileSnpGeno()
{
    long i, n;
    n = (long) ancSnpSmpGenos.size();
    for (i = 0; i < n; i++) {
        delete ancSnpSmpGenos[i];
    }
    ancSnpSmpGenos.clear();
    ancSnpSnpIds.clear();
}

char BedFileSnpGeno::GetCompAllele(char a)
{
    char c = '0';
    if      (a == 'A') c = 'T';
    else if (a == 'T') c = 'A';
    else if (a == 'G') c = 'C';
    else if (a == 'C') c = 'G';

    return c;
}

char* BedFileSnpGeno::RecodeBedSnpGeno(char *snpBedGenos, int numBytes, bool swap)
{
    char *snpGenos = new char[numSamples]; // char only takes one byte
    for (int i = 0; i < numSamples; i++) snpGenos[i] = 3;

    int smpNo = 0;
    int byteNo = 0;

    for (byteNo = 0; byteNo < numBytes; byteNo++) {
        char genoByte = snpBedGenos[byteNo];
        //int val = int(genoByte);

        for (int byteSmpNo = 0; byteSmpNo < 4; byteSmpNo++) {
            int bit1Pos = byteSmpNo * 2;
            int bit2Pos = bit1Pos + 1;;

            int bit1 = genoByte & baseNums[bit1Pos];
            int bit2 = genoByte & baseNums[bit2Pos];

            int intGeno = 3;
            /* Original code
            if      ( bit1 &&  bit2) intGeno = 2;
            else if (!bit1 &&  bit2) intGeno = 1;
            else if (!bit1 && !bit2) intGeno = 0;
            if (swap) {
                if      (intGeno == 0) intGeno = 2;
                else if (intGeno == 2) intGeno = 0;
            }
            */
            if      ( bit1 &&  bit2) intGeno = 1;
            else if (!bit1 && !bit2) intGeno = 0;
            if (swap) {
                if      (intGeno == 0) intGeno = 1;
                else if (intGeno == 1) intGeno = 0;
            }

            if (smpNo < numSamples) snpGenos[smpNo] = intGeno;
            smpNo++;
        }
    }

    return snpGenos;
}

bool BedFileSnpGeno::ReadGenotypesFromBedFile(int print)
{
    bool hasErr = false;

    char header[2];
    char mode[1];

    ifstream bedFilePtr (bedFile, ios::in | ios::binary);

    long snpNumBytes = (numSamples - 1) / 4 + 1;
    long expFileLen = snpNumBytes * numBimSnps + 3;

    bedFilePtr.seekg (0, bedFilePtr.end);
    long fileLen = bedFilePtr.tellg();
    bedFilePtr.seekg(0, bedFilePtr.beg);

    bedFilePtr.read (header, 2);
    bedFilePtr.read (mode, 1);

    if (header[0] != 108 || header[1] != 27) {
      if (print) Rprintf("ERROR: File is not a valid PLINK bed file!\n");
      hasErr = true;
      Rprintf("ERROR with input genotype file");
    }
    else if (mode[0] != 1) {
      if (print) Rprintf("ERROR: File is not in SNP mode!\n");
      hasErr = true;
      Rprintf("ERROR with input genotype file");
    }

    if (fileLen != expFileLen) {
      if (print) {
        Rprintf("ERROR: Number of genotypes in bed file doesn't match fam and bim File!\n");
        Rprintf("Fam file has %d samples.\n", numSamples); 
        Rprintf("Bim file has %d SNPs.\n", numBimSnps);
        Rprintf("Each SNP should have %ld bytes\n", snpNumBytes);
        Rprintf("Expected total %ld bytes.\n", expFileLen);
        Rprintf("Bed file has %ld bytes.\n", fileLen);
      }
      hasErr = true;
      Rprintf("ERROR with input genotype file");
    }

    if (hasErr) return hasErr;
    if (print) Rprintf("Reading genotypes from %s\n", bedFile.c_str());

    //char buff[snpNumBytes]; // Reusable memory to keep the genotypes
    char *buff = new char[snpNumBytes];
    int bimAncSnpNo = 0;

    for (int i = 0; i < numBimSnps; i++) {
        bedFilePtr.read (buff, snpNumBytes);
        int ancSnpId = bimSnps->GetAncSnpIdGivenBimSnpPos(i);
        int match = bimSnps->GetAlleleMatchGivenBimSnpPos(i);
        bool swap = match ==  2 || match == -2 ? true : false;

        if (ancSnpId >= 0) {
            char* snpGenoStr = new char[snpNumBytes];
            for (int j = 0; j < snpNumBytes; j++) snpGenoStr[j] = buff[j];
            //ASSERT(bimAncSnpNo < numAncSnps, "bim ancestry SNP ID " << bimAncSnpNo << " not less than " << numAncSnps << "\n");

            char *snpSmpGeno = RecodeBedSnpGeno(snpGenoStr, snpNumBytes, swap);

            ancSnpSmpGenos.push_back(snpSmpGeno);
            ancSnpSnpIds.push_back(ancSnpId);

            bimAncSnpNo++;
        }
    }
    delete[] buff;
    bedFilePtr.close();
    numBimAncSnps = bimAncSnpNo;

    if (print) {
      Rprintf("Read genotypes of %d Ancestry SNPs from total %d SNPs.\n",
              bimAncSnpNo, numBimSnps);
      Rprintf("Bed file has genotypes of %d SNPs.\n", numBimSnps); 
      Rprintf("Read genotypes of %d ancestry SNPs for %d samples.\n",
              numBimAncSnps, numSamples);
    }

    return 0;
}


int BedFileSnpGeno::GetSnpGenoInt(bool b1, bool b2)
{
    int g = 3;

    if (!b1 && !b2) {
        g = 0;
    }
    else if (b1 && b2) {
        g = 1;
    }

    return g;
}

void BedFileSnpGeno::ShowSummary()
{
    Rprintf("\n");
    Rprintf("Total %d samples\n", numSamples);
    Rprintf("Total %d Ancestry SNPs\n", numAncSnps);
    Rprintf("Total %d Ancestry SNPs in bim file\n\n", numBimAncSnps);
}
