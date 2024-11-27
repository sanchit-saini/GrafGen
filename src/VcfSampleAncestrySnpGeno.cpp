#include "VcfSampleAncestrySnpGeno.h"

VcfSampleAncestrySnpGeno::VcfSampleAncestrySnpGeno(string file, AncestrySnps *aSnps, int singleChr, int prophFlag)
{
    onechr   = singleChr;
    prophage = prophFlag;

    ancSnps = aSnps;
    vcfFile = file;

    vcfSamples = {};
    //vcfAncSnpGtRefs = {};
    vcfAncSnpGtAlts = {};
    vcfAncSnpChrs = {};
    vcfAncSnpPoss = {};
    //vcfAncSnpSnps = {};
    vcfAncSnpRefs = {};
    vcfAncSnpAlts = {};
    //vcfRsIdAncSnpIds = {};
    vcfGb37AncSnpIds = {};
    //vcfGb38AncSnpIds = {};
    vcfAncSnpCodedGenos = {};
    vcfAncSnpIds = {};

    totAncSnps = ancSnps->GetNumAncestrySnps();
    numSamples = 0;
    totVcfSnps = 0;
    putativeAncSnps = 0;
    numRsIdAncSnps = 0;
    numGb37AncSnps = 0;
    numGb38AncSnps = 0;
    numVcfAncSnps = 0;
}

VcfSampleAncestrySnpGeno::~VcfSampleAncestrySnpGeno()
{
    vcfSamples.clear();
    vcfAncSnpChrs.clear();
    vcfAncSnpPoss.clear();
    //vcfAncSnpSnps.clear();
    vcfAncSnpRefs.clear();
    vcfAncSnpAlts.clear();
    //vcfRsIdAncSnpIds.clear();
    vcfGb37AncSnpIds.clear();
    //vcfGb38AncSnpIds.clear();
    vcfAncSnpIds.clear();

    DeleteAncSnpCodedGenos();
    DeleteAncSnpGtValues();
}

void VcfSampleAncestrySnpGeno::DeleteAncSnpCodedGenos()
{
    long n = (long) vcfAncSnpCodedGenos.size();
    for (int i = 0; i < n; i++) {
        if (vcfAncSnpCodedGenos[i]) delete vcfAncSnpCodedGenos[i];
    }
    vcfAncSnpCodedGenos.clear();
}

void VcfSampleAncestrySnpGeno::DeleteAncSnpGtValues()
{
    long n = (long) vcfAncSnpGtAlts.size();
    for (int i = 0; i < n; i++) {
        //if (!vcfAncSnpGtRefs[i].empty()) vcfAncSnpGtRefs[i].clear();
        if (!vcfAncSnpGtAlts[i].empty()) vcfAncSnpGtAlts[i].clear();
    }
    //vcfAncSnpGtRefs.clear();
    vcfAncSnpGtAlts.clear();
}

bool VcfSampleAncestrySnpGeno::ReadDataFromFile(int print)
{
    // Allocate memory
    char *colValue = new char[WORDLEN+1];
    char *buffer   = new char[BUFFERLEN+1];

    if (print) Rprintf("Reading data from file %s\n", vcfFile.c_str());
    gzFile file = gzopen (vcfFile.c_str(), "r");

    int lineNo = 0;
    int numVcfSnps = 0;

    vector<int> ancSnpIds;
    vector<string> chromosomes;
    vector<string> positions;
    vector<string> ids;
    vector<string> refs;
    vector<string> alts;
    vector<vector<int>> snpGenoVals;

    colValue[0] = 0;
    int valPos = 0;
    bool fileDone = false;
    bool hasHeadRow = false;
    int buffNo = 0;
    int vcfColNo = 0;

    string gtyStr, chrStr, posStr, snpStr, refStr, altStr;

    vector<string> snpGts;

    while (!fileDone) {
        int err;
        int bytesRead;
        bytesRead = gzread(file, buffer, BUFFERLEN);
        buffer[bytesRead] = '\0';

        if (bytesRead < BUFFERLEN - 1) {
            if (gzeof(file)) {
                fileDone = true;
            }
            else {
                const char * errorString;
                errorString = gzerror(file, & err);
                if (err) {
                    if (print) Rprintf("Error: %s.\n", errorString);
                    Rprintf("ERROR reading binary file");
                }
            }
        }

        int buffPos = 0;
        while (buffPos < bytesRead) {

            //bool isNewLine = false;
            if (buffPos >= BUFFERLEN) Rprintf("ERROR BUFF 0");
            if (buffer[buffPos] == '\t' || buffer[buffPos] == '\n') {

                //if (buffer[buffPos] == '\n') isNewLine = true;
                if (valPos >= WORDLEN) Rprintf("ERROR 1");
                colValue[valPos] = 0;
                if      (vcfColNo == 0) {
                    chrStr = string(colValue);
                    if (!onechr) chromosomes.push_back(chrStr);
                }
                else if (vcfColNo == 1) {
                    posStr = string(colValue);
                    positions.push_back(posStr);
                }
                else if (vcfColNo == 2) {
                    //snpStr = string(colValue);
                    //ids.push_back(snpStr);
                }
                else if (vcfColNo == 3) {
                    refStr = string(colValue);
                    refs.push_back(refStr);
                }
                else if (vcfColNo == 4) {
                    altStr = string(colValue);
                    alts.push_back(altStr);
                }
                else if (vcfColNo == 8) {
                    gtyStr = string(colValue);
                }
                else if (vcfColNo > 8)  {
                    string gtStr = string(colValue);
                    snpGts.push_back(gtStr);
                }
                valPos = 0;
                colValue[0] = 0;
                vcfColNo++;
            }
            else {
                if (buffPos >= BUFFERLEN) Rprintf("ERROR BUFF 1");
                if (buffer[buffPos] != '\r' && valPos < WORDLEN-1) {
                    colValue[valPos] = buffer[buffPos];
                    valPos++;
                }
                else {
                    if (valPos >= WORDLEN) Rprintf("ERROR 2");
                    colValue[valPos] = 0;
                }
            }
            if (buffPos >= BUFFERLEN) Rprintf("ERROR BUFF 2");
            if (buffer[buffPos] == '\n') {
                vcfColNo = 0;
                lineNo++;

                bool isGt = false;
                if (gtyStr.size() > 1 && gtyStr[0] == 'G' && gtyStr[1] == 'T') isGt = true;
                int numCols = snpGts.size();

                if (strcmp(chrStr.c_str(), "#CHROM") == 0) {
                    if (numSamples > 0) {
                        assert(numCols == numSamples);
                    }
                    else {
                        if (numCols > 0) {
                            for (int smpNo = 0; smpNo < numCols; smpNo++) vcfSamples.push_back(snpGts[smpNo]);
                            numSamples = numCols;
                            if (print) Rprintf("Vcf file has %d samples\n", numSamples);
                        }
                        else {
                            if (print) Rprintf("ERROR: vcf file %s doesn't include samples!\n",
                                    vcfFile.c_str());
                            Rprintf("ERROR with VCF file");
                            fileDone = true;
                            return false;
                        }
                    }

                    hasHeadRow = true;
                    snpGts.clear();
                }
                else if (chrStr[0] && chrStr[0] != '#') {
                    if (!hasHeadRow) {
                        Rprintf("ERROR: didn't find #CHROM row in vcf file");
                        return false;
                    }

                    if (numCols != numSamples) {
                      if (print) Rprintf("ERROR at line # %d \n", lineNo);
                      Rprintf("ERROR with VCF file");
                      return false;
                    }
                    int chr, gb37SnpId;

                    //if (!onechr) chr = GetChromosomeFromString(chrStr.c_str());
                    if (!onechr) {
                      chr = ancSnps->GetChrNumFromChrStr((char *) chrStr.c_str());
                      // if chr < 0, then gb37SnpId will be <= -1 below
                    }
                    //int rsNum = GetRsNumFromString(snpStr.c_str());
                    int pos = 0;
                    try { pos = stoi(posStr); }
                    catch (exception &err) { pos = 0; }

                    /*
                    int rsSnpId = ancSnps->FindSnpIdGivenRs(rsNum);
                    int gb37SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos, 37);
                    int gb38SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos, 38);
                    */
                    if (onechr) {
                      gb37SnpId = ancSnps->FindSnpIdGivenPos(pos);
                    } else {
                      gb37SnpId = ancSnps->FindSnpIdGivenChrPos(chr, pos);
                    }
                    if (isGt && (gb37SnpId > -1)) {
                        putativeAncSnps++;
                        //if (rsSnpId > -1)   numRsIdAncSnps++;
                        if (gb37SnpId > -1) numGb37AncSnps++;
                        //if (gb38SnpId > -1) numGb38AncSnps++;

                        if (!onechr) vcfAncSnpChrs.push_back(chr);
                        vcfAncSnpPoss.push_back(pos);
                        //vcfAncSnpSnps.push_back(snpStr);
                        vcfAncSnpRefs.push_back(refStr);
                        vcfAncSnpAlts.push_back(altStr);

                        //vcfRsIdAncSnpIds.push_back(rsSnpId);
                        vcfGb37AncSnpIds.push_back(gb37SnpId);
                        //vcfGb38AncSnpIds.push_back(gb38SnpId);

                        //vector<char> gtRefs;
                        vector<char> gtAlts;
                        for (int i = 0; i < numCols; i++) {
                            string gtStr = snpGts[i];
                            int gtLen = gtStr.size();
                            char gAlt = 0;
                            /*
                            if (gtLen > 2 && (gtStr[1] == '|' || gtStr[1] == '/')) {
                                gRef = gtStr[0];
                                gAlt = gtStr[2];
                            }
                            */
                            if (gtLen == 1) gAlt = gtStr[0];

                            //gtRefs.push_back(gRef);
                            gtAlts.push_back(gAlt);
                        }

                        //vcfAncSnpGtRefs.push_back(gtRefs);
                        vcfAncSnpGtAlts.push_back(gtAlts);
                    }

                    numVcfSnps++;
                    snpGts.clear();
                }

                if (print && (lineNo % 1000000 == 0)) {
                  Rprintf("Checked %d lines. Found %d lines with ancestry SNPs\n",
                           lineNo, putativeAncSnps);
                }
            }

            buffPos++;
        }

        buffNo++;
    }

    lineNo++;
    delete[] buffer;
    delete[] colValue;

    if (print) Rprintf("Done. Checked %d lines. Found %d lines with ancestry SNPs\n",
              lineNo, putativeAncSnps);
    gzclose (file);
    totVcfSnps += numVcfSnps;

    return true;
}

void VcfSampleAncestrySnpGeno::RecodeSnpGenotypes()
{
    /*
    ancSnpType = AncestrySnpType::RSID;
    int maxVcfAncSnps = numRsIdAncSnps;

    if (numGb37AncSnps > maxVcfAncSnps) {
        ancSnpType = AncestrySnpType::GB37;
        maxVcfAncSnps = numGb37AncSnps;
    }

    if (numGb38AncSnps > maxVcfAncSnps) {
        ancSnpType = AncestrySnpType::GB38;
        maxVcfAncSnps = numGb38AncSnps;
    }
    */
    ancSnpType = AncestrySnpType::GB37;
    //int maxVcfAncSnps = numGb37AncSnps;


    int saveSnpNo = 0; // Putative SNPs saved after vcf file was read
    int ancSnpNo = 0;  // Final list of ancestry SNPs to be used for ancestry inference

    for (saveSnpNo = 0; saveSnpNo < putativeAncSnps; saveSnpNo++) {
        /*
        int ancSnpId = -1;
        if      (ancSnpType == AncestrySnpType::RSID) ancSnpId = vcfRsIdAncSnpIds[saveSnpNo];
        else if (ancSnpType == AncestrySnpType::GB37) ancSnpId = vcfGb37AncSnpIds[saveSnpNo];
        else if (ancSnpType == AncestrySnpType::GB38) ancSnpId = vcfGb38AncSnpIds[saveSnpNo];
        */
        int ancSnpId = vcfGb37AncSnpIds[saveSnpNo];

        if (ancSnpId > -1) {
            char eRef = ancSnps->snps[ancSnpId].ref;
            char eAlt = ancSnps->snps[ancSnpId].alt;

            string vcfRef = vcfAncSnpRefs[saveSnpNo];
            string vcfAlt = vcfAncSnpAlts[saveSnpNo];

            int expRefIdx = -1;
            int expAltIdx = -1;

            if (prophage) {
              CompareAncestrySnpAlleles_prophage(vcfRef, vcfAlt, eRef, eAlt, &expRefIdx, &expAltIdx);
            } else {
              CompareAncestrySnpAlleles(vcfRef, vcfAlt, eRef, eAlt, &expRefIdx, &expAltIdx);
            }

            if (expRefIdx > -1 && expAltIdx > -1) {
                char* smpGenos = new char[numSamples];
                //vector<char> gtRefVal = vcfAncSnpGtRefs[saveSnpNo];
                vector<char> gtAltVal = vcfAncSnpGtAlts[saveSnpNo];

                // Genotypes for the samples consist of (0=ref) or the alternate allele id number (1, 2, ...)
                for (int smpNo = 0; smpNo < numSamples; smpNo++) {
                    //int refGval = gtRefVal[smpNo] - '0';
                    int altGval = gtAltVal[smpNo] - '0';
                    char geno = RecodeGenotypeGivenIntegers(expRefIdx, expAltIdx, altGval);
                    smpGenos[smpNo] = geno; // geno will be 0 or 1 (or 3 if missing)
                }
                vcfAncSnpIds.push_back(ancSnpId);
                vcfAncSnpCodedGenos.push_back(smpGenos);
                ancSnpNo++;
            }
        }
    }

    DeleteAncSnpGtValues();
}

void VcfSampleAncestrySnpGeno::CompareAncestrySnpAlleles(const string refStr, const string altsStr,
const char eRef, const char eAlt, int* expRefIdx, int* expAltIdx)
{
    // Let index of the ref alleles in vcf be 0, and indices of alts be 1, 2, 3, ...
	// Check them to find which one is the expected ref and and which is alt for the Ancestry SNP
    *expRefIdx = -1;
    *expAltIdx = -1;

    int refLen = refStr.size();
    if (refLen == 1) {
        char ref = refStr[0];
        char fRef = FlipAllele(ref);  // Flipped ref allele

        if (ref == eRef || fRef == eRef) {
            *expRefIdx = 0; 
        }
        else if (ref == eAlt || fRef == eAlt) {
            *expAltIdx = 0;
        }

        // If ref seems to be flipped, then all the alts are flipped
        bool flip = false;
        if (fRef == eRef || fRef == eAlt) {
            flip = true;
        }

        if (*expRefIdx > -1 || *expAltIdx > -1) {
            vector<string> altWords = SplitString(altsStr, ",");
            int numAlts = altWords.size();

            for (int altNo = 0; altNo < numAlts; altNo++) {
                string altWord = altWords[altNo];
                int alleleIdx = altNo + 1;  // allele index starts from 0 = ref, then 1 = first alt, ...
	      if (altWord.size() == 1) {
                    char alt = altWord[0];
                    if (flip)  alt = FlipAllele(alt);

                    if      (alt == eRef && *expRefIdx < 0) *expRefIdx = alleleIdx;
                    else if (alt == eAlt && *expAltIdx < 0) *expAltIdx = alleleIdx;
                }
            }
        } // end if (*expRefIdx ...)
    }
}

void VcfSampleAncestrySnpGeno::CompareAncestrySnpAlleles_prophage(const string refStr, const string altsStr,
const char eRef, const char eAlt, int* expRefIdx, int* expAltIdx)
{
    // The REF allele must match - this is required for prophage data.
    // For the ALT allele, an asterisk will be treated as a missing allele. 
    // Any other allele is the variant allele.
    // Currently, eAlt is not used
    *expRefIdx = -1;
    *expAltIdx = -1;

    int refLen = refStr.size();
    if (refLen == 1) {
        char ref = refStr[0];

        if (ref == eRef) {
          *expRefIdx = 0; 
        } else {
          return;  // REFs must match !!!
        }
        
        // Check ALT
        vector<string> altWords = SplitString(altsStr, ",");
        int numAlts = altWords.size();

        for (int altNo = 0; altNo < numAlts; altNo++) {
          string altWord = altWords[altNo];
          int alleleIdx = altNo + 1;  // allele index starts from 0 = ref, then 1 = first alt, ...
	  if (altWord.size() == 1) {
            char alt = altWord[0];
            if ((alt == 'A') || (alt == 'C') || (alt == 'G') || (alt == 'T')) {
              *expAltIdx = alleleIdx;
              return;
            }
          }
        }
    }
}


int VcfSampleAncestrySnpGeno::RecodeGenotypeGivenString(const int expRefIdx, const int expAltIdx, const string genoStr)
{
    int genoInt = 3;

    // genoStr should just be a single integer 0=ref, or k for which alternate allele 
    if (genoStr.size() == 1 && expRefIdx > -1 && expAltIdx > -1) {
      int gNum = genoStr[0] - '0';
      genoInt = RecodeGenotypeGivenIntegers(expRefIdx, expAltIdx, gNum);
    }

    //cout << "exp " << expRefIdx << expAltIdx << " geno " << genoStr << " => " << genoInt << "\n";

    return genoInt;
}

int VcfSampleAncestrySnpGeno::RecodeGenotypeGivenIntegers(const int expRefIdx, const int expAltIdx, const int gNum)
{
    int genoInt = 3; // number of alt alleles, valid counts are 0, 1

    // Genotype is valid only if both alleles are valid, i.e., same as one of the expected alleles
    int numValidGenos = 0;
    int numAlts = 0;

    if (gNum == expRefIdx || gNum == expAltIdx) {
        numValidGenos++;
        if (gNum == expAltIdx) numAlts++;
    }

    if (numValidGenos == 1) genoInt = numAlts;

    return genoInt;
}

void VcfSampleAncestrySnpGeno::ShowSummary()
{
    Rprintf("Total %d ancestry SNPs used by GrafPop\n", totAncSnps);
    Rprintf("Number of samples found in the vcf file: %d\n", numSamples);
    Rprintf("Total %d ancestry SNPs found from %d SNPs\n", putativeAncSnps, totVcfSnps);

    //cout << "\n#RSID Ancs: " << numRsIdAncSnps << "\n"
    //<< "#GB37 Ancs: " << numGb37AncSnps << "\n"
    //<< "#GB38 Ancs: " << numGb38AncSnps << "\n";
}
