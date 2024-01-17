#include "AncestrySnps.h"

AncestrySnp::AncestrySnp(int id, int g37, char a1, char a2, float* refPops, float *vtxPops)
{
    snpId = id;
    posG37 = g37;
    ref = a1;
    alt = a2;

    for (int i = 0; i < numRefPops; i++) refPopAfs[i] = refPops[i];
    for (int i = 0; i < numVtxPops; i++) vtxPopAfs[i] = vtxPops[i];
}

AncestrySnps::AncestrySnps()
{
    // see order in ReadAncestrySnpsFromFile

    snps = {};
}

AncestrySnps::~AncestrySnps()
{
    snps.clear();
    pos37ToAncSnpId.clear();
}

int AncestrySnps::ReadAncestrySnpsFromFile(string ancSnpFile, int print)
{
    double popExpPfSums[numRefPops];
    double popExpPaSums[numRefPops];
    double popExpPeSums[numRefPops];

    for (int popId = 0; popId < numRefPops; popId++) {
        popExpPeSums[popId] = 0;
        popExpPfSums[popId] = 0;
        popExpPaSums[popId] = 0;
    }

    int lineLen = 5000;
    char snpLine[5000];

    FILE *ifp = fopen(ancSnpFile.c_str(), "r");
    if (!ifp) error("ERROR: Couldn't open genotype file");

    int lineNo = 0;
    bool fileIsValid = true;

    int numSnps = 0;
    int pos;
    char a1, a2;
    float vtEur, vtAfr, vtEas;
    float rfAfr, rfAfrAm, rfEurAm, rfAkl, rfEAsia;
    float rfEurasia, rfSEur, rfNEur, rfZAF;

    while (fgets(snpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        if (lineNo == 0) {
            if (snpLine[0] != 'p' || snpLine[1] != 'o' || snpLine[2] != 's') {
                fileIsValid = false;
            }
            if (!fileIsValid) return(1);
        }
        else {
            float* refPopAfs = new float[numRefPops];
            float* vtxPopAfs = new float[numVtxPops];

            sscanf(snpLine, "%d %c %c %f %f %f %f %f %f %f %f %f %f %f %f",
            &pos, &a1, &a2, &vtAfr, &vtEas, &vtEur, 
            &rfAfr, &rfAfrAm, &rfEurAm, &rfAkl, &rfEAsia, 
            &rfEurasia, &rfSEur, &rfNEur, &rfZAF);

            refPopAfs[0] = rfNEur;
            refPopAfs[1] = rfAfr;
            refPopAfs[2] = rfEAsia;
            refPopAfs[3] = rfAfrAm; 
            refPopAfs[4] = rfEurAm; 
            refPopAfs[5] = rfSEur;
            refPopAfs[6] = rfEurasia;
            refPopAfs[7] = rfAkl;
            refPopAfs[8] = rfZAF;

            vtxPopAfs[0] = vtEur;
            vtxPopAfs[1] = vtAfr;
            vtxPopAfs[2] = vtEas;

            /*
            int rc = parseSnpLine(snpLine, nrefPop, &pos, &a1, &a2, 
                                  vtxPopAfs, refPopAfs);
            if (rc) {
              Rprintf("ERROR on line %d of ancestry SNP file", lineNo);
              continue;
            }
            */

            AncestrySnp ancSnp(numSnps, pos, a1, a2, refPopAfs, vtxPopAfs);

            snps.push_back(ancSnp);

            pos37ToAncSnpId[pos] = numSnps;

            double pev = refPopAfs[0];
            double pfv = refPopAfs[1];
            double pav = refPopAfs[2];

            double qev = 1 - pev;
            double qfv = 1 - pfv;
            double qav = 1 - pav;

            for (int vtxId = 0; vtxId < 3; vtxId++) {
                double pv  = vtxPopAfs[vtxId];
                double qv  = 1 - pv;

                // Change values here for one allele
                double aaPev = log(pev);
                double bbPev = log(qev);

                double aaPfv = log(pfv);
                double bbPfv = log(qfv);

                double aaPav = log(pav);
                double bbPav = log(qav);
                
                double eGd = aaPev * pv + bbPev * qv;
                double fGd = aaPfv * pv + bbPfv * qv;
                double aGd = aaPav * pv + bbPav * qv;

                vtxExpGenoDists[vtxId][0][numSnps] = eGd;
                vtxExpGenoDists[vtxId][1][numSnps] = fGd;
                vtxExpGenoDists[vtxId][2][numSnps] = aGd;

                popExpPeSums[vtxId] += eGd;
                popExpPfSums[vtxId] += fGd;
                popExpPaSums[vtxId] += aGd;
            }

            delete[] refPopAfs;
            delete[] vtxPopAfs;

            numSnps++;
        }

        lineNo++;
    }
    fclose(ifp);

    //ASSERT(numSnps == numAncSnps, "numSnps = " << numAncSnps << ".\n");

    for (int vtxId = 0; vtxId < 3; vtxId++) {
        vtxPopExpGds[vtxId].e = -1 * popExpPeSums[vtxId]/numSnps;
        vtxPopExpGds[vtxId].f = -1 * popExpPfSums[vtxId]/numSnps;
        vtxPopExpGds[vtxId].a = -1 * popExpPaSums[vtxId]/numSnps;
    }

    if (print) Rprintf("Read %d ancestry SNPs from file %s\n", numSnps, ancSnpFile.c_str());
    
    return(0);
}

int AncestrySnps::parseSnpLine(char *snpLine, int nrefPop, int *pos, char *a1, char *a2, 
                               float *vtxPopAfs, float *refPopAfs) {

  int ret = 0;
  char *token;

  token = strtok(snpLine, "\t");
  if (token == NULL) return(1);
  *pos = stoi(string(token));
  token = strtok(NULL, "\t");
  if (token == NULL) return(2);
  *a1 = *token;
  token = strtok(NULL, "\t");
  if (token == NULL) return(3);
  *a2 = *token;
  for (int i=0; i<3; i++) {
    token = strtok(NULL, "\t");
    if (token == NULL) return(4+i);
    vtxPopAfs[i] = atof(token);
  }
  for (int i=0; i<nrefPop; i++) {
    token = strtok(NULL, "\t");
    if (token == NULL) return(4+i);
    refPopAfs[i] = atof(token);
  }

  return(ret);
}

int AncestrySnps::FindSnpIdGivenPos(int pos)
{
    int snpId       = -1;
    long int chrPos = (long) pos;
    int build       = 37;

    if (build == 37) {
        if (pos37ToAncSnpId.find(chrPos) != pos37ToAncSnpId.end()) {
            snpId = pos37ToAncSnpId[chrPos];
        }
    }
    

    return snpId;
}

AncestrySnp AncestrySnps::GetAncestrySnp(int snpId)
{
    return snps[snpId];
}

