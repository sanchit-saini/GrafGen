#include "AncestrySnps.h"

AncestrySnp::AncestrySnp(int id, int ch, int g37, char a1, char a2, float* refPops, 
                        float *vtxPops, int nref)
{
    snpId      = id;
    chr        = ch;
    posG37     = g37;
    ref        = a1;
    alt        = a2;
    numRefPops = nref;

    refPopAfs = fvec(nref);
    for (int i = 0; i < numRefPops; i++) refPopAfs[i] = refPops[i];
    for (int i = 0; i < numVtxPops; i++) vtxPopAfs[i] = vtxPops[i];
}

AncestrySnps::AncestrySnps(int nref, int singleChr, char **map, int maplen, int nancsnps)
{
    // see order in ReadAncestrySnpsFromFile
    numRefPops = nref;
    onechr     = singleChr;
    snps       = {};
    chrmap     = map;  // passed in from R code
    chrmaplen  = maplen;
    
    vtxExpGenoDists = alloc_3d_double(numVtxPops, numVtxPops, nancsnps);
}

AncestrySnps::~AncestrySnps()
{
    snps.clear();
    pos37ToAncSnpId.clear();
    free_3d_double(vtxExpGenoDists, numVtxPops, numVtxPops);
}

int AncestrySnps::ReadAncestrySnpsFromFile(string ancSnpFile, int print)
{ 
    int nref = numRefPops;
    dvec popExpPfSums(nref);
    dvec popExpPaSums(nref);
    dvec popExpPeSums(nref);

    for (int popId = 0; popId < nref; popId++) {
        popExpPeSums[popId] = 0;
        popExpPfSums[popId] = 0;
        popExpPaSums[popId] = 0;
    }

    int lineLen = 5000;
    char snpLine[5000];

    FILE *ifp = fopen(ancSnpFile.c_str(), "r");
    if (!ifp) Rprintf("ERROR: Couldn't open genotype file\n");

    int lineNo = 0;
    bool fileIsValid = true;

    int numSnps = 0;
    int chr, pos;
    char a1, a2;

    while (fgets(snpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        if (lineNo == 0) {
            if (snpLine[0] != '.' || snpLine[1] != '.' || snpLine[2] != '.') {
                fileIsValid = false;
            }
            if (!fileIsValid) return(1);
        }
        else {
            float* refPopAfs = new float[nref];
            float* vtxPopAfs = new float[numVtxPops];
            
            int rc = parseSnpLine(snpLine, nref, &chr, &pos, &a1, &a2, 
                                  vtxPopAfs, refPopAfs);
            if (rc) {
              Rprintf("ERROR on line %d of ancestry SNP file", lineNo);
              continue;
            }
           
            AncestrySnp ancSnp(numSnps, chr, pos, a1, a2, refPopAfs, vtxPopAfs, nref);

            snps.push_back(ancSnp);

            if (onechr) {
              pos37ToAncSnpId[pos] = numSnps;
            } else {
              long int chrPos37 = (long)chr * 1000000000 + pos;
              pos37ToAncSnpId[chrPos37] = numSnps;
            }

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

int AncestrySnps::parseSnpLine(char *snpLine, int nrefPop, int *chr, int *pos, char *a1, char *a2, 
                               float *vtxPopAfs, float *refPopAfs) {

  int ret = 0;
  char *token;
  int ch;

  token = strtok(snpLine, "\t");
  if (token == NULL) return(-1);
  if (onechr) {
    //*chr = stoi(string(token));  // chr will be coded 0 from R code
    *chr  = 0;
  } else {
    ch = GetChrNumFromChrStr(token); 
    if (ch < 0) {
      Rprintf("ERROR with chr = %s\n", token);
      //error("1"); 
    }
    *chr = ch;
  }
  token = strtok(NULL, "\t");
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

int AncestrySnps::FindSnpIdGivenChrPos(int chr, int pos)
{
    int snpId = -1;

    long int chrPos = long(chr) * 1000000000 + pos;

    if (pos37ToAncSnpId.find(chrPos) != pos37ToAncSnpId.end()) {
      snpId = pos37ToAncSnpId[chrPos];
    }

    return snpId;
}


AncestrySnp AncestrySnps::GetAncestrySnp(int snpId)
{
    return snps[snpId];
}

int AncestrySnps::GetChrNumFromChrStr(char *chrstr)
{
    int ret = -1;
    for (int i=0; i<chrmaplen; i++) {
      if (strcmp(chrstr, chrmap[i]) == 0) {
        ret = i;
        break;
      }
    }
    if (ret < 0) {
      Rprintf("ERROR with chr = %s\n", chrstr);
    }
    return ret;
}

