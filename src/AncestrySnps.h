#ifndef ANCESTRY_SNPS_H
#define ANCESTRY_SNPS_H

#include "Util.h"

class AncestrySnp
{
public:
    int snpId;
    //int rs;
    int chr;
    int posG37;
    //int posG38;
    int numRefPops;
    char ref;
    char alt;
    float vtxPopAfs[numVtxPops]; // E, F, A, or EUR, AFR, EAS
    fvec refPopAfs; // EUR, AFA, EAS, AFR-near, EUR-near, AKL, ZAF

public:
    AncestrySnp(int, int, int, char, char, float*, float*, int);
};

class AncestrySnps
{
public:
    int numRefPops;
    int onechr;
    char **chrmap;
    int chrmaplen;

    map<long int, int> pos37ToAncSnpId;

public:
    AncestrySnps(int, int, char **, int, int);
    ~AncestrySnps();
    vector<AncestrySnp> snps;
    // For each SNP, keeps the expected genetic distance from the 3 vertices to the 3 ref population
    //double vtxExpGenoDists[numVtxPops][numVtxPops][numAncSnps];
    double ***vtxExpGenoDists; // numAncSnps is no longer fixed

    // Vertex genetic distances summed up using all ancestry SNPs
    GenoDist vtxPopExpGds[numVtxPops];

    int ReadAncestrySnpsFromFile(string, int);
    int FindSnpIdGivenPos(int);

    int FindSnpIdGivenRs(int);
    int FindSnpIdGivenChrPos(int, int);
    AncestrySnp GetAncestrySnp(int);
    void SetVertexExpecteGeneticDists();
    int GetNumAncestrySnps() { return snps.size(); };
    void ShowAncestrySnps();
    int parseSnpLine(char *, int , int *, int *, char *, char *, float *, float *);
    int GetChrNumFromChrStr(char *);
};

#endif
