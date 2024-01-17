#ifndef ANCESTRY_SNPS_H
#define ANCESTRY_SNPS_H

#include "Util.h"

class AncestrySnp
{
public:
    int snpId;
    //int rs;
    //int chr;
    int posG37;
    //int posG38;
    char ref;
    char alt;
    float vtxPopAfs[numVtxPops]; // E, F, A, or EUR, AFR, EAS
    float refPopAfs[numRefPops]; // EUR, AFA, EAS, AFR-near, EUR-near, AKL, ZAF

public:
    AncestrySnp(int, int, char, char, float*, float*);
};

class AncestrySnps
{
    map<long int, int> pos37ToAncSnpId;

public:
    AncestrySnps();
    ~AncestrySnps();
    vector<AncestrySnp> snps;
    // For each SNP, keeps the expected genetic distance from the 3 vertices to the 3 ref population
    double vtxExpGenoDists[numVtxPops][numVtxPops][numAncSnps];
    // Vertex genetic distances summed up using all ancestry SNPs
    GenoDist vtxPopExpGds[numVtxPops];

    //string refPopNames[numRefPops];

    int ReadAncestrySnpsFromFile(string, int);
    int FindSnpIdGivenPos(int);

    int FindSnpIdGivenRs(int);
    int FindSnpIdGivenChrPos(int, int, int);
    AncestrySnp GetAncestrySnp(int);
    void SetVertexExpecteGeneticDists();
    int GetNumAncestrySnps() { return snps.size(); };
    void ShowAncestrySnps();
    int parseSnpLine(char *, int , int *, char *, char *, float *, float *);
};

#endif
