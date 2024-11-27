#ifndef SAMPLE_GENO_ANCESTRY_H
#define SAMPLE_GENO_ANCESTRY_H

#include <fstream>
#include "Util.h"
#include "AncestrySnps.h"
#include "FamFileSamples.h"
#include "SampleGenoDist.h"

class GenoSample
{
public:
    string name;
    string father;
    string mother;
    int sex;
    int numRefPops;

    // Ancestry results calculated from genotypes
    int numAncSnps;
    bool ancIsSet;
    float gd1, gd2, gd3;
    float ePct, fPct, aPct;   // Ancestry (EUR, AFR, EAS) components of the sample
    fvec gdist;

public:
    GenoSample(string, int);
    void SetAncestryScores(int, float, float, float, dvec, float, float, float, bool);
};


class SampleGenoAncestry
{
public:
    int numRefPops;

private:
    int numSamples;
    int numAncSmps;

    int minAncSnps;
    int totAncSnps;
    int numAncSnps;
    int numThreads;                // Number of threads for parallel computing

    AncestrySnps *ancSnps;
    SampleGenoDist *vtxExpGd0;    // Genetic distances from 3 vertices to ref populations when all SNPs have genotypes

public:
    vector<GenoSample> samples;
    vector<int> *ancSnpIds;
    vector<char*> *ancSnpCodedGenos; // Use char, instead of int, to save space

    SampleGenoAncestry(AncestrySnps*, int, int=100);
    ~SampleGenoAncestry();

    void SetGenoSamples(const vector<string>&);
    void SetGenoSamples(const vector<FamSample>&);
    int SaveAncestryResults(string, int);
    void SetAncestryPvalues(int, int);
    void SetSnpGenoData(vector<int>*, vector<char*>*);
    void SetNumThreads(int);
    void InitPopPvalues();

    int GetNumSamples() { return numSamples; };
    int GetNumAncSamples() { return numAncSmps; };
    bool HasEnoughAncestrySnps(int numSnps) { return numSnps >= minAncSnps; }

    void ShowSummary();
};

#endif
