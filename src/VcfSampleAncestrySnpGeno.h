#ifndef VCF_SAMPLE_ANCESTRY_SNP_GENO_H
#define VCF_SAMPLE_ANCESTRY_SNP_GENO_H

#include <zlib.h>
#include <errno.h>
#include "Util.h"
#include "AncestrySnps.h"

//#define BUFFERLEN 0x0010
#define BUFFERLEN 1234
#define WORDLEN 10000

class VcfSampleAncestrySnpGeno
{
public:
    int onechr;
    int prophage;

private:
    string vcfFile;
    AncestrySnps *ancSnps;

    // For all the following arrays:
    // One record for each putative ancestry SNP (checked using pos)

    // Saves the one allele id number for each sample (0=ref, 1=first alternate allele, etc) 
    //vector<vector<char>> vcfAncSnpGtRefs;
    vector<vector<char>> vcfAncSnpGtAlts;

    vector<int> vcfAncSnpChrs;      // chr value from The CHROM string
    vector<int> vcfAncSnpPoss;      // pos value from  POS string
    //vector<string> vcfAncSnpSnps;   // The ID string
    vector<string> vcfAncSnpRefs;   // The REF string
    vector<string> vcfAncSnpAlts;   // The ALT string
    //vector<int> vcfRsIdAncSnpIds;   // Ancestry SNP ID derived using RS ID
    vector<int> vcfGb37AncSnpIds;   // Ancestry SNP ID derived using Build 37 chr + pos
    //vector<int> vcfGb38AncSnpIds;   // Ancestry SNP ID derived using Build 38 chr + pos

    int totAncSnps;
    int numSamples;
    int totVcfSnps;
    int putativeAncSnps;
    int numRsIdAncSnps;
    int numGb37AncSnps;
    int numGb38AncSnps;
    int numVcfAncSnps;
    AncestrySnpType ancSnpType;

    void CompareAncestrySnpAlleles(const string, const string, const char, const char, int*, int*);
    void CompareAncestrySnpAlleles_prophage(const string, const string, const char, const char, int*, int*);
    int RecodeGenotypeGivenString(const int, const int, const string);
    int RecodeGenotypeGivenIntegers(const int, const int, const int);

public:
    // Each genotype (per SNP and sample) is coded with number of alts, i.e., 0 = ref, 1 = alternate allele, 3 = unknown
    // Each row in the array is the genotypes of one SNP, coded as  an array of integers for all samples
    // Ancestry SNP ID of each row is saved in the array of vcfAncSnpIds
    vector<string> vcfSamples;
    vector<int> vcfAncSnpIds;
    vector<char*> vcfAncSnpCodedGenos; // Use char, instead of int, to save space

    VcfSampleAncestrySnpGeno(string, AncestrySnps*, int, int);
    ~VcfSampleAncestrySnpGeno();

    int GetNumVcfSnps() { return totVcfSnps; };
    int GetNumSamples() { return numSamples; };
    int GetNumVcfAncestrySnps() { return numVcfAncSnps; };
    bool ReadDataFromFile(int);
    void RecodeSnpGenotypes();

    void ShowSummary();
    void DeleteAncSnpGtValues();
    void DeleteAncSnpCodedGenos();
};


#endif
