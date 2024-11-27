#include "FamFileSamples.h"

using namespace std;


FamSample::FamSample(string smp, string dad, string mom, int gender)
{
    name = smp;
    father = dad;
    mother = mom;
    sex = gender;
}


FamFileSamples::FamFileSamples(string file)
{
    filename = file;
    numFamSmps = 0;
    numMales = 0;
    numFemales = 0;

    ReadSamplesFromFile();
}

bool FamFileSamples::Summarize()
{
    bool success = false;
    int arraySize = samples.size();

    if (arraySize == numFamSmps) {
        for (int i = 0; i < numFamSmps; i++) {
            FamSample smp = samples[i];
            if (smp.sex == 1) {
                numMales++;
            }
            else if (smp.sex == 2) {
                numFemales++;
            }
        }

        success = true;
    }

    return success;
}

void FamFileSamples::ShowSummary()
{
    //bool success = false;
    Rprintf("Total %d samples in fam file %s .\n",
            numFamSmps, filename.c_str());

}

int FamFileSamples::ReadSamplesFromFile()
{
    int numFileSmps = 0;

    int lineLen = 300;
    char fpLine[300];

    FILE *ifp = fopen(filename.c_str(), "r");
    if (!ifp) Rprintf("ERROR: Couldn't open fam file");

    int lineNo = 0;
    bool fileIsValid = true;

    int smpSex, pheno;
    char famId[80], smpId[80], dadId[80], momId[80];

    while (fgets(fpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        sscanf(fpLine, "%s %s %s %s %d %d", famId, smpId, dadId, momId, &smpSex, &pheno);

        if (!smpSex) smpSex = 0;

        FamSample sample(smpId, dadId, momId, smpSex);
        samples.push_back(sample);
        numFileSmps++;
        
        lineNo++;
    }
    fclose(ifp);

    numFamSmps = numFileSmps;

    Summarize();

    return numFileSmps;
}

