#include "GrafPop.h"

extern "C" {

void C_main(char **cargs, int *iargs)
{
    int DEBUG = iargs[IARG_DEBUG];
    int print = iargs[IARG_PRINT];
    SampleGenoAncestry *smpGenoAnc = NULL;

    AncestrySnps *ancSnps = new AncestrySnps();
    if (DEBUG) Rprintf("Begin: ReadAncestrySnpsFromFile\n");
    ancSnps->ReadAncestrySnpsFromFile(cargs[CARG_ANCSNPFILE], print);
    if (DEBUG) Rprintf("End: ReadAncestrySnpsFromFile\n");

    //ancSnps->ShowAncestrySnps();

    int totAncSnps = ancSnps->GetNumAncestrySnps();
    int minAncSnps = iargs[IARG_MINANCSNPS];
    int fileType   = iargs[IARG_INFILE_TYPE];

    smpGenoAnc = new SampleGenoAncestry(ancSnps, minAncSnps);

    if (fileType == VCF_FILE) {
        VcfSampleAncestrySnpGeno *vcfGeno = new VcfSampleAncestrySnpGeno(cargs[CARG_INFILE], ancSnps);
        bool dataRead = vcfGeno->ReadDataFromFile(print);
        if (!dataRead) {
            error("Failed to read genotype data");
        }
        if (print) vcfGeno->ShowSummary();
        vcfGeno->RecodeSnpGenotypes();

        int numAncSnps = vcfGeno->vcfAncSnpIds.size();
        //int numVcfSmps = vcfGeno->GetNumSamples();

        if (smpGenoAnc->HasEnoughAncestrySnps(numAncSnps)) {
            smpGenoAnc->SetGenoSamples(vcfGeno->vcfSamples);
            smpGenoAnc->SetSnpGenoData(&vcfGeno->vcfAncSnpIds, &vcfGeno->vcfAncSnpCodedGenos);
        }
        else {
           error("Too few genotyped ancestry SNPs");
        }
    }
    else if (fileType == BED_FILE) {
        string bfile   = std::string(cargs[CARG_INFILE]); 
        string bedFile = bfile + ".bed";
        string bimFile = bfile + ".bim";
        string famFile = bfile + ".fam";

        FamFileSamples *famSmps = new FamFileSamples(famFile);
        if (print) famSmps->ShowSummary();

        smpGenoAnc->SetGenoSamples(famSmps->samples);
        //int numSmps = smpGenoAnc->GetNumSamples();

        BimFileAncestrySnps *bimSnps = new BimFileAncestrySnps(totAncSnps);
        bimSnps->ReadAncestrySnpsFromFile(bimFile, ancSnps, print);
        int numBimAncSnps = bimSnps->GetNumBimAncestrySnps();
        if (print) bimSnps->ShowSummary();

        if (smpGenoAnc->HasEnoughAncestrySnps(numBimAncSnps)) {
            BedFileSnpGeno *bedGenos = new BedFileSnpGeno(bedFile, ancSnps, bimSnps, famSmps);
            bool hasErr = bedGenos->ReadGenotypesFromBedFile(print);
            if (hasErr) error("ERROR reading genotype file");
            if (print) bedGenos->ShowSummary();

            smpGenoAnc->SetSnpGenoData(&bedGenos->ancSnpSnpIds, &bedGenos->ancSnpSmpGenos);
        }
        else {
            error("Too few genotyped ancestry SNPs");
        }
    }

    smpGenoAnc->SetAncestryPvalues(0, print);
    string outputFile = std::string(cargs[CARG_OUTFILE]);
    smpGenoAnc->SaveAncestryResults(outputFile, print);

    return ;
}

}
