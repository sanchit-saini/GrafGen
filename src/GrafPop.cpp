#include "GrafPop.h"

extern "C" {

void C_main(char **cargs, int *iargs, char **chrmap)
{
    int DEBUG     = iargs[IARG_DEBUG];
    int print     = iargs[IARG_PRINT];
    int nref      = iargs[IARG_NREFPOP];
    int onechr    = iargs[IARG_ONECHR];
    int chrmaplen = iargs[IARG_CHRMAPLEN];
    int nancsnps  = iargs[IARG_NANCSNPS];
    int prophage  = iargs[IARG_PROPHAGE];
    if (DEBUG) Rprintf("nref=%d, onechr=%d, chrmaplen=%d, nancsnps=%d, prophage=%d\n",
                       nref, onechr, chrmaplen, nancsnps, prophage);
    if (onechr && chrmaplen) Rprintf("onechr and chrmaplen");
    SampleGenoAncestry *smpGenoAnc = NULL;
    AncestrySnps *ancSnps = new AncestrySnps(nref, onechr, chrmap, chrmaplen, nancsnps);
    if (DEBUG) Rprintf("Begin: ReadAncestrySnpsFromFile\n");
    ancSnps->ReadAncestrySnpsFromFile(cargs[CARG_ANCSNPFILE], print);
    if (DEBUG) Rprintf("End: ReadAncestrySnpsFromFile\n");
    //ancSnps->ShowAncestrySnps();

    int totAncSnps = ancSnps->GetNumAncestrySnps();
    int minAncSnps = iargs[IARG_MINANCSNPS];
    int fileType   = iargs[IARG_INFILE_TYPE];

    smpGenoAnc = new SampleGenoAncestry(ancSnps, nref, minAncSnps);

    if (fileType == VCF_FILE) {
        VcfSampleAncestrySnpGeno *vcfGeno = new VcfSampleAncestrySnpGeno(cargs[CARG_INFILE], ancSnps, onechr, prophage);
        if (DEBUG) Rprintf("Begin: ReadDataFromFile\n");
        bool dataRead = vcfGeno->ReadDataFromFile(print);
        if (DEBUG) Rprintf("End: ReadDataFromFile\n");
        if (!dataRead) {
            Rprintf("Failed to read genotype data");
        }
        if (print) vcfGeno->ShowSummary();
        vcfGeno->RecodeSnpGenotypes();

        int numAncSnps = vcfGeno->vcfAncSnpIds.size();
        if (DEBUG) Rprintf("numAncSnps = %d\n", numAncSnps);
        //int numVcfSmps = vcfGeno->GetNumSamples();
        if (smpGenoAnc->HasEnoughAncestrySnps(numAncSnps)) {
            smpGenoAnc->SetGenoSamples(vcfGeno->vcfSamples);
            smpGenoAnc->SetSnpGenoData(&vcfGeno->vcfAncSnpIds, &vcfGeno->vcfAncSnpCodedGenos);
        }
        else {
           Rprintf("ERROR: Too few genotyped ancestry SNPs\n");
           Rprintf("Too few genotyped ancestry SNPs");
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

        BimFileAncestrySnps *bimSnps = new BimFileAncestrySnps(totAncSnps, onechr);
        bimSnps->ReadAncestrySnpsFromFile(bimFile, ancSnps, print);
        int numBimAncSnps = bimSnps->GetNumBimAncestrySnps();
        if (print) bimSnps->ShowSummary();

        if (smpGenoAnc->HasEnoughAncestrySnps(numBimAncSnps)) {
            BedFileSnpGeno *bedGenos = new BedFileSnpGeno(bedFile, ancSnps, bimSnps, famSmps);
            bool hasErr = bedGenos->ReadGenotypesFromBedFile(print);
            if (hasErr) Rprintf("ERROR reading genotype file");
            if (print) bedGenos->ShowSummary();

            smpGenoAnc->SetSnpGenoData(&bedGenos->ancSnpSnpIds, &bedGenos->ancSnpSmpGenos);
        }
        else {
            Rprintf("Too few genotyped ancestry SNPs");
        }
    }

    smpGenoAnc->SetAncestryPvalues(0, print);
    string outputFile = std::string(cargs[CARG_OUTFILE]);
    smpGenoAnc->SaveAncestryResults(outputFile, print);

    return ;
}

}
