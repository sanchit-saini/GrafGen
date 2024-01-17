#include "SampleGenoAncestry.h"

GenoSample::GenoSample(string smp)
{
    name = smp;
    father = "";
    mother = "";
    sex = 0;

    numAncSnps = 0;
    ancIsSet = false;
}

void GenoSample::SetAncestryScores(int numSnps, float d1, float d2, float d3, double *dist, float e, float f, float a, bool isAnc)
{
    numAncSnps = numSnps;
    ancIsSet = isAnc;
    gd1 = d1;
    gd2 = d2;
    gd3 = d3;
    for (int i=0; i < numRefPops; i++) gdist[i] = (float) dist[i];
    ePct = e;
    fPct = f;
    aPct = a;
}

SampleGenoAncestry::SampleGenoAncestry(AncestrySnps *aSnps, int minSnps)
{
    ancSnps = aSnps;
    if (minSnps) minAncSnps = minSnps;
    else         minAncSnps = 100;
    numSamples = 0;
    numAncSnps = 0;
    totAncSnps = ancSnps->GetNumAncestrySnps();

    ancSnpIds = NULL;
    ancSnpCodedGenos = NULL;

    samples = {};

    numThreads = 1;
    vtxExpGd0 = new SampleGenoDist(&aSnps->vtxPopExpGds[0], &aSnps->vtxPopExpGds[1],
    &aSnps->vtxPopExpGds[2], &aSnps->vtxPopExpGds[0]);
    vtxExpGd0->TransformAllDists();
    vtxExpGd0->CalculateBaryCenters();
}

SampleGenoAncestry::~SampleGenoAncestry()
{
    delete vtxExpGd0;
    samples.clear();
}

void SampleGenoAncestry::SetNumThreads(int threads)
{
    numThreads = threads;
}

void SampleGenoAncestry::SetGenoSamples(const vector<string> &smps)
{
    if (!smps.empty()) {
        numSamples = smps.size();

        for (int i = 0; i < numSamples; i++) {
            GenoSample genoSmp = GenoSample(smps[i]);
            samples.push_back(genoSmp);
        }
    }
    numAncSmps = 0;
}

void SampleGenoAncestry::SetGenoSamples(const vector<FamSample> &smps)
{
    if (!smps.empty()) {
        numSamples = smps.size();

        for (int i = 0; i < numSamples; i++) {
            GenoSample genoSmp = GenoSample(smps[i].name);
            samples.push_back(genoSmp);
        }
    }
    numAncSmps = 0;
}

void SampleGenoAncestry::SetSnpGenoData(vector<int> *snpIds, vector<char*> *snpCodedGenos)
{
    ancSnpIds = snpIds;
    ancSnpCodedGenos = snpCodedGenos;
    numAncSnps = ancSnpIds->size();
}

int SampleGenoAncestry::SaveAncestryResults(string outFile, int print)
{
    int numSaveSmps = 0;
    for (int i = 0; i < numSamples; i++) {
        if (samples[i].ancIsSet) numSaveSmps++;
    }

    if (numSaveSmps < 1) {
      if (print) {
        Rprintf("NOTE: None of the %d samples have enough genotypes for ancestry inference.",
                numSamples);
        Rprintf("No ancestry results were generated.\n");
      }
      return numSaveSmps;
    }

    //string vtxTitle = "Positions (x, y, z coordinates) of the three vertices";
    //vtxExpGd0->ShowPositions(vtxTitle);

    FILE *ifp = fopen(outFile.c_str(), "w");
    if(ifp) {
       
        fprintf(ifp, "# Positions of the three vertices\n");
        fprintf(ifp, "#\n");

        fprintf(ifp, "#          x       y      z\n");

        fprintf(ifp, "# F: \t%5.4f  %5.4f %5.4f\n", 
                vtxExpGd0->fPt.x, vtxExpGd0->fPt.y, vtxExpGd0->fPt.z);

        fprintf(ifp, "# A: \t%5.4f  %5.4f %5.4f\n", 
                vtxExpGd0->aPt.x, vtxExpGd0->aPt.y, vtxExpGd0->aPt.z);

        fprintf(ifp, "# E: \t%5.4f  %5.4f %5.4f\n", 
                 vtxExpGd0->ePt.x, vtxExpGd0->ePt.y, vtxExpGd0->ePt.z);
        fprintf(ifp, "#\n");

        fprintf(ifp, "%s\t%s\tGD1_x\tGD2_y\tGD3_z\tE_percent\tF_percent\tA_percent\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
                "Sample", "N_SNPs", REFPOP0, REFPOP1, REFPOP2, REFPOP3, REFPOP4, REFPOP5, REFPOP6, REFPOP7, REFPOP8);

        for (int i = 0; i < numSamples; i++) {
            GenoSample smp = samples[i];
            if (!smp.ancIsSet) continue;

            fprintf(ifp, "%s\t%d\t%7.6f\t%7.6f\t%7.6f\t%6.2f\t%6.2f\t%6.2f\t%7.6f\t%7.6f\t%7.6f\t%7.6f\t%7.6f\t%7.6f\t%7.6f\t%7.6f\t%7.6f\n",
                    smp.name.c_str(), smp.numAncSnps,
                    smp.gd1, smp.gd2, smp.gd3, smp.ePct, smp.fPct, smp.aPct,
                    smp.gdist[0], smp.gdist[1], smp.gdist[2], smp.gdist[3],
                    smp.gdist[4], smp.gdist[5], smp.gdist[6], smp.gdist[7], smp.gdist[8]);
        }
    }
    else {
      if (print) Rprintf("ERROR: Can't open %s for writing!\n", outFile.c_str());
      error("ERROR opening output file");
      return 0;
    }

    fclose(ifp);
    //Rprintf("Saved population results of %d samples to %s.\n", numSaveSmps, outFile);

    return numSaveSmps;
}

void SampleGenoAncestry::SetAncestryPvalues(int thNo, int print)
{
    //bool debug = numAncSnps < 0 && numSamples < 0 ? true : false;

    int meanThSmps = int(numSamples / numThreads);
    int rmSmps = numSamples % numThreads;
    int chkThSmps = meanThSmps;
    if (thNo+1 <= rmSmps) chkThSmps = meanThSmps + 1;

    int stSmp = thNo * chkThSmps;
    if (thNo >= rmSmps) stSmp = thNo * chkThSmps + rmSmps;
    int edSmp = stSmp + chkThSmps - 1;

    int smpCnt = 0;
    for (int smpNo = stSmp; smpNo <= edSmp; smpNo++) {
        // Calculate 16 scores for each sample, based on the SNPs with genotypes, i.e.,
        // 9 expected genetics distances from the 3 vertices to the first 3 reference populations, and
        // 7 genetic distances from the sample to the 7 referene populations

        //if (debug && thNo == 0) cout << "Smp No. " << smpNo << ". smp " << samples[smpNo].name << " ";

        double popPvalues[numRefPops];
        double popMeanPvals[numRefPops];

        for (int popId = 0; popId < numRefPops; popId++) {
            popPvalues[popId] = 0;   // The raw log p-value
            popMeanPvals[popId] = 0; // Genetic distancs from the sample to each ref population
        }

        // Assuming the three vertex populations have allele freqs for all Ancestry SNPs
        int numGenoSnps = 0;
        double vtxExpPeSums[3];
        double vtxExpPfSums[3];
        double vtxExpPaSums[3];

        // Initialize expected values for each vertex population (E, F, A)
        for (int vtxId = 0; vtxId < 3; vtxId++) {
            vtxExpPfSums[vtxId] = 0;
            vtxExpPeSums[vtxId] = 0;
            vtxExpPaSums[vtxId] = 0;
        }

        // Assuming some reference populations might not have allele freqs for all Ancestry SNPs,
        // which is not ture for the current 5 reference populations being used
        int refPopSnps[numRefPops]; // Counts of SNPs with freqs for each ref population
        for (int popId = 0; popId < numRefPops; popId++) {
            refPopSnps[popId] = 0;
        }

        //bool hasRefPv = false;
        int snpNo = 0;
        for (snpNo = 0; snpNo < numAncSnps; snpNo++) {
            int geno = (*ancSnpCodedGenos)[snpNo][smpNo];
            int ancSnpId = (*ancSnpIds)[snpNo];

            // Alt allele freq p of the 3 vertices
            //double v0p = ancSnps->snps[ancSnpId].vtxPopAfs[0];
            //double v1p = ancSnps->snps[ancSnpId].vtxPopAfs[1];
            //double v2p = ancSnps->snps[ancSnpId].vtxPopAfs[2];

            if (geno > -1 && geno < 3) {
                /*
                if (debug && 0) {
                    cout << "    Snp No. " << snpNo << " ID " << ancSnpId
                    << " pos " << ancSnps->snps[ancSnpId].posG37
                    << " ref " << ancSnps->snps[ancSnpId].ref
                    << " alt " << ancSnps->snps[ancSnpId].alt
                    << " v0 " << ancSnps->snps[ancSnpId].vtxPopAfs[0]
                    << " p0 " << ancSnps->snps[ancSnpId].refPopAfs[0]
                    << " geno " << geno << "\n";
                }
                */
                for (int popId = 0; popId < numRefPops; popId++) {
                    double pv = ancSnps->snps[ancSnpId].refPopAfs[popId]; // reference population allele freq p value

                    if (pv > 0 && pv < 1) {
                        double qv = 1 - pv;

                        // Change values for one allele
                        double bbPv = log(qv);
                        double aaPv = log(pv);
                        if      (geno == 1) {
                            popPvalues[popId] += bbPv;
                        }
                        else if (geno == 0) {
                            popPvalues[popId] += aaPv;
                        }

                        /* Original code 
                        double bbPv = log(qv) * 2;
                        double abPv = log(pv * qv * 2);
                        double aaPv = log(pv) * 2;
                        if      (geno == 2) {
                            popPvalues[popId] += bbPv;
                        }
                        else if (geno == 1) {
                            popPvalues[popId] += abPv;
                        }
                        else if (geno == 0) {
                            popPvalues[popId] += aaPv;
                        }
                        */

                        //hasRefPv = true;
                        refPopSnps[popId]++;
                    }
                }

                // Calculate expected P values for each of the vertices E, F, A
                for (int vtxId = 0; vtxId < 3; vtxId++) {
                    vtxExpPeSums[vtxId] += ancSnps->vtxExpGenoDists[vtxId][0][ancSnpId];
                    vtxExpPfSums[vtxId] += ancSnps->vtxExpGenoDists[vtxId][1][ancSnpId];
                    vtxExpPaSums[vtxId] += ancSnps->vtxExpGenoDists[vtxId][2][ancSnpId];
                }

                numGenoSnps++;
            }
        }

        for (int popId = 0; popId < numRefPops; popId++) {
            if (refPopSnps[popId] > 0) {
                popMeanPvals[popId] = -1 * popPvalues[popId]/refPopSnps[popId];
            }
        }

        float gd1 = 0, gd2 = 0, gd3 = 0;
        float ePct = 0, fPct = 0, aPct = 0;
        bool hasAncGeno = false;

        if (numGenoSnps >= minAncSnps) {
            GenoDist smpDist;
            GenoDist vtxExpDists[numVtxPops];

            smpDist.e = popMeanPvals[0];
            smpDist.f = popMeanPvals[1];
            smpDist.a = popMeanPvals[2];

            for (int vtxId = 0; vtxId < numVtxPops; vtxId++) {
                vtxExpDists[vtxId].e  = -1 * vtxExpPeSums[vtxId]/numGenoSnps;
                vtxExpDists[vtxId].f  = -1 * vtxExpPfSums[vtxId]/numGenoSnps;
                vtxExpDists[vtxId].a  = -1 * vtxExpPaSums[vtxId]/numGenoSnps;
            }

            // Calculate GD and ancestry components using the raw scores
            SampleGenoDist *smpGd = new SampleGenoDist(&vtxExpDists[0], &vtxExpDists[1], &vtxExpDists[2], &smpDist);
            smpGd->TransformAllDists();
            smpGd->CalculateBaryCenters();

            // Show rotated x, y, z values as GD1, GD2, GD3
            gd1 = smpGd->eWt * vtxExpGd0->ePt.x + smpGd->fWt * vtxExpGd0->fPt.x + smpGd->aWt * vtxExpGd0->aPt.x;
            gd2 = smpGd->eWt * vtxExpGd0->ePt.y + smpGd->fWt * vtxExpGd0->fPt.y + smpGd->aWt * vtxExpGd0->aPt.y;
            gd3 = smpGd->sPt.z;

            double ejWt = smpGd->eWt > 0 ? smpGd->eWt : 0;
            double fjWt = smpGd->fWt > 0 ? smpGd->fWt : 0;
            double ajWt = smpGd->aWt > 0 ? smpGd->aWt : 0;
            double totWt = fjWt + ejWt + ajWt;
            ePct = ejWt * 100 / totWt;
            fPct = fjWt * 100 / totWt;
            aPct = ajWt * 100 / totWt;

            hasAncGeno = true;
            numAncSmps++;
            delete smpGd;
        }

        samples[smpNo].SetAncestryScores(numGenoSnps, gd1, gd2, gd3, popMeanPvals, ePct, fPct, aPct, hasAncGeno);

        smpCnt++;
        if (print && (thNo == 0) && (smpCnt % 100 == 0))
          Rprintf("Calculated scores for %d of %d samples\n", 
                  smpCnt, chkThSmps);
    }
}

/*
void SampleGenoAncestry::ShowSummary()
{
    cout << "Num samples: " << numSamples << "\n";

    cout << "Tot ancestry SNPs: " << totAncSnps << "\n";
    cout << "Num ancestry SNPs in dataset: " << numAncSnps << "\n";

    for (int i = 0; i < numAncSnps; i++) {
        int ancSnpId = (*ancSnpIds)[i];
        cout << "No. " << i << ": " << ancSnpId << ": ";

        for (int j = 0; j < numSamples; j++) {
            if  (j < 20) cout << int((*ancSnpCodedGenos)[i][j]) << " ";
        }
        cout << "\n";

        if (i > 20) break;
    }
}
*/
