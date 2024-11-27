
#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <map>
#include <vector>
#include <unistd.h>
#include <R.h>
#include <R_ext/Print.h>
#include "zlib.h"

//#include <Rmath.h>
//#include <R_ext/Memory.h>
//#include <Rinternals.h>
//#include <R_ext/Rdynload.h>


const double pi = 3.1415926;

// bi-allele, MAF=0.01, R2=0 
//static const int numAncSnps = 143733;

static const int numVtxPops = 3;

#define REFPOP0 "hpgpEurope"
#define REFPOP1 "hpgpAfrica"
#define REFPOP2 "hpgpAsia"
#define REFPOP3 "hpgpAfroamerica"
#define REFPOP4 "hpgpEuroamerica"
#define REFPOP5 "hpgpMediterranea"
#define REFPOP6 "hpgpEurasia"
#define REFPOP7 "hpgpAklavik86-like"
#define REFPOP8 "hpgpAfrica-distant"

#define CARG_INFILE 0
#define CARG_OUTFILE 1
#define CARG_ANCSNPFILE 2
#define IARG_INFILE_TYPE 0
#define IARG_MINANCSNPS 1
#define IARG_DEBUG 2
#define IARG_PRINT 3
#define IARG_NANCSNPS 4
#define IARG_NREFPOP 5
#define IARG_ONECHR 6
#define IARG_CHRMAPLEN 7
#define IARG_PROPHAGE 8

#define BED_FILE 0
#define VCF_FILE 1

using namespace std;

typedef vector<float> fvec;
typedef vector<double> dvec;
typedef vector<int> ivec;

enum class AncestrySnpType
{
    RSID = 0,
    GB37 = 1,
};

enum class GenoDatasetType
{
    NOT_EXISTS  = 0,
    IS_PLINK    = 1,
    IS_PLINK_GZ = 2,
    IS_VCF      = 3,
    IS_VCF_GZ   = 4,
    IS_OTHER    = 5
};

// Define Genetic Distances to the three reference populations
struct GenoDist
{
    double e; // To European
    double f; // To African
    double a; // To East Asian
};

// A 3-D point in space
struct Point
{
    double x;
    double y;
    double z;
};

char FlipAllele(char);
int AmbigAlleles(char, char);
vector<string> SplitString(const string&, const string&);
int GetChromosomeFromString(const char*);

double** alloc_2d_double(int, int);
double*** alloc_3d_double(int, int, int);
void free_2d_double(double**, int);
void free_3d_double(double***, int, int);

#endif
