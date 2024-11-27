#include "Util.h"

char FlipAllele(char allele)
{
    char flipAllele = '0';

    switch(allele) {
        case 'A': flipAllele = 'T'; break;
        case 'T': flipAllele = 'A'; break;
        case 'G': flipAllele = 'C'; break;
        case 'C': flipAllele = 'G'; break;
    }

    return flipAllele;
}

int AmbigAlleles(char a1, char a2)
{
  int ret = 0;
  
  if (((a1 == 'A') && (a2 == 'T')) || ((a1 == 'T') && (a2 == 'A')) ||
     ((a1 == 'C') && (a2 == 'G')) || ((a1 == 'G') && (a2 == 'C'))) ret = 1;

  return(ret);
}

vector<string> SplitString(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;

    while (pos < str.size() && prev < str.size()) {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.size();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.size();
    }

    return tokens;
}

int GetChromosomeFromString(const char* chrStr)
{
    int chrNum = 0;
    int i = 0;

    if (strlen(chrStr) > 3 &&
        (chrStr[0] == 'c' || chrStr[0] == 'C') &&
        (chrStr[1] == 'h' || chrStr[1] == 'H') &&
        (chrStr[2] == 'r' || chrStr[2] == 'R') ) {
        i = 3;
    }

    while(chrStr[i] != 0) {
        int num = chrStr[i] - '0';

        if (num >= 0 && num < 10) {
            chrNum = chrNum * 10 + num;
        }
        else {
            chrNum = 0;
            break;
        }

        i++;
    }

    if (chrNum < 1 || chrNum > 22) chrNum = 0;

    return chrNum;
}

double** alloc_2d_double(int n1, int n2)
{
  double **ret;
  int i;
  ret = new double*[n1];
  if (ret == NULL) Rprintf("Out of memory in alloc_2d_double");
  for(int i = 0; i < n1; i++) {
    ret[i] = new double[n2];
    if (ret[i] == NULL) Rprintf("Out of memory in alloc_2d_double");
  }
  return(ret);
}

double*** alloc_3d_double(int n1, int n2, int n3)
{
  double ***ret;
  int i;
  ret = new double**[n1];
  if (ret == NULL) Rprintf("Out of memory in alloc_3d_double");
  for(int i = 0; i < n1; i++) {
    ret[i] = alloc_2d_double(n2, n3);
  }
  return(ret);
}

void free_2d_double(double **mat, int nr) {

  for(int i = 0; i < nr; i++) delete[] mat[i];
  delete[] mat;

}

void free_3d_double(double ***mat, int n1, int n2) {

  for(int i = 0; i < n1; i++) free_2d_double(mat[i], n2);
  delete[] mat;

}

