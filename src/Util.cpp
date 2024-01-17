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

