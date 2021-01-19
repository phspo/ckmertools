#include "KmersWrapper.h"
#include<cmath>

std::map<std::string, int> KmersWrapper::get_hamming_distances(std::string kmer) {
    return KmersWrapper::hamming_distance_matrix[kmer];
}

int KmersWrapper::get_hamming_distance(std::string kmer1,std::string kmer2){
    return KmersWrapper::hamming_distance_matrix[kmer1][kmer2];
}


KmersWrapper::KmersWrapper(std::map<std::string, std::map<std::string, int>> hd, Json::Value oc, Json::Value ec, std::unordered_set<std::string> v, std::unordered_set<std::string> o, std::string itype, float kmerError) {
    hamming_distance_matrix = hd;
    observedCounts = oc;
    expectedCounts = ec;
    V = v;
    O = o;
    itersetType = itype;
    iterset = getIterset();
    int max_hd = 5;
    float *a = (float*) malloc(sizeof(float) * max_hd);
    pre_compute_hd_probabilities(kmerError, max_hd);
};

enum ItersetOptions {
    InvalidType,
    OType,
    VType,
    OuVType,
    OnVType
};

ItersetOptions resolveItersetOption(std::string &input) {
    if( input == "O" ) return OType;
    if( input == "V" ) return VType;
    if( input == "OuV" ) return OuVType;
    if( input == "OnV" ) return OnVType;
    return InvalidType;
}

void KmersWrapper::pre_compute_hd_probabilities(float kmerError, int max_hd) {
    int length_kmer = (*begin(O)).length();
    a[0] = 1;
    for (int i = 1; i < max_hd; i++)
    {
        //TODO: replace with mutation rate
        // i is hamming distance
        a[i] = std::pow(kmerError,i)*std::pow((1 - kmerError),(length_kmer - i));
    }
}

float KmersWrapper::get_computed_probability(int i) {
    return a[i];
}

std::unordered_set<std::string> KmersWrapper::getIterset() {
    std::set<std::string> ordO(O.begin(), O.end());
    std::set<std::string> ordV(V.begin(), V.end());

    switch( resolveItersetOption(itersetType) )
    {
        case OType: {
            return O;
        }
        case VType: {
            return V;
        }
        case OuVType: {
            // O Union V
            std::unordered_set <std::string> OuV;
            std::set_union(ordO.begin(), ordO.end(), ordV.begin(), ordV.end(), std::inserter(OuV, OuV.begin()));
            return OuV;
        }
        case OnVType: {
            // O Intersect V
            std::unordered_set <std::string> OnV;
            std::set_intersection(ordO.begin(), ordO.end(), ordV.begin(), ordV.end(), std::inserter(OnV, OnV.begin()));
            return OnV;
        }
        default: {
            return O;
        }
    }
} 