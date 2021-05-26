#include "KmersWrapper.h"
#include<cmath>
#include "parsing.h"
#include <iostream>


std::map<std::string, int> KmersWrapper::get_hamming_distances(std::string kmer) {
    //if (hamming_distance_matrix.find(kmer) != hamming_distance_matrix.end()) {
    return KmersWrapper::hamming_distance_matrix[kmer];
    //}
    //return KmersWrapper::ohamming_distance_matrix[kmer]; 
}

int KmersWrapper::get_hamming_distance(std::string kmer1,std::string kmer2){
    return KmersWrapper::hamming_distance_matrix[kmer1][kmer2];
}


KmersWrapper::KmersWrapper(std::string hammingdist, std::string kmersindex, std::string ukmersindex, std::string observed, std::string expected, std::string itype, float kmerError, std::string ohammingdist, std::string okmersindex) {
    hamming_distance_matrix = parsing::get_hammingdistances(hammingdist, kmersindex);
    std::cout <<"Init ohamming_distance_matrix \n";
    ohamming_distance_matrix = parsing::get_hammingdistancesO(ohammingdist, okmersindex, kmersindex);
    observedCounts = parsing::readDictionary(observed);
    expectedCounts = parsing::readDictionary(expected);
    V = parsing::get_V(expectedCounts);
    U = parsing::get_U(ukmersindex);
    O = parsing::get_O(observedCounts);
    itersetType = itype;
    iterset = getIterset();
    a = (float *) malloc(5*sizeof(float));
    pre_compute_hd_probabilities(kmerError, 5);
};

KmersWrapper::~KmersWrapper() {
    delete a;
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
    a[0] = (float) 1;
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