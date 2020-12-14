#include "KmersWrapper.h"

std::map<std::string, int> KmersWrapper::get_hamming_distances(std::string kmer) {
    return KmersWrapper::hamming_distance_matrix[kmer];
}

int KmersWrapper::get_hamming_distance(std::string kmer1,std::string kmer2){
    return KmersWrapper::hamming_distance_matrix[kmer1][kmer2];
}


KmersWrapper::KmersWrapper(std::map<std::string, std::map<std::string, int>> hd, Json::Value oc, Json::Value ec, std::unordered_set<std::string> v, std::unordered_set<std::string> o, std::string itype) {
    hamming_distance_matrix = hd;
    observedCounts = oc;
    expectedCounts = ec;
    V = v;
    O = o;
    itersetType = itype;
    iterset = getIterset();
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