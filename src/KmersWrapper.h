// KmersWrapper.h

#include <iostream>
#include <map>
#include <string>
#include <iterator>
#include "jsoncpp/json/json.h"
#include <unordered_set>
#include <set>
#include <algorithm>

class KmersWrapper {
    private:
    std::unordered_set<std::string> getIterset();
    public:
    std::map<std::string, std::map<std::string, int>> hamming_distance_matrix;
    Json::Value observedCounts;
    Json::Value expectedCounts;
    std::unordered_set<std::string> V;
    std::unordered_set<std::string> O;
    std::unordered_set<std::string> iterset;
    std::string itersetType;

    KmersWrapper(std::map<std::string, std::map<std::string, int>> hd, Json::Value oc, Json::Value ec, std::unordered_set<std::string> v, std::unordered_set<std::string> o, std::string itype);
    std::map<std::string, int> get_hamming_distances(std::string kmer);
    int get_hamming_distance(std::string kmer1,std::string kmer2);
};