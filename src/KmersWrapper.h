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
    void pre_compute_hd_probabilities(float kmerError, int max_hd);
    float* l;

    public:
    std::map<std::string, std::map<std::string, int>> hamming_distance_matrix;
    std::map<std::string, std::map<std::string, int>> ohamming_distance_matrix;
    Json::Value observedCounts;
    Json::Value expectedCounts;
    std::unordered_set<std::string> U;
    std::unordered_set<std::string> V;
    std::unordered_set<std::string> O;
    std::unordered_set<std::string> iterset;
    std::string itersetType;

    KmersWrapper(std::string hammingdist, std::string kmersindex, std::string ukmersindex, std::string observed, std::string expected, std::string itype, float kmerError, std::string ohammingdist, std::string okmersindex);
    ~KmersWrapper();
    std::map<std::string, int> get_hamming_distances(std::string kmer);
    int get_hamming_distance(std::string kmer1,std::string kmer2);
    float get_computed_probability(int i);
};