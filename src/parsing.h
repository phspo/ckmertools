#ifndef C_KMERTOOLS_PARSING_H
#define C_KMERTOOLS_PARSING_H

#include <iostream>
#include <fstream>
#include "jsoncpp/json/json.h"
#include <boost/algorithm/string.hpp>
#include<jsoncpp/json/writer.h>
#include <unordered_set>

#include "cnpy.h"

namespace parsing {
    Json::Value readDictionary(const std::string &filePath);
    void writeDictionary(const std::map<std::string,double> &map,const std::string &filePath);
    std::map<std::string, std::map<std::string, int>> get_hammingdistances(std::string distances_file, std::string kmers_idx_file);
    std::map<std::string, std::map<std::string, int>> parsing::get_hammingdistancesO(std::string distances_file, std::string kmers_idx_fileO, std::string kmers_idx_fileV);
    std::unordered_set<std::string> get_V(Json::Value &expectedCounts);
    std::unordered_set<std::string> get_O(Json::Value &observedCounts);
};


#endif //C_KMERTOOLS_PARSING_H
