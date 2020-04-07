#ifndef C_KMERTOOLS_PARSING_H
#define C_KMERTOOLS_PARSING_H

#include <iostream>
#include <fstream>
#include "jsoncpp/json/json.h"
#include <boost/algorithm/string.hpp>
#include<jsoncpp/json/writer.h>

namespace parsing {
    Json::Value readDictionary(const std::string &filePath);
    void writeDictionary(const std::map<std::string,double> &map,const std::string &filePath);
};


#endif //C_KMERTOOLS_PARSING_H
