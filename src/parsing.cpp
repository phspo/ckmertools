#include "parsing.h"

Json::Value parsing::readDictionary(const std::string &filePath) {
    std::ifstream ifs(filePath);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    ifs.close();
    return obj;
}

void parsing::writeDictionary(const std::map<std::string,double> &map,const std::string &filePath){
    Json::Value output;

    for (auto it=map.begin(); it!=map.end(); ++it){
        output[it->first] = Json::Value(it->second);
    }


    std::ofstream ofs(filePath);
    Json::StyledWriter writer;
    ofs << writer.write(output);
    ofs.close();
}
