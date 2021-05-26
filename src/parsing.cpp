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

// distances_file = "V_kmer_distances.npz";     kmers_idx_file="V_kmers.json"
std::map<std::string, std::map<std::string, int>> parsing::get_hammingdistances(std::string distances_file, std::string kmers_idx_file) {
    cnpy::npz_t M = cnpy::npz_load(distances_file);
    std::vector<uint8_t> M_data = M["data"].as_vec<uint8_t>();
    std::vector<int> M_i = M["col"].as_vec<int>();
    std::vector<int> M_j = M["row"].as_vec<int>();
    std::vector<int> M_shape = M["shape"].as_vec<int>();

    Json::Value V_kmers_index_json = parsing::readDictionary(kmers_idx_file);
    std::vector<std::string> V_kmers_index;
    int idx = 0;
    for(Json::Value::const_iterator kmer=V_kmers_index_json.begin(); kmer!=V_kmers_index_json.end(); ++kmer, ++ idx ) {
        V_kmers_index.push_back(kmer->asString());
    }
    std::map<std::string, std::map<std::string, int>> hamming_distances;
    for (int i = 0; i < M_data.size(); i++) {
        std::string id_x = V_kmers_index[M_i[i]];
        std::string id_y = V_kmers_index[M_j[i]];
        hamming_distances[id_x][id_y] = M_data[i];
    }
    // example for hd
    //std::cout << "TEST:" << std::to_string(hamming_distances["TTTTTGCCAGGCTTGTTGTTGTCTTCTTTACCAGGCTT"]["TTTTTGCCAGGCTTGTTATTGTCTTCTTTGCCAGGCTT"]) << std::endl;
    return hamming_distances;
}

// distances_file = "V_kmer_distances.npz";     kmers_idx_file="V_kmers.json"
std::map<std::string, std::map<std::string, int>> parsing::get_hammingdistancesO(std::string distances_file, std::string kmers_idx_fileO, std::string kmers_idx_fileV) {
    cnpy::npz_t M = cnpy::npz_load(distances_file);
    std::vector<uint8_t> M_data = M["data"].as_vec<uint8_t>();
    std::vector<int> M_i = M["col"].as_vec<int>();
    std::vector<int> M_j = M["row"].as_vec<int>();
    std::vector<int> M_shape = M["shape"].as_vec<int>();
    std::cout <<'begin indexing \n';
    Json::Value V_kmers_index_json = parsing::readDictionary(kmers_idx_fileV);
    Json::Value O_kmers_index_json = parsing::readDictionary(kmers_idx_fileO);
    std::vector<std::string> V_kmers_index;
    int idx = 0;
    for(Json::Value::const_iterator kmer=V_kmers_index_json.begin(); kmer!=V_kmers_index_json.end(); ++kmer, ++ idx ) {
        V_kmers_index.push_back(kmer->asString());
    }
    std::vector<std::string> O_kmers_index;
    idx = 0;
    for(Json::Value::const_iterator kmer=O_kmers_index_json.begin(); kmer!=O_kmers_index_json.end(); ++kmer, ++ idx ) {
        O_kmers_index.push_back(kmer->asString());
    }
    std::cout <<'begin hd write \n';
    std::map<std::string, std::map<std::string, int>> hamming_distances;
    for (int i = 0; i < M_data.size(); i++) {
        std::string id_x = O_kmers_index[M_i[i]];
        std::string id_y = V_kmers_index[M_j[i]];
        hamming_distances[id_x][id_y] = M_data[i];
    }
    std::cout <<'done hd write \n';
    // example for hd
    //std::cout << "TEST:" << std::to_string(hamming_distances["TTTTTGCCAGGCTTGTTGTTGTCTTCTTTACCAGGCTT"]["TTTTTGCCAGGCTTGTTATTGTCTTCTTTGCCAGGCTT"]) << std::endl;
    return hamming_distances;
}

std::unordered_set<std::string> parsing::get_V(Json::Value &expectedCounts) {
    std::unordered_set<std::string> V;
    for (Json::Value::const_iterator spaType = expectedCounts.begin(); spaType != expectedCounts.end(); ++spaType) {
        if (spaType->getMemberNames().size() > 0) {
            for (Json::Value::const_iterator kmer = (*spaType).begin(); kmer != (*spaType).end(); ++kmer) {
                V.insert(kmer.key().asString());
            }
        }
    }
    return V;
}

std::unordered_set<std::string> parsing::get_O(Json::Value &observedCounts) {
    std::unordered_set<std::string> O;
    for (Json::Value::const_iterator kmer = observedCounts.begin(); kmer != observedCounts.end(); ++kmer) {
        O.insert(kmer.key().asString());
    }
    return O;
}
