// KmersWrapper.h

#include <iostream>
#include <map>
#include <string>
#include <iterator>

class KmersWrapper {
    public:
    std::map<std::string, std::map<std::string, int>> hamming_distance_matrix;
    
    std::map<std::string, int> get_hamming_distances(std::string kmer);
    int get_hamming_distance(std::string kmer1,std::string kmer2);
};