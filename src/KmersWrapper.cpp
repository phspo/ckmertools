#include "KmersWrapper.h"

std::map<std::string, int> KmersWrapper::get_hamming_distances(std::string kmer) {
    return KmersWrapper::hamming_distance_matrix[kmer];
}

int KmersWrapper::get_hamming_distance(std::string kmer1,std::string kmer2){
    return KmersWrapper::hamming_distance_matrix[kmer1][kmer2];
}
