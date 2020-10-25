
#ifndef C_KMERTOOLS_PROBABILISTIC_H
#define C_KMERTOOLS_PROBABILISTIC_H


#include <unordered_set>
#include "jsoncpp/json/json.h"
#include <unordered_set>
#include <algorithm>
#include "logsystem.h"
#include <boost/math/special_functions/factorials.hpp>
#include <cmath>
#include <iostream>


namespace probabilistic {

    struct CoverageBasedResult{
        double likelihood;
        double errorLikelihood;
    };

    struct GenerativeResult{
        double likelihood;
    };

    double poisson_pmf(const double k, const double lambda);

    int hamming_distance(std::string s1, std::string s2);

    CoverageBasedResult calculateLikelihoodCoverageBased(
            int threadID,
            const std::shared_ptr<Json::Value> observedCountsPointer,
            const Json::Value &expectedCounts,
            const float &kmerError,
            const std::string spaTypeName,
            const int deviationCutoff,
            const std::shared_ptr<std::unordered_set<std::string>> OPointer,
            const std::shared_ptr<std::unordered_set<std::string>> itersetPointer
    );

    GenerativeResult calculateLikelihoodGenerative(
            int threadID,
            const std::shared_ptr<Json::Value> observedCountsPointer,
            const Json::Value &sequenceProfile,
            const float &baseErrorRate,
            const double* hdLikelihoods
            );
};


#endif //C_KMERTOOLS_PROBABILISTIC_H
