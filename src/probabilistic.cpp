#include <iostream>
#include "probabilistic.h"
#include<cmath>

//const bool SKIP_ERRORS = true;
const bool SINGLE_ERROR_TERM = false;

double probabilistic::poisson_pmf(const double k, const double lambda)
{
    return k * log(lambda) - lgamma(k + 1.0) - lambda;
}

int probabilistic::hamming_distance(std::string s1, std::string s2){
    int count = 0;
    for (int i=0; i < s1.length(); i++){
        count += (s1[i] != s2[i]);
    }
    return count;
}

float * pre_compute_hd_probabilities(const float &normalizer, const int &length_kmer, const float kmerError) {
    int max_hd = 5;
    float *a = (float*) malloc(sizeof(float) * max_hd);
    a[0] = 1;
    for (int i = 1; i < max_hd; i++)
    {
        //TODO: replace with mutation rate
        // i is hamming distance
        a[i] = normalizer*std::pow(kmerError,i)*std::pow((1 - kmerError),(length_kmer - i));
    }
    return a;
}

float get_expected_count(std::unordered_set<std::string> &Si, std::map<std::string, std::map<std::string, int>> &hd, float a[], std::string kmer, const Json::Value &expectedCounts, bool isExpectedKmer) {
    //Use expectation value if available
    float expectedCount = 0;
    if (isExpectedKmer){
        expectedCount = expectedCounts.get(kmer,-1).asFloat();
        if (expectedCount == -1){
            BOOST_LOG_TRIVIAL(fatal) << kmer << " was not found in the expected kmers but should be there, aborting! \n";
            throw std::exception();
        }
    } else{
        // expectedCount = expectedDefaultValue; //Default value for non-expected but observed kmers is the previously calculated value
        // iterrate through all spatype kmers add a^hd * (1-a)^(len-hd) when hd small enough
        
        if(hd[kmer].size() < Si.size()) {
            std::map<std::string, int>::iterator it;
            for ( it = hd[kmer].begin(); it != hd[kmer].end(); it++) {
                // if kmer in Si
                std::string target_kmer = it->first;
                int hd = it->second;
                expectedCount += a[hd]*expectedCounts.get(target_kmer,0).asFloat();
            }
        } else {
            for (std::unordered_set<std::string>::const_iterator sikmer = Si.begin(); sikmer != Si.end(); sikmer++){
                if ( hd[kmer].count(*sikmer) > 0 ) {
                    // a^hd * (1-a)^(len-hd) * |sikmer| when hd small enough
                    expectedCount += a[hd[kmer][*sikmer]]*expectedCounts.get(*sikmer,0).asFloat();
                }
            }
        }
    }
    return expectedCount;
}

int get_observation_count(std::string kmer, Json::Value &observedCounts) {
    //Use observation value
    int observedCount = observedCounts.get(kmer,-1).asInt();
    if (observedCount == -1){
        observedCount = 0;
        //BOOST_LOG_TRIVIAL(fatal) << *kmer << " was not found in the observed kmers but should be there, aborting! \n";
        //throw std::exception();
    }
    return observedCount;
}

bool a_subset_of_b(std::unordered_set<std::string> &a, std::unordered_set<std::string> &b) {
    for(std::unordered_set<std::string>::const_iterator kmer=a.begin(); kmer!=a.end(); ++kmer) {
        if (b.find(*kmer) == b.end()){ //not found
            return false;
        }
    }
    return true;
}

probabilistic::CoverageBasedResult error_result(probabilistic::CoverageBasedResult result) {
    result.likelihood = NAN;
    result.errorLikelihood = NAN;
    return result;
}


probabilistic::CoverageBasedResult probabilistic::calculateLikelihoodCoverageBased(
        const int threadID,
        const std::shared_ptr<Json::Value> observedCountsPointer,
        const Json::Value &expectedCounts,
        const float &kmerError,
        const std::string spaTypeName,
        const int deviationCutoff,
        const std::shared_ptr<std::unordered_set<std::string>> OPointer,
        const std::shared_ptr<std::unordered_set<std::string>> itersetPointer,
        std::shared_ptr<std::map<std::string, std::map<std::string, int>>> hammingDistancedPointer
        ){

    //std::cout << "Pointer: " << observedCountsPointer << std::endl;
    Json::Value observedCounts = *observedCountsPointer.get();

    std::map<std::string, std::map<std::string, int>> hd = *hammingDistancedPointer.get();
    //std::cout << "Deviation Cutoff: " << deviationCutoff << std::endl;

    //Create empty result object
    probabilistic::CoverageBasedResult result;
    result.likelihood  = 0.0;
    result.errorLikelihood = 0.0;

    // O = observedKmers
    std::unordered_set<std::string> O = *OPointer;
    // Iterationset = O or V or both?
    std::unordered_set<std::string> iterset = *itersetPointer.get();

    // Si = expectedKmers
    std::unordered_set<std::string> Si;
    for(Json::Value::const_iterator kmer=expectedCounts.begin(); kmer!=expectedCounts.end(); ++kmer) {
        Si.insert(kmer.key().asString());
    }

    //Calculate set of assumed error kmers
    std::unordered_set<std::string> assumedErrorKmers;
    for (std::unordered_set<std::string>::const_iterator kmer = O.begin(); kmer != O.end(); kmer++)
    {
        if (Si.find(*kmer) != Si.end()){ //equiv to kmer is expected
            //we don't want this
        }
        else{
            assumedErrorKmers.insert(*kmer);
        }
    }

    // TODO: WHAT HAPPEND HERE? REPLACE O WITH itersetPointer ?
    //Sanity Check: If an expected k-mer is not observed at all we discard this type instantly
    if (!a_subset_of_b(Si, O)) { //not found
        return error_result(result);
    }

    //Calculate default value for expected counts
    int sumOfObservedCounts = 0;
    for(Json::Value::const_iterator kmer=observedCounts.begin(); kmer!=observedCounts.end(); ++kmer) {
        sumOfObservedCounts += get_observation_count(kmer.key().asString(), observedCounts);
    }

    //TODO: check if normalizer correct
    float normalizer = sumOfObservedCounts*kmerError/ Si.size();
    // kmer length
    int length_kmer = (*begin(Si)).length();
    float* a = pre_compute_hd_probabilities(normalizer, length_kmer, kmerError);

    // |O|*e/|Uo|
    // float expectedDefaultValue = sumOfObservedCounts * kmerError / assumedErrorKmers.size();
    // BOOST_LOG_TRIVIAL(info) << spaTypeName << "\t" << expectedDefaultValue << "\n";
    unsigned int expectedErrors = sumOfObservedCounts * kmerError;
    unsigned int observedErrors = 0;

    //Calculate likelihoods
    for(std::unordered_set<std::string>::const_iterator kmer=iterset.begin(); kmer!=iterset.end(); ++kmer) {

        bool isExpectedKmer = (Si.find(*kmer) != Si.end());
        if(!isExpectedKmer) {
            observedErrors += 1;
        }
        int observedCount = get_observation_count(*kmer, observedCounts);
        float expectedCount = get_expected_count(Si, hd, a, *kmer, expectedCounts, isExpectedKmer);

        if (SINGLE_ERROR_TERM && !isExpectedKmer){
            //kmer is an error and not taken into account as a single term
        }
        else{
            //if the deviation cutoff is used we decide here whether or not to drop the spa type immediately
            if (deviationCutoff != -1){
                if (abs(observedCount-expectedCount) >= deviationCutoff){
                    return error_result(result);
                }
            }

            if(expectedCount == 0) {
                // really unlikely meaning: log(pmf(0, 0.000001)) approx log(1) = 0
                continue;
            }
            //Assuming poisson distribution, thus: Pr(X = k) = (lambda^k * e^-lambda) / k!
            //long double lnLikelihood = log(std::pow(expectedCount,observedCount) * exp(-expectedCount) / boost::math::factorial<long double>(observedCount));
            double lnLikelihood = poisson_pmf(observedCount,expectedCount);

            if (lnLikelihood == 0 or lnLikelihood != lnLikelihood){
                BOOST_LOG_TRIVIAL(fatal) << "lnLikelihood of 0 or NaN was calculated ... (" << lnLikelihood << ") \n" << "kmerID: " << *kmer <<" observed: " << observedCount << " expected: " << expectedCount << "\n";
                throw std::exception();
            }
            else{
                //Add value to the total (because we are in log-space we add to multiply the probabilities)
                if (!isExpectedKmer){ //if this is not an expected kmer, we also add it to the probability of the unexpected kmers
                    result.errorLikelihood += lnLikelihood;
                }
                result.likelihood += lnLikelihood;
                if (result.likelihood == 0){
                    BOOST_LOG_TRIVIAL(fatal) << "lnLikelihood of 0 resulted after addition of " << lnLikelihood << " (" << result.likelihood << ") \n";
                    throw std::exception();
                }
            }
        }

    }

    if (SINGLE_ERROR_TERM){
        double lnLikelihood = poisson_pmf(observedErrors,expectedErrors);
        result.likelihood += lnLikelihood;
        result.errorLikelihood += lnLikelihood;
    }

    return result;
}


probabilistic::GenerativeResult probabilistic::calculateLikelihoodGenerative(
        const int threadID,
        const std::shared_ptr<Json::Value> observedCountsPointer,
        const Json::Value &sequenceProfile,
        const float &baseErrorRate,
        const double* hdLikelihoods)
         {

    Json::Value observedCounts = *observedCountsPointer.get();

    GenerativeResult result;

    double nrOfSequenceKmers = 0; //Must be double to perform real division

    for (Json::Value::const_iterator kmer = sequenceProfile.begin(); kmer != sequenceProfile.end(); ++kmer) {
        nrOfSequenceKmers += kmer->asInt();
    }

    double chanceOfDrawingKmer = 1 / nrOfSequenceKmers;

    //Switch to log space
    double lnLikelihood = 0;

    for (Json::Value::const_iterator kmerObs = observedCounts.begin(); kmerObs != observedCounts.end(); ++kmerObs) {
        long double observationProbability = 0;
        int observationCount = kmerObs->asInt();
        for (Json::Value::const_iterator kmerGen = sequenceProfile.begin(); kmerGen != sequenceProfile.end(); ++kmerGen) {
            int factor = kmerGen->asInt();
            int hammingDistance = probabilistic::hamming_distance(kmerObs.key().asString(),kmerGen.key().asString());

            double hdLikelihood = hdLikelihoods[hammingDistance];
            observationProbability += factor*hdLikelihood;
            //std::cout << observationProbability << "/" << hdLikelihood << "/" << factor << "\n";
        }
        observationProbability *= chanceOfDrawingKmer;
        observationProbability = pow(observationProbability,observationCount);
        lnLikelihood += log(observationProbability);
        if (lnLikelihood >= 0 or lnLikelihood != lnLikelihood){
            BOOST_LOG_TRIVIAL(fatal) << "LnLikelihood of > 0 or NaN was calculated ... Precision insufficient? [lnLikelihood: " << lnLikelihood << " obsProbability: " <<observationProbability << "]\n";
            throw std::exception();
        }
    }


    result.likelihood = lnLikelihood;
    return result;
}
