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

float get_expected_count(std::unordered_set<std::string> &Si, std::shared_ptr<KmersWrapper> kmer_wrap_ptr, std::string kmer, const Json::Value &expectedCounts, bool isExpectedKmer, float normalizer) {
    //Use expectation value if available
    BOOST_LOG_TRIVIAL(info) << "get_expected_count \n";
    float expectedCount = 0;
    if (isExpectedKmer){
        expectedCount = expectedCounts.get(kmer,-1).asFloat();
        if (expectedCount == -1){
            BOOST_LOG_TRIVIAL(fatal) << kmer << " was not found in the expected kmers but should be there, aborting! \n";
            throw std::exception();
        }
    }
    else {
        // expectedCount = expectedDefaultValue; //Default value for non-expected but observed kmers is the previously calculated value
        // iterrate through all spatype kmers add a^hd * (1-a)^(len-hd) when hd small enough
        expectedCount = 0.01;
        int foundmaxexpectedcount = 0;
        std::map<std::string, int> hd_kmer = (*kmer_wrap_ptr.get()).get_hamming_distances(kmer);
        if(hd_kmer.size() < Si.size()) {
            std::map<std::string, int>::iterator it;
            for ( it = hd_kmer.begin(); it != hd_kmer.end(); it++) {
                std::string target_kmer = it->first;
                // if kmer in Si    
                if(Si.find(target_kmer) != Si.end()) {    
                    int hd = it->second;                                                        // = hamming distance(kmer, target_kmer)
                    float e_i = expectedCounts.get(target_kmer,0).asFloat();                    // = |target_kmer|
                    float a_hd = ((*kmer_wrap_ptr.get()).get_computed_probability(hd));         // = a^hd * (1-a)^(len-hd)

                    expectedCount += a_hd*e_i*normalizer;
                    //int current = a_hd*e_i;
                    //if(foundmaxexpectedcount < current) {
                    //    foundmaxexpectedcount = current;
                    //}
                }
            }
        } else {
            for (std::unordered_set<std::string>::const_iterator sikmer = Si.begin(); sikmer != Si.end(); sikmer++){
                if ( hd_kmer.count(*sikmer) > 0 ) {
                    // a^hd * (1-a)^(len-hd) * |sikmer| when hd small enough
                    int hd = hd_kmer[*sikmer];                                              // = hamming distance(kmer, sikmer)
                    float e_i = expectedCounts.get(*sikmer,0).asFloat();                    // = |sikmer|
                    float a_hd = ((*kmer_wrap_ptr.get()).get_computed_probability(hd));     // = a^hd * (1-a)^(len-hd)

                    expectedCount += a_hd*e_i*normalizer;
                    // int current = a_hd*e_i;
                    // if(foundmaxexpectedcount < current) {
                    //     foundmaxexpectedcount = current;
                    // }
                }
            }
        }
        //expectedCount += foundmaxexpectedcount*normalizer;
        //std::cout << expectedCount << " result; " << normalizer << " normalizer;\n";
    }

    BOOST_LOG_TRIVIAL(info) << "get_expected_count END \n";
    return expectedCount;
}

int get_observation_count(std::string kmer, Json::Value &observedCounts) {
    BOOST_LOG_TRIVIAL(info) << "get_observation_count \n";
    //Use observation value
    int observedCount = observedCounts.get(kmer,-1).asInt();
    if (observedCount == -1){
        observedCount = 0;
        //BOOST_LOG_TRIVIAL(fatal) << *kmer << " was not found in the observed kmers but should be there, aborting! \n";
        //throw std::exception();
    }
    BOOST_LOG_TRIVIAL(info) << "get_observation_count END \n";
    return observedCount;
}

bool a_subset_of_b(std::unordered_set<std::string> &a, std::unordered_set<std::string> &b) {
    BOOST_LOG_TRIVIAL(info) << "a_subset_of_b \n";
    for(std::unordered_set<std::string>::const_iterator kmer=a.begin(); kmer!=a.end(); ++kmer) {
        if (b.find(*kmer) == b.end()){ //not found
            BOOST_LOG_TRIVIAL(fatal) << *kmer << " was not found in O but was in Si \n";
            BOOST_LOG_TRIVIAL(fatal) << *(b.begin()) << " example O entry \n";
            BOOST_LOG_TRIVIAL(fatal) << b.size() << " length O \n";
            BOOST_LOG_TRIVIAL(fatal) << a.size() << " length Si \n";
            return false;
        }
    }
    BOOST_LOG_TRIVIAL(info) << "a_subset_of_b END \n";
    return true;
}

std::unordered_set<std::string> relative_complement(std::unordered_set<std::string> &a, std::unordered_set<std::string> &b) {
    std::set<std::string> orda(a.begin(), a.end());
    std::set<std::string> ordb(b.begin(), b.end());
    std::unordered_set <std::string> anb;
    std::unordered_set <std::string> adiffanb;
    std::set_intersection(orda.begin(), orda.end(), ordb.begin(), ordb.end(), std::inserter(anb, anb.begin()));
    std::set_difference(orda.begin(), orda.end(), ordb.begin(), ordb.end(), std::inserter(adiffanb, adiffanb.begin()));
    return adiffanb;
}

probabilistic::CoverageBasedResult probabilistic::calculateLikelihoodCoverageBased(
            int threadID,
            const std::shared_ptr<KmersWrapper> kmer_wrap_ptr,
            const Json::Value expectedCounts,
            const float kmerError,
            const std::string spaTypeName,
            const int deviationCutoff
        ){
    //std::cout << "Deviation Cutoff: " << deviationCutoff << std::endl;

    // BOOST_LOG_TRIVIAL(info) << "kmerwrap: " << kmer_wrap_ptr.get() << ", kmer_wrap_ptr.get() \n";
    // BOOST_LOG_TRIVIAL(info) << "kmerwrap: " << &(*kmer_wrap_ptr.get()) << ", &(*kmer_wrap_ptr.get()) \n";
    // BOOST_LOG_TRIVIAL(info) << "kmerwrap: " << &(kmer_wrap) << ", &(kmer_wrap) \n";
    // BOOST_LOG_TRIVIAL(info) << "kmerwrap: " << &((*kmer_wrap_ptr.get()).hamming_distance_matrix) << ", &((*kmer_wrap_ptr.get()).hamming_distance_matrix) \n";
    // BOOST_LOG_TRIVIAL(info) << "kmerwrap: " << &(kmer_wrap.hamming_distance_matrix) << ", &(kmer_wrap.hamming_distance_matrix) \n";

    //Create empty result object
    probabilistic::CoverageBasedResult result;
    result.likelihood  = 0.0;
    result.errorLikelihood = 0.0;
    BOOST_LOG_TRIVIAL(info) << "starting calculateLikelihoodCoverageBased \n";
    // Si = expectedKmers
    std::unordered_set<std::string> Si;
    for(Json::Value::const_iterator kmer=expectedCounts.begin(); kmer!=expectedCounts.end(); ++kmer) {
        Si.insert(kmer.key().asString());
    }

    //Calculate set of assumed error kmers
    std::unordered_set<std::string> assumedErrorKmers;
    for (std::unordered_set<std::string>::const_iterator kmer = (*kmer_wrap_ptr.get()).O.begin(); kmer != (*kmer_wrap_ptr.get()).O.end(); kmer++)
    {
        if (Si.find(*kmer) != Si.end()){ //equiv to kmer is expected
            //we don't want this
        }
        else{
            assumedErrorKmers.insert(*kmer);
        }
    }

    // TODO: REMOVE??
    //Sanity Check: If an expected k-mer is not observed at all we discard this type instantly
    //if (not (a_subset_of_b(Si, (*kmer_wrap_ptr.get()).O))) { //not found
    //    result.likelihood = NAN;
    //    result.errorLikelihood = NAN;
    //    return result;
    //}
    // Ignore spatypes wich dont have enough kmers with O in common
    if (relative_complement(Si, (*kmer_wrap_ptr.get()).O).size()>3) { //not found
        result.likelihood = NAN;
        result.errorLikelihood = NAN;
        return result;
    }
    //Calculate default value for expected counts
    int sumOfObservedCounts = 0;
    for(Json::Value::const_iterator kmer=(*kmer_wrap_ptr.get()).observedCounts.begin(); kmer!=(*kmer_wrap_ptr.get()).observedCounts.end(); ++kmer) {
        sumOfObservedCounts += get_observation_count(kmer.key().asString(), (*kmer_wrap_ptr.get()).observedCounts);
    }

    //TODO: check if normalizer correct
    float old_normalizer =sumOfObservedCounts * kmerError / assumedErrorKmers.size();
    float default_normalizer = sumOfObservedCounts*kmerError/Si.size();
    float normalizer = 1;
    // |O|*e/|Uo|
    // float expectedDefaultValue = sumOfObservedCounts * kmerError / assumedErrorKmers.size();
    // BOOST_LOG_TRIVIAL(info) << spaTypeName << "\t" << expectedDefaultValue << "\n";
    unsigned int expectedErrors = sumOfObservedCounts * kmerError;
    unsigned int observedErrors = 0;

    BOOST_LOG_TRIVIAL(info) << "Beginning for loop \n";
    //Calculate likelihoods
    for(std::unordered_set<std::string>::const_iterator kmer=(*kmer_wrap_ptr.get()).iterset.begin(); kmer!=(*kmer_wrap_ptr.get()).iterset.end(); ++kmer) {

        bool isExpectedKmer = (Si.find(*kmer) != Si.end());
        if(!isExpectedKmer) {
            observedErrors += 1;
        }
        int observedCount = get_observation_count(*kmer, (*kmer_wrap_ptr.get()).observedCounts);
        float expectedCount = get_expected_count(Si, kmer_wrap_ptr, *kmer, expectedCounts, isExpectedKmer, normalizer);

        if (SINGLE_ERROR_TERM && !isExpectedKmer){
            //kmer is an error and not taken into account as a single term
        }
        else{
            //if the deviation cutoff is used we decide here whether or not to drop the spa type immediately
            if (deviationCutoff != -1){
                if (abs(observedCount-expectedCount) >= deviationCutoff){
                    BOOST_LOG_TRIVIAL(fatal) << observedCount << " observedCount," << expectedCount << " expectedCount, difference was to big in spatype " << spaTypeName << "\n";
                    result.likelihood = NAN;
                    result.errorLikelihood = NAN;
                    return result;
                }
            }

            if(expectedCount == 0) {
                // really unlikely meaning: log(pmf(0, 0.000001)) approx log(1) = 0
                continue;
            }
            //Assuming poisson distribution, thus: Pr(X = k) = (lambda^k * e^-lambda) / k!
            //long double lnLikelihood = log(std::pow(expectedCount,observedCount) * exp(-expectedCount) / boost::math::factorial<long double>(observedCount));
            BOOST_LOG_TRIVIAL(info) << "poisson_pmf \n";
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
