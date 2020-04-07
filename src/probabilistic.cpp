
#include <iostream>
#include "probabilistic.h"

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


probabilistic::CoverageBasedResult probabilistic::calculateLikelihoodCoverageBased(
        const int threadID,
        const Json::Value &observedCounts,
        const Json::Value &expectedCounts,
        const float &kmerError,
        const std::string spaTypeName,
        const int deviationCutoff
        ){

    //std::cout << "Deviation Cutoff: " << deviationCutoff << std::endl;

    //Create empty result object
    probabilistic::CoverageBasedResult result;
    result.likelihood  = 0.0;
    result.errorLikelihood = 0.0;

    //Define sets of observed and expected kmers
    std::unordered_set<std::string> observedKmers;
    for(Json::Value::const_iterator kmer=observedCounts.begin(); kmer!=observedCounts.end(); ++kmer) {
        observedKmers.insert(kmer.key().asString());
    }

    std::unordered_set<std::string> expectedKmers;
    for(Json::Value::const_iterator kmer=expectedCounts.begin(); kmer!=expectedCounts.end(); ++kmer) {
        expectedKmers.insert(kmer.key().asString());
    }

    //Calculate set of assumed error kmers
    std::unordered_set<std::string> assumedErrorKmers;
    for (std::unordered_set<std::string>::const_iterator kmer = observedKmers.begin(); kmer != observedKmers.end(); kmer++)
    {
        if (expectedKmers.find(*kmer) != expectedKmers.end()){ //equiv to kmer is expected
            //we don't want this
        }
        else{
            assumedErrorKmers.insert(*kmer);
        }
    }


    //Sanity Check: If an expected k-mer is not observed at all we discard this type instantly
    for(std::unordered_set<std::string>::const_iterator kmer=expectedKmers.begin(); kmer!=expectedKmers.end(); ++kmer) {
        if (observedKmers.find(*kmer) == observedKmers.end()){ //not found
            result.likelihood = NAN;
            result.errorLikelihood = NAN;
            return result;
        }
    }

    //Calculate default value for expected counts
    int sumOfObservedCounts = 0;
    for(Json::Value::const_iterator kmer=observedCounts.begin(); kmer!=observedCounts.end(); ++kmer) {
        sumOfObservedCounts += kmer->asInt();
    }

    float expectedDefaultValue = sumOfObservedCounts * kmerError / assumedErrorKmers.size();

    unsigned int expectedErrors = sumOfObservedCounts * kmerError;

    unsigned int observedErrors = 0;

    BOOST_LOG_TRIVIAL(info) << spaTypeName << "\t" << expectedDefaultValue << "\n";

    //Calculate likelihoods
    for(std::unordered_set<std::string>::const_iterator kmer=observedKmers.begin(); kmer!=observedKmers.end(); ++kmer) {



        float expectedCount = expectedDefaultValue; //Default value for non-expected but observed kmers is the previously calculated value


        //Use observation value
        int observedCount = observedCounts.get(*kmer,-1).asInt();
        if (observedCount == -1){
            BOOST_LOG_TRIVIAL(fatal) << *kmer << " was not found in the observed kmers but should be there, aborting! \n";
            throw std::exception();
        }

        bool isExpectedKmer = false;

        //Use expectation value if available
        if (expectedKmers.find(*kmer) != expectedKmers.end()){
            expectedCount = expectedCounts.get(*kmer,-1).asFloat();
            isExpectedKmer = true;
            if (expectedCount == -1){
                BOOST_LOG_TRIVIAL(fatal) << *kmer << " was not found in the expected kmers but should be there, aborting! \n";
                throw std::exception();
            }
        }
        else{
            observedErrors += 1;
        }

        if (SINGLE_ERROR_TERM && !isExpectedKmer){
            //kmer is an error and not taken into account as a single term
        }
        else{
            //if the deviation cutoff is used we decide here whether or not to drop the spa type immediately
            if (deviationCutoff != -1){
                if (abs(observedCount-expectedCount) >= deviationCutoff){
                    result.likelihood = NAN;
                    result.errorLikelihood = NAN;
                    return result;
                }
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
        const Json::Value &observedCounts,
        const Json::Value &sequenceProfile,
        const float &baseErrorRate,
        const double* hdLikelihoods)
         {


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
