#include <iostream>
#include <memory>
#include "logsystem.h"
#include "parsing.h"
#include "probabilistic.h"
#include <boost/program_options.hpp>
#include <unordered_set>
#include "ctpl/ctpl.h"
#include <stdint.h>
#include "cnpy.h"
// https://github.com/rogersce/cnpy


using namespace boost::program_options;


enum ItersetOptions {
    InvalidType,
    OType,
    VType,
    OuVType,
    OnVType
};

ItersetOptions resolveItersetOption(std::string input) {
    if( input == "O" ) return OType;
    if( input == "V" ) return VType;
    if( input == "OuV" ) return OuVType;
    if( input == "OnV" ) return OnVType;
    return InvalidType;
}


std::shared_ptr<std::unordered_set<std::string>> chooseIterset(std::shared_ptr<std::unordered_set<std::string>> OPointer, Json::Value expectedCounts, std::string itersetType) {
    //Define sets of all Spatype kmers
    std::unordered_set <std::string> V;
    for (Json::Value::const_iterator spaType = expectedCounts.begin(); spaType != expectedCounts.end(); ++spaType) {
        if (spaType->getMemberNames().size() > 0) {
            for (Json::Value::const_iterator kmer = (*spaType).begin(); kmer != (*spaType).end(); ++kmer) {
                V.insert(kmer.key().asString());
            }
        }
    }
    auto VPointer = std::make_shared<std::unordered_set<std::string>>(V);
    std::unordered_set<std::string> O = *OPointer;

    /////////////////////////////////////// TEST SETS
    std::set<std::string> ordO(O.begin(), O.end());
    std::set<std::string> ordV(V.begin(), V.end());
    // O Union V
    std::unordered_set <std::string> OuV;
    std::set_union(ordO.begin(), ordO.end(), ordV.begin(), ordV.end(), std::inserter(OuV, OuV.begin()));

    // O Intersect V
    std::unordered_set <std::string> OnV;
    std::set_intersection(ordO.begin(), ordO.end(), ordV.begin(), ordV.end(), std::inserter(OnV, OnV.begin()));

    //std::cout << "Length O: " << O.size() << std::endl;
    //std::cout << "Length V: " << V.size() << std::endl;
    //std::cout << "Length OuV: " << OuV.size() << std::endl;
    //std::cout << "Length OnV: " << OnV.size() << std::endl;
    /////////////////////////////////////////////


    switch( resolveItersetOption(itersetType) )
    {
        case OType: {
            return OPointer;
        }
        case VType: {
            return VPointer;
        }
        case OuVType: {
            auto iterset = std::make_shared<std::unordered_set<std::string>>(OuV);
            return iterset;
        }
        case OnVType: {
            auto iterset = std::make_shared<std::unordered_set<std::string>>(OnV);
            return iterset;
        }
        default: {
            return OPointer;
        }
    }
}

// distances_file = "V_kmer_distances.npz";     kmers_idx_file="V_kmers.json"
std::shared_ptr<std::map<std::string, std::map<std::string, int>>> get_hammingdistances(std::string distances_file, std::string kmers_idx_file) {
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
    auto hd = std::make_shared<std::map<std::string, std::map<std::string, int>>>(hamming_distances);
    return hd;
} 

// e.g. ./c_kmertools --e "expected_counts.json" --c "1" --m 0 --o "alignment.counts.json" --kmererror 0.1 --d 1 --target "likelihoods_cov.json" --unexpected "unexpected_likelihoods_cov.json" --log "likelihoods_cov.log" --itersetType "O" --hammingdist "V_kmer_distances.npz" --kmersindex "V_kmers.json"
 int main(int argc, char* argv[]) {
    try
    {
        options_description desc{"Options"};
        desc.add_options()
                ("help,h", "Help screen")
                ("log,l", value<std::string>(), "Log file")
                ("expected,e",value<std::string>(),"Expected counts json file")
                ("profiles,p",value<std::string>(),"Spa type kmer profiles json file")
                ("observed,o",value<std::string>()->required(),"Observed counts json file")
                ("kmererror,r",value<float>(),"Kmer error rate")
                ("baseerror,b",value<float>(),"Base error rate")
                ("k",value<int>(),"k for kmers")
                ("cores,c",value<int>(),"Number of cores that the software may use")
                ("deviationcutoff,d",value<int>(),"Discard any type with kmer differences exceeding the cutoff value")
                ("m",value<int>()->required(),"Mode, [0] for coverage-based, [1] for generative")
                ("target,t",value<std::string>()->required(),"Target file")
                ("unexpected,u",value<std::string>(),"Unexpected kmers file")
                ("itersetType,i",value<std::string>(),"Iterset Type O, V, OuV or OnV")
                ("hammingdist,h",value<std::string>(),"path to the hammingdistance npz file")
                ("kmersindex,j",value<std::string>(),"path to the kmers json file")
        ;
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        if (vm.count("help"))
            std::cout << desc << '\n';
        else{
            if (vm.count("log")){
                logging::initLogging(vm["log"].as<std::string>());
                BOOST_LOG_TRIVIAL(info) << "Log system initialized ...";
            }
        }

        //Coverave-Based Mode
        if (vm["m"].as<int>() == 0){

            Json::Value expectedCounts = parsing::readDictionary(vm["expected"].as<std::string>());
            Json::Value observedCounts = parsing::readDictionary(vm["observed"].as<std::string>());
            auto observedCountsPointer = std::make_shared<Json::Value>(observedCounts);
            std::shared_ptr<std::map<std::string, std::map<std::string, int>>> hamming_distance = get_hammingdistances(vm["hammingdist"].as<std::string>(), vm["kmersindex"].as<std::string>());

            float kmerError = vm["kmererror"].as<float>();

            std::map<std::string,double> likelihoods;
            std::map<std::string,double> unexpectedKmerLikelihoods;


            ctpl::thread_pool p(vm.count("cores") ? vm["cores"].as<int>() : 1 );

            std::vector<std::future< probabilistic::CoverageBasedResult>> results(expectedCounts.size());


            /// HANDLE SETS /////////////////////
            std::unordered_set <std::string> O;
            for (Json::Value::const_iterator kmer = observedCounts.begin(); kmer != observedCounts.end(); ++kmer) {
                O.insert(kmer.key().asString());
            }
            std::string itersetType = vm["itersetType"].as<std::string>();
            auto OPointer = std::make_shared<std::unordered_set<std::string>>(O);
            std::shared_ptr<std::unordered_set<std::string>> itersetPointer = chooseIterset(OPointer,expectedCounts, itersetType);
            ////////////////////////////


            int idx = 0;
            //Distribute tasks
            for(Json::Value::const_iterator spaType=expectedCounts.begin(); spaType!=expectedCounts.end(); ++spaType, ++ idx) {
                if (spaType->getMemberNames().size() > 0){
                    int deviationCutoff =  vm.count("deviationcutoff")  ? vm["deviationcutoff"].as<int>()  :  -1;
                    results[idx] = p.push(probabilistic::calculateLikelihoodCoverageBased,observedCountsPointer,*spaType,kmerError,spaType.key().asString(),deviationCutoff,OPointer,itersetPointer);
                }
                else{
                    BOOST_LOG_TRIVIAL(info) << "No expected k-mers found for spa-type: " << spaType.key().asString() << ", maybe the type is too small? \n";
                }
            }
            idx = 0;
            //Fetch results
            for(Json::Value::const_iterator spaType=expectedCounts.begin(); spaType!=expectedCounts.end(); ++spaType, ++idx) {
                if (spaType->getMemberNames().size() > 0){
                    probabilistic::CoverageBasedResult result = results[idx].get();
                    likelihoods.insert(std::pair<std::string,long double>(spaType.key().asString(),result.likelihood));
                    unexpectedKmerLikelihoods.insert(std::pair<std::string,long double>(spaType.key().asString(),result.likelihood));
                    //std::cout << result.likelihood << "/" << result.errorLikelihood << "\n";
                }
                else{
                    //Ignore, warning was given during job distribution phase
                }
            }
            parsing::writeDictionary(likelihoods,vm["target"].as<std::string>());
            if (vm.count("unexpected")){
                parsing::writeDictionary(unexpectedKmerLikelihoods,vm["unexpected"].as<std::string>());
            }
        }
        else if (vm["m"].as<int>() == 1) {
            Json::Value sequenceProfiles = parsing::readDictionary(vm["profiles"].as<std::string>());
            Json::Value observedCounts = parsing::readDictionary(vm["observed"].as<std::string>());
            
            auto observedCountsPointer = std::make_shared<Json::Value>(observedCounts);

            std::map<std::string,double> likelihoods;

            float baseErrorRate = vm["baseerror"].as<float>();
            int k = vm["k"].as<int>();

            //Pre-Calculate all Hamming-Distance probabilities
            double* hdLikelihoods =  new double[k+1];
            for (int i=0;i<k+1;i++){
                double intactBasesProbability = pow((1-baseErrorRate),(k-i));
                double alteredBasesProbability = pow(baseErrorRate,i);
                hdLikelihoods[i] = intactBasesProbability * alteredBasesProbability;
            }

            ctpl::thread_pool p(vm.count("cores") ? vm["cores"].as<int>() : 1 );

            std::vector<std::future< probabilistic::GenerativeResult>> results(sequenceProfiles.size());

            int idx = 0;
            //Distribute tasks
            for(Json::Value::const_iterator spaType=sequenceProfiles.begin(); spaType!=sequenceProfiles.end(); ++spaType,++idx) {
                //std::cout << "Calculating spa Type: " << spaType.key().asString()+ "\n";
                if (spaType->getMemberNames().size() > 0){
                    //probabilistic::GenerativeResult result = probabilistic::calculateLikelihoodGenerative(idx,observedCounts,*spaType,baseErrorRate,hdLikelihoods);
                    results[idx] = p.push(probabilistic::calculateLikelihoodGenerative,observedCountsPointer,*spaType,baseErrorRate,hdLikelihoods);
                }
                else{
                    BOOST_LOG_TRIVIAL(info) << "No expected k-mers found for spa-type: " << spaType.key().asString() << ", maybe the type is too small? \n";
                }
            }
            idx = 0;
            //Collect results
            for(Json::Value::const_iterator spaType=sequenceProfiles.begin(); spaType!=sequenceProfiles.end(); ++spaType, ++idx) {
                //std::cout << "Calculating spa Type: " << spaType.key().asString()+ "\n";
                if (spaType->getMemberNames().size() > 0){
                    probabilistic::GenerativeResult result = results[idx].get();
                    likelihoods.insert(std::pair<std::string,long double>(spaType.key().asString(),result.likelihood));
                    //std::cout << result.likelihood << "\n";
                }
                else{
                    //Handled before
                }
            }

            delete hdLikelihoods;

            parsing::writeDictionary(likelihoods,vm["target"].as<std::string>());
        }
        else{
            BOOST_LOG_TRIVIAL(error) << "Illegal Mode ... \n";
        }
    }
    catch (const error &ex)
    {
        BOOST_LOG_TRIVIAL(error) << ex.what() << '\n';
    }
}
