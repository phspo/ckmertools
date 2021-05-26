#include <iostream>
#include <memory>
#include "logsystem.h"
#include "parsing.h"
#include "probabilistic.h"

#include <boost/program_options.hpp>
#include "ctpl/ctpl.h"
#include <stdint.h>
// https://github.com/rogersce/cnpy


using namespace boost::program_options;

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
                ("unexpected,u",value<std::string>()->required(),"Unexpected kmers file")
                ("itersetType,i",value<std::string>()->required(),"Iterset Type O, V, OuV or OnV")
                ("hammingdist,h",value<std::string>()->required(),"path to the hammingdistance npz file")
                ("ukmersindex,j",value<std::string>(),"path to the kmers json file")
                ("vkmersindex,x",value<std::string>()->required(),"path to the kmers json file")
                ("ohammingdist,y",value<std::string>()->required(),"path to the hammingdistance npz file")
                ("okmersindex,z",value<std::string>()->required(),"path to the kmers json file")
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
        int mode = vm["m"].as<int>();
        //Coverave-Based Mode
        if (mode == 0){
            BOOST_LOG_TRIVIAL(info) << "Running coverage based \n";
            Json::Value expectedCounts = parsing::readDictionary(vm["expected"].as<std::string>());

            float kmerError = vm["kmererror"].as<float>();
            std::cout <<"Init KmersWrapper \n";
            std::string ukmersindex = "";
            if (vm.count("ukmersindex")) {
                ukmersindex = vm["ukmersindex"].as<std::string>();
            }
            std::shared_ptr<KmersWrapper> kmer_wrap_ptr = std::make_shared<KmersWrapper>(vm["hammingdist"].as<std::string>(),
            vm["vkmersindex"].as<std::string>(),
            ukmersindex,
            vm["observed"].as<std::string>(),
            vm["expected"].as<std::string>(), 
            vm["itersetType"].as<std::string>(),
            kmerError,
            vm["ohammingdist"].as<std::string>(),
            vm["okmersindex"].as<std::string>()
            );

            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << kmer_wrap_ptr.get() << ", kmer_wrap_ptr.get() \n";
            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << &(*kmer_wrap_ptr.get()) << ", &(*kmer_wrap_ptr.get()) \n";
            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << &(kmer_wrap) << ", &(kmer_wrap) \n";
            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << &((*kmer_wrap_ptr.get()).hamming_distance_matrix) << ", &((*kmer_wrap_ptr.get()).hamming_distance_matrix) \n";
            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << &(kmer_wrap.hamming_distance_matrix) << ", &(kmer_wrap.hamming_distance_matrix) \n";
            

            std::map<std::string,double> likelihoods;
            std::map<std::string,double> unexpectedKmerLikelihoods;

            BOOST_LOG_TRIVIAL(info) << "threadpools \n";
            ctpl::thread_pool p(vm.count("cores") ? vm["cores"].as<int>() : 1 );

            std::vector<std::future< probabilistic::CoverageBasedResult>> results(expectedCounts.size());
            // one thread variant:
            //std::vector<probabilistic::CoverageBasedResult> results(expectedCounts.size()); 
            int idx = 0;
            //Distribute tasks
            BOOST_LOG_TRIVIAL(info) << "disribute \n";
            for(Json::Value::const_iterator spaType=expectedCounts.begin(); spaType!=expectedCounts.end(); ++spaType, ++ idx) {
                if (spaType->getMemberNames().size() > 0){
                    int deviationCutoff = vm.count("deviationcutoff")  ? vm["deviationcutoff"].as<int>()  :  -1;
                    BOOST_LOG_TRIVIAL(info) << "Pushed " << spaType.key().asString() << " to threadpool \n";
                    results[idx] = p.push(probabilistic::calculateLikelihoodCoverageBased,kmer_wrap_ptr,*spaType,kmerError,spaType.key().asString(),deviationCutoff);
                    // one thread variant:
                    //results[idx] = probabilistic::calculateLikelihoodCoverageBased(1,kmer_wrap_ptr,*spaType,kmerError,spaType.key().asString(),deviationCutoff);
                    BOOST_LOG_TRIVIAL(info) << "Done " << spaType.key().asString() << " threadpool \n";
                }
                else{
                    BOOST_LOG_TRIVIAL(info) << "No expected k-mers found for spa-type: " << spaType.key().asString() << ", maybe the type is too small? \n";
                }
            }
            BOOST_LOG_TRIVIAL(info) << "DONE ALL \n";
            idx = 0;
            //Fetch results
            for(Json::Value::const_iterator spaType=expectedCounts.begin(); spaType!=expectedCounts.end(); ++spaType, ++idx) {
                if (spaType->getMemberNames().size() > 0){
                    probabilistic::CoverageBasedResult result = results[idx].get();
                    // one thread variant:
                    //probabilistic::CoverageBasedResult result = results[idx];
                    likelihoods.insert(std::pair<std::string,long double>(spaType.key().asString(),result.likelihood));
                    unexpectedKmerLikelihoods.insert(std::pair<std::string,long double>(spaType.key().asString(),result.likelihood));
                    BOOST_LOG_TRIVIAL(info) << "PARSED " << spaType.key().asString() << "\n";
                    //std::cout << result.likelihood << "/" << result.errorLikelihood << "\n";
                }
                else{
                    //Ignore, warning was given during job distribution phase
                }
            }
            BOOST_LOG_TRIVIAL(info) << "DONE PARSE \n";
            parsing::writeDictionary(likelihoods,vm["target"].as<std::string>());
            if (vm.count("unexpected")){
                parsing::writeDictionary(unexpectedKmerLikelihoods,vm["unexpected"].as<std::string>());
            }
            kmer_wrap_ptr.reset();
            BOOST_LOG_TRIVIAL(info) << "FINISHED \n";
            return 0;
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
        BOOST_LOG_TRIVIAL(info) << "ERROR DETECTED \n";
        BOOST_LOG_TRIVIAL(error) << ex.what() << '\n';
    }
    BOOST_LOG_TRIVIAL(info) << "FINISHED \n";
}
