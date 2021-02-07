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
        int mode = vm["m"].as<int>();
        //Coverave-Based Mode
        if (mode == 0){
            BOOST_LOG_TRIVIAL(info) << "Running coverage based \n";
            Json::Value expectedCounts = parsing::readDictionary(vm["expected"].as<std::string>());

            float kmerError = vm["kmererror"].as<float>();
            
            std::shared_ptr<KmersWrapper> kmer_wrap_ptr = std::make_shared<KmersWrapper>(vm["hammingdist"].as<std::string>(),
            vm["kmersindex"].as<std::string>(),
            vm["observed"].as<std::string>(),
            vm["expected"].as<std::string>(), 
            vm["itersetType"].as<std::string>(),
            kmerError);

            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << kmer_wrap_ptr.get() << ", kmer_wrap_ptr.get() \n";
            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << &(*kmer_wrap_ptr.get()) << ", &(*kmer_wrap_ptr.get()) \n";
            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << &(kmer_wrap) << ", &(kmer_wrap) \n";
            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << &((*kmer_wrap_ptr.get()).hamming_distance_matrix) << ", &((*kmer_wrap_ptr.get()).hamming_distance_matrix) \n";
            // BOOST_LOG_TRIVIAL(info) << "INITIAL PTR for kmerwrap: " << &(kmer_wrap.hamming_distance_matrix) << ", &(kmer_wrap.hamming_distance_matrix) \n";
            

            std::map<std::string,double> likelihoods;
            std::map<std::string,double> unexpectedKmerLikelihoods;


            //ctpl::thread_pool p(vm.count("cores") ? vm["cores"].as<int>() : 1 );

            //std::vector<std::future< probabilistic::CoverageBasedResult>> results(expectedCounts.size());
            std::vector<probabilistic::CoverageBasedResult> results(expectedCounts.size()); 
            int idx = 0;
            //Distribute tasks
            for(Json::Value::const_iterator spaType=expectedCounts.begin(); spaType!=expectedCounts.end(); ++spaType, ++ idx) {
                if (spaType->getMemberNames().size() > 0){
                    int deviationCutoff = vm.count("deviationcutoff")  ? vm["deviationcutoff"].as<int>()  :  -1;
                    BOOST_LOG_TRIVIAL(info) << "Pushed " << spaType.key().asString() << " to threadpool \n";
                    //results[idx] = p.push(probabilistic::calculateLikelihoodCoverageBased,kmer_wrap_ptr,*spaType,kmerError,spaType.key().asString(),deviationCutoff);
                    results[idx] = probabilistic::calculateLikelihoodCoverageBased(1,kmer_wrap_ptr,*spaType,kmerError,spaType.key().asString(),deviationCutoff);
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
                    BOOST_LOG_TRIVIAL(info) << "PARSING result \n";
                    probabilistic::CoverageBasedResult result = results[idx]; //.get();
                    BOOST_LOG_TRIVIAL(info) << "PARSED result \n";
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
            // stop from deallocing kmer_wrap_ptr?
            BOOST_LOG_TRIVIAL(info) << "DONE IF \n";
            // return early?
            BOOST_LOG_TRIVIAL(info) << "RESET \n";
            kmer_wrap_ptr.reset();
            BOOST_LOG_TRIVIAL(info) << "FINISHED \n";
            return 0;
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
