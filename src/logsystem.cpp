
#include "logsystem.h"

void logging::initLogging(std::string logFilePath){
    boost::log::add_file_log(logFilePath,boost::log::keywords::auto_flush = true);

    boost::log::core::get()->set_filter
            (
                    boost::log::trivial::severity >= boost::log::trivial::info
            );
}
