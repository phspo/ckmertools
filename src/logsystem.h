
#ifndef C_KMERTOOLS_LOGSYSTEM_H
#define C_KMERTOOLS_LOGSYSTEM_H

#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/expressions.hpp>
#include <string>

namespace logging {
    void initLogging(std::string logFilePath);
};


#endif //C_KMERTOOLS_LOGSYSTEM_H
