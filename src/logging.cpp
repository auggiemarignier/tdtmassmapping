#include "logging.hpp"

Logger::Logger() : filename(nullptr), log_file(NULL){};

void Logger::open_log_file()
{
    if (filename !=  nullptr)
    {
        log_file = fopen(filename, "w");
        if (log_file == NULL)
        {
            fprintf(stderr, "Failed to open log file %s\n", filename);
        }
    }
    else
    {
        log_file = stderr;
    }
}
