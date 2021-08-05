#include "logging.hpp"

#include <stdarg.h>
#include <time.h>

Logger::Logger() : filename(nullptr), log_file(stderr){};

Logger::Logger(const char *filename) : filename(filename), log_file(NULL)
{
    open_log_file();
    write_log(INFO, "Beginning logs\n");
};

Logger::~Logger()
{
    write_log(INFO, "Ending logs\n");
}

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

void Logger::write_log(logger_level_t level, const char *fmt, ...)
{
    if (level == ERROR)
        fprintf(log_file, "ERROR:\t%s[%d] ", __FILE__, __LINE__);
    else if (level == WARNING)
        fprintf(log_file, "WARNING:\t");
    else if (level == INFO)
        fprintf(log_file, "INFO:\t");
    else
        fprintf(log_file, "DEBUG:\t%s:%s[%d] ", __FILE__, __FUNCTION__, __LINE__);

    va_list args;
    va_start(args, fmt);
    vfprintf(log_file, fmt, args);
    va_end(args);

    fprintf(log_file, "\n");
}

