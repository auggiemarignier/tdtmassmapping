#include "logging.hpp"

#include <stdarg.h>

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

void Logger::log_error(const char *fmt, ...)
{
    fprintf(log_file, "ERROR:\t%s[%d] ", __FILE__, __LINE__);

    va_list args;
    va_start(args, fmt);
    vfprintf(log_file, fmt, args);
    va_end(args);

    fprintf(log_file, "\n");
}

void Logger::log_warning(const char *fmt, ...)
{
    fprintf(log_file, "WARNING:\t");

    va_list args;
    va_start(args, fmt);
    vfprintf(log_file, fmt, args);
    va_end(args);

    fprintf(log_file, "\n");
}

void Logger::log_info(const char *fmt, ...)
{
    fprintf(log_file, "INFO:\t");

    va_list args;
    va_start(args, fmt);
    vfprintf(log_file, fmt, args);
    va_end(args);

    fprintf(log_file, "\n");
}

void Logger::log_debug(const char *fmt, ...)
{
    fprintf(log_file, "DEBUG:\t%s:%s:[%d] ", __FILE__, __FUNCTION__, __LINE__);

    va_list args;
    va_start(args, fmt);
    vfprintf(log_file, fmt, args);
    va_end(args);

    fprintf(log_file, "\n");
}