#include "logging.hpp"

#include <stdarg.h>
#include <iostream>

static char *timestamp();


Logger& Logger::Get()
{
    static Logger instance;
    return instance;
}

void Logger::open_log_file()
{
    if (filename != nullptr)
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
        fprintf(log_file, "WARNING: ");
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

static char *timestamp()
{
    const char *TIME_FORMAT = "%Y-%m-%d %H:%M:%S";
    static char buffer[32];
    time_t tmp;
    struct tm *t;

    tmp = time(NULL);
    t = localtime(&tmp);
    strftime(buffer, sizeof(buffer), TIME_FORMAT, t);

    return buffer;
}
