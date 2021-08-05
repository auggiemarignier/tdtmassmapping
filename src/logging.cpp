#include "logging.hpp"

#include <stdarg.h>

static char *timestamp();

Logger::Logger() : filename(nullptr), log_file(stderr)
{
    write_log(INFO, "Beginning logs %s", timestamp());
};

Logger::Logger(const char *filename) : filename(filename), log_file(NULL)
{
    open_log_file();
    write_log(INFO, "Beginning logs %s", timestamp());
};

Logger::~Logger()
{
    write_log(INFO, "Ending logs %s", timestamp());
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
