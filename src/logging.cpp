#include "logging.hpp"

#include <stdarg.h>
#include <string>

static char *timestamp();

Logger &Logger::Get()
{
    static Logger instance;
    return instance;
}

Logger::~Logger()
{
    close_log();
    log_open = false;
}

void Logger::open_log_imp(std::string filename)
{
    const char *c = const_cast<char *>(filename.c_str());
    open_log_imp(c);
}

void Logger::open_log_imp(const char *filename)
{
    if (filename != nullptr)
    {
        log_file = fopen(filename, "a");
        if (log_file == NULL)
        {
            fprintf(stderr, "Failed to open log file %s\n", filename);
        }
    }
    else
        log_file = stderr;
    log_open = true;
    fprintf(log_file, "Starting logs %s\n\n", timestamp());
}

void *Logger::write_log_imp(logger_level_t level, const char *sourcefile, const char *function, int lineno, const char *fmt, va_list args)
{
    if (!log_open) // Nothing ever gets logged
        return nullptr;

    switch (level)
    {
    case ERROR:
        fprintf(log_file, "ERROR:\t%s[%d] ", sourcefile, lineno);
        break;
    case WARNING:
        fprintf(log_file, "WARNING: ");
        break;
    case INFO:
        fprintf(log_file, "INFO:\t");
        break;
    case DEBUG:
        fprintf(log_file, "DEBUG:\t%s:%s[%d] ", sourcefile, function, lineno);
        break;
    default:
        break;
    }

    vfprintf(log_file, fmt, args);
    fprintf(log_file, "\n");
    return nullptr;
}

void Logger::close_log_imp()
{
    if (log_open)
    {
        fprintf(log_file, "\nClosing logs %s\n", timestamp());
        if (log_file != stderr)
            fclose(log_file);
    }
}

void Logger::flush_imp()
{
    fflush(log_file);
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
