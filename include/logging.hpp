#pragma once

#include <string>
#include <time.h>
#include <stdarg.h>

#ifdef ERROR
#undef ERROR
#endif
#ifdef WARNING
#undef WARNING
#endif
#ifdef INFO
#undef INFO
#endif
#ifdef DEBUG
#undef DEBUG
#endif

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define ERROR(fmt, ...) Logger::write_log(ERROR, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define WARNING(fmt, ...) Logger::write_log(WARNING, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define INFO(fmt, ...) Logger::write_log(INFO, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define DEBUG(fmt, ...) Logger::write_log(DEBUG, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG(fmt, ...) Logger::write_log(BLANK, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

typedef enum
{
    ERROR = 0,
    WARNING = 1,
    INFO = 3,
    DEBUG = 4,
    BLANK = 5
} logger_level_t;

class Logger
{
public:
    Logger(const Logger &) = delete;
    Logger operator=(const Logger &) = delete;

    static Logger &Get();

    ~Logger();

    static void *write_log(logger_level_t level, const char *sourcefile, const char *function, int lineno, const char *fmt, ...)
    {
        va_list args;
        va_start(args, fmt);
        Logger::Get().write_log_imp(level, sourcefile, function, lineno, fmt, args);
        va_end(args);
    }

    static void open_log(const char *filename)
    {
        Logger::Get().open_log_imp(filename);
    }

    static void close_log()
    {
        Logger::Get().close_log_imp();
    }

private:
    Logger(){};
    void open_log_imp(const char *filename);
    void close_log_imp();

    void *write_log_imp(logger_level_t level, const char *sourcefile, const char *function, int lineno, const char *fmt, va_list arg);

    FILE *log_file;
    bool log_open = false;
};