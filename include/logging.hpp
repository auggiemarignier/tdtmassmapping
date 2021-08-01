#include <string>

extern "C"
{
#include "slog.h"
};

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


#define ERROR(fmt, ...) slog(SLOG_ERROR, __FILENAME__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define WARNING(fmt, ...) slog(SLOG_WARNING, "", "", 0, fmt, ##__VA_ARGS__)
#define INFO(fmt, ...) slog(SLOG_INFO, "", "", 0, fmt, ##__VA_ARGS__)
#define DEBUG(fmt, ...) slog(SLOG_DEBUG, __FILENAME__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

typedef enum
{
    ERROR = 0,
    WARNING = 1,
    INFO = 3,
    DEBUG = 4
} logger_level_t;

struct Logger
{
    Logger();

    void open_log_file();
    void write_log(logger_level_t level, const char *fmt, ...);
    void close_log();

    void log_error(const char *fmt, ...);
    void log_warning(const char *fmt, ...);
    void log_info(const char *fmt, ...);
    void log_debug(const char *fmt, ...);

    const char *filename;
    FILE *log_file;
};