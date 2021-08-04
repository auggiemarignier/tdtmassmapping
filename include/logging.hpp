#include <string>

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define ERROR(logger, fmt, ...) logger.write_log(ERROR, fmt, ##__VA_ARGS__)
#define WARNING(logger, fmt, ...) logger.write_log(WARNING, fmt, ##__VA_ARGS__)
#define INFO(logger, fmt, ...) logger.write_log(INFO, fmt, ##__VA_ARGS__)
#define DEBUG(logger, fmt, ...) logger.write_log(DEBUG, fmt, ##__VA_ARGS__)

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

    const char *filename;
    FILE *log_file;
};