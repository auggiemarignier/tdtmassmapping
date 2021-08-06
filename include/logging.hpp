#include <string>
#include <time.h>

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define ERROR(fmt, ...) Logger::Get().write_log(ERROR, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define WARNING(fmt, ...) Logger::Get().write_log(WARNING, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define INFO(fmt, ...) Logger::Get().write_log(INFO, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define DEBUG(fmt, ...) Logger::Get().write_log(DEBUG, __FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

typedef enum
{
    ERROR = 0,
    WARNING = 1,
    INFO = 3,
    DEBUG = 4
} logger_level_t;

class Logger
{
public:
    Logger(const Logger &) = delete;
    Logger operator=(const Logger &) = delete;

    static Logger &Get();

    void write_log(logger_level_t level, const char *sourcefile, const char *function, int lineno, const char *fmt, ...);
    void open_log_file(const char *filename);
    void close_log();

private:
    Logger(){};
    FILE *log_file;
};