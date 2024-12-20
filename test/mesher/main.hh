#include <vector>
#include <string>

static inline bool is_key(const char *arg)
{
    if (*arg != '-') return false;
    ++arg; // eat 1st dash
    if (*arg == '\0') return false;
    if (isdigit(*arg)) return false;
    if (isalpha(*arg)) return true;
    ++arg; // eat 2nd dash
    return isalnum(*arg);
}

static inline bool has_key(const char **begin, const char **end, const char *key)
{
    return std::find(begin, end, std::string { key }) != end;
}

static inline const char **find_key_next(const char **begin, const char **end, const char *key)
{
    const char **iter = std::find(begin, end, std::string { key });
    if (iter != end && ++iter != end) return iter;
    return end;
}

static inline const char **find_next_key(const char **begin, const char **end)
{
    for (const char **iter = begin; iter != end; ++iter)
        if (is_key(*iter)) return iter;
    return end;
}

static inline bool has_value(const char **begin, const char **end, const char *key)
{
    const char **iter = std::find(begin, end, std::string { key });
    return iter != end && ++iter != end && !is_key(*iter);
}

static inline const char *get_value(const char **begin, const char **end, const char *key)
{
    const char **iter = std::find(begin, end, std::string { key });
    if (iter != end && ++iter != end) return *iter;
    return nullptr;
}

static std::vector<const char*> get_values(const int argc, const char **argv, const char *key)
{
    const char **begin = argv, **end = argv + argc;
    begin = find_key_next(begin, end, key);
    end   = find_next_key(begin, end);
    std::vector<const char*> vals {};
    for (const char **iter = begin; iter != end; ++iter)
        vals.push_back(*iter);
    return vals;
}
