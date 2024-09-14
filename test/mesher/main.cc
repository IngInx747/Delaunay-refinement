#include "math_def.hh"
#include "mesh.hh"
#include "mesher.hh"
#include "mesh.io.hh"

using namespace OpenMesh;

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

int main(const int argc, const char **argv)
{
    int err {};
    std::string filename, prefix, path;
    if (argc < 2) { printf("No file provided.\n"); return 1; }

    double min_angle = 26.0;
    double max_length = 0;
    double max_area = 0;

    const char *arg {};
    if ((arg = get_value(argv, argv + argc, "--min-angle")))  { min_angle = atof(arg); }
    if ((arg = get_value(argv, argv + argc, "--max-length"))) { max_length = atof(arg); }
    if ((arg = get_value(argv, argv + argc, "--max-area")))   { max_area = atof(arg); }

    filename.append(argv[1]);
    prefix = filename.substr(0, filename.find_last_of("."));
    path = filename.substr(0, filename.find_last_of("/\\"));

    std::vector<Vec2> vs {};
    std::vector<Int2> es {};
    std::vector<Vec2> ss {};
    if (read_poly(vs, es, ss, filename.c_str()) != 0)
    { printf("Cannot load poly file.\n"); return 1; }

    TriMesh mesh;
    getOrMakeProperty<Mh, std::string>(mesh, var_m_name())() = prefix;
    getOrMakeProperty<Mh, std::string>(mesh, var_m_path())() = path;

    err = triangulate(vs, es, mesh);
    mesh.delete_isolated_vertices(); // remove dups
    save_mesh(mesh, (prefix + ".CDT.mesh").c_str());
    if (err) { printf("Segments intersecting.\n"); return err; }

    hide_exterior_region(mesh, ss);
    save_mesh(mesh, (prefix + ".CDT.mesh").c_str());

    err = refine(mesh, radian(min_angle));
    save_mesh(mesh, (prefix + ".refined.mesh").c_str());
    if (err) { printf("Refinement failed.\n"); return err; }

    return err;
}