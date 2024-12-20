#include "math_def.hh"
#include "mesh.hh"
#include "mesher.hh"
#include "mesh.io.hh"
#include "main.hh"

using namespace OpenMesh;

int main(const int argc, const char **argv)
{
    int err {};
    std::string filename, prefix, path;
    if (argc < 2) { printf("No file provided.\n"); return 1; }

    double min_angle = 26.0;
    double max_length = 1e5;
    double max_area = .5e10;

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

    err = refine(mesh, radian(min_angle), max_length, max_area);
    save_mesh(mesh, (prefix + ".refined.mesh").c_str());
    if (err) { printf("Refinement failed.\n"); return err; }

    return err;
}