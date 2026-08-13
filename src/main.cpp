#include "boolean_ops.hpp"
#include "gmsh_interop.hpp"
#include <iostream>
#include <filesystem>
#include <string>
namespace fs = std::filesystem;

static std::string joinMsh(const fs::path& folder, const char* name) {
    fs::path p = folder;
    p /= name;
    return p.string();
}


void testMeshDifference(bem::TriangleMesh<3>& meshA,
                         bem::TriangleMesh<3>& meshB,
                         const std::string& folder)
{
    auto A = meshA;
    auto B = meshB;
    meshDifference(A, B);
    saveMSH(A, joinMsh(folder, "difference1.msh"));
    saveMSH(B, joinMsh(folder, "difference2.msh"));
}

void testMeshUnion(bem::TriangleMesh<3>& meshA,
                    bem::TriangleMesh<3>& meshB,
                    const std::string& folder)
{
    auto A = meshA;
    auto B = meshB;
    meshUnion(A, B);
    saveMSH(A, joinMsh(folder, "union1.msh"));
    saveMSH(B, joinMsh(folder, "union2.msh"));
}

void testMeshIntersect(bem::TriangleMesh<3>& meshA,
                        bem::TriangleMesh<3>& meshB,
                        const std::string& folder)
{
    auto A = meshA;
    auto B = meshB;
    meshIntersect(A, B);
    saveMSH(A, joinMsh(folder, "intersect1.msh"));
    saveMSH(B, joinMsh(folder, "intersect2.msh"));
}

void testMeshCombine(bem::TriangleMesh<3>& meshA,
                      bem::TriangleMesh<3>& meshB,
                      const std::string& folder,
                    bool removeTouchingSurfaces = false)
{
    auto A = meshA;
    auto B = meshB;
    meshCombine(A, B, removeTouchingSurfaces);
    saveMSH(A, joinMsh(folder, "combine1.msh"));
    saveMSH(B, joinMsh(folder, "combine2.msh"));
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <inFolder> <outFolder>\n";
        return 1;
    }
    std::string inFolder  = argv[1];
    std::string outFolder = argv[2];
    std::string mshAPath = inFolder + "/untitled1.msh";
    std::string mshBPath = inFolder + "/untitled2.msh";
    bem::TriangleMesh<3> meshA = loadMSH(mshAPath);
    bem::TriangleMesh<3> meshB = loadMSH(mshBPath);
    testMeshDifference(meshA, meshB, outFolder);
    testMeshUnion(meshA, meshB, outFolder);
    testMeshIntersect(meshA, meshB, outFolder);
    //testMeshCombine(meshA, meshB, outFolder, true);
    return 0;
}

