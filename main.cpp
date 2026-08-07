#include "mesh_types.hpp"
#include "math_utils.hpp"
#include "bem_interop.hpp"
#include "bvh_collisions.hpp"
#include "boolean_ops.hpp"
#include "gmsh_interop.hpp"
#include <gmsh.h>
#include <gmsh.h>
#include <iostream>
#include <filesystem>
#include <string>
namespace fs = std::filesystem;
bool verifyMeshConformity(const CollisionContext& ctx, const MeshData& meshA, const MeshData& meshB) {
    auto checkMeshCuts = [&](const MeshData& mesh, 
                             const std::unordered_map<size_t, PolyLine>& NCcoords,
                             const std::unordered_map<size_t, std::vector<PolyLine>>& Ccoords,
                             const std::string& meshLabel) -> bool 
    {
        SpatialGrid3D grid(ctx.eps);
        for (const Vec3& node : mesh.nodes) {
            grid.getOrAdd(node);
        }

        auto isNodeInMesh = [&](const Vec3& pt) -> bool {
            double epsSq = ctx.eps * ctx.eps * 4.0;
            for (const Vec3& node : mesh.nodes) {
                Vec3 d = pt - node;
                if (dot(d, d) <= epsSq) {
                    return true;
                }
            }
            return false;
        };

        size_t totalCutPointsChecked = 0;
        size_t missingPoints = 0;

        auto checkPolyline = [&](const PolyLine& poly) {
            for (const auto& seg : poly) {
                totalCutPointsChecked++;
                if (!isNodeInMesh(seg.first)) {
                    //std::cout << "[" << meshLabel << " FAIL] Missing cut point: (" << seg.first[0] << ", " << seg.first[1] << ", " << seg.first[2] << ")\n";
                    missingPoints++;
                }
                Vec3 diff = seg.first - seg.second;
                if (dot(diff, diff) > ctx.eps * ctx.eps) {
                    totalCutPointsChecked++;
                    if (!isNodeInMesh(seg.second)) {
                        //std::cout << "[" << meshLabel << " FAIL] Missing cut point: (" << seg.second[0] << ", " << seg.second[1] << ", " << seg.second[2] << ")\n";
                        missingPoints++;
                    }
                }
            }
        };

        for (const auto& [_, cuts] : NCcoords) checkPolyline(cuts);
        for (const auto& [_, polyList] : Ccoords) {
            for (const auto& poly : polyList) checkPolyline(poly);
        }

        if (missingPoints > 0) {
            std::cout << "[" << meshLabel << " VERIFIER] Failed: " << missingPoints << " / " 
                      << totalCutPointsChecked << " cut points missing!\n";
            return false;
        }

        std::cout << "[" << meshLabel << " VERIFIER] Success: All " << totalCutPointsChecked 
                  << " cut endpoints map to nodes!\n";
        return true;
    };

    bool okA = checkMeshCuts(meshA, ctx.NCAcoords, ctx.CAcoords, "Mesh A");
    bool okB = checkMeshCuts(meshB, ctx.NCBcoords, ctx.CBcoords, "Mesh B");
    return okA && okB;
}

static std::string joinMsh(const fs::path& folder, const char* name) {
    fs::path p = folder;
    p /= name;
    return p.string();
}


void testMeshDifference(bem::TriangleMesh<3>& meshA,
                         bem::TriangleMesh<3>& meshB,
                         const std::string& folder,
                         bool cleanDegenerate = false)
{
    // Work on copies so tests don't affect each other
    auto A = meshA;
    auto B = meshB;

    CollisionContext originalCtx =
        detectCollisions(extractMeshData(A), extractMeshData(B));

    // difference1 / difference2
    meshDifference(A, B, cleanDegenerate);

    saveMSH(A, joinMsh(folder, "difference1.msh"));
    saveMSH(B, joinMsh(folder, "difference2.msh"));

    MeshData finalMeshA = extractMeshData(A);
    MeshData finalMeshB = extractMeshData(B);
    verifyMeshConformity(originalCtx, finalMeshA, finalMeshB);
}

void testMeshUnion(bem::TriangleMesh<3>& meshA,
                    bem::TriangleMesh<3>& meshB,
                    const std::string& folder,
                    bool cleanDegenerate = false)
{
    auto A = meshA;
    auto B = meshB;

    CollisionContext originalCtx =
        detectCollisions(extractMeshData(A), extractMeshData(B));

    // union1 / union2
    meshUnion(A, B, cleanDegenerate);

    saveMSH(A, joinMsh(folder, "union1.msh"));
    saveMSH(B, joinMsh(folder, "union2.msh"));

    MeshData finalMeshA = extractMeshData(A);
    MeshData finalMeshB = extractMeshData(B);
    verifyMeshConformity(originalCtx, finalMeshA, finalMeshB);
}

void testMeshIntersect(bem::TriangleMesh<3>& meshA,
                        bem::TriangleMesh<3>& meshB,
                        const std::string& folder,
                        bool cleanDegenerate = false)
{
    auto A = meshA;
    auto B = meshB;

    CollisionContext originalCtx =
        detectCollisions(extractMeshData(A), extractMeshData(B));

    // intersect1 / intersect2
    meshIntersect(A, B, cleanDegenerate);

    saveMSH(A, joinMsh(folder, "intersect1.msh"));
    saveMSH(B, joinMsh(folder, "intersect2.msh"));

    MeshData finalMeshA = extractMeshData(A);
    MeshData finalMeshB = extractMeshData(B);
    verifyMeshConformity(originalCtx, finalMeshA, finalMeshB);
}

void testMeshCombine(bem::TriangleMesh<3>& meshA,
                      bem::TriangleMesh<3>& meshB,
                      const std::string& folder)
{
    auto A = meshA;
    auto B = meshB;

    CollisionContext originalCtx =
        detectCollisions(extractMeshData(A), extractMeshData(B));

    // combine1 / combine2
    meshCombine(A, B, false, false);

    saveMSH(A, joinMsh(folder, "combine1.msh"));
    saveMSH(B, joinMsh(folder, "combine2.msh"));

    MeshData finalMeshA = extractMeshData(A);
    MeshData finalMeshB = extractMeshData(B);
    verifyMeshConformity(originalCtx, finalMeshA, finalMeshB);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <inFolder> <outFolder>\n";
        return 1;
    }
    std::string inFolder  = argv[1];
    std::string outFolder = argv[2];
    std::string mshAPath = inFolder + "/untitled.msh";
    std::string mshBPath = inFolder + "/untitled2.msh";
    bem::TriangleMesh<3> meshA = loadMSH(mshAPath);
    bem::TriangleMesh<3> meshB = loadMSH(mshBPath);
    testMeshDifference(meshA, meshB, outFolder);
    testMeshUnion(meshA, meshB, outFolder);
    testMeshIntersect(meshA, meshB, outFolder);
    testMeshCombine(meshA, meshB, outFolder);
    return 0;
}

