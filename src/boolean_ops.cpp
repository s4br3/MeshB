#include "boolean_ops.hpp"
#include "bem_interop.hpp"
#include "bvh_collisions.hpp"
#include "mesh_types.hpp"
#include "triangulation.hpp"
#include "classify.hpp"
#include "mesh_clean.hpp"
#include <cstddef>
#include <queue>
#include <utility>
std::vector<bool> getRemovalMask(const MeshData& mesh, const MeshData& targetMesh, const BVH& targetBVH, double eps, BoolOp op) {
    const size_t N = mesh.triangles.size();
    std::vector<bool> removeMask(N, false);
    if (N == 0) return removeMask;
    auto getEdgeKey = [eps](const Vec3& a, const Vec3& b) {
        Vec3 mid;
        mid.x = (a.x + b.x) * 0.5;
        mid.y = (a.y + b.y) * 0.5;
        mid.z = (a.z + b.z) * 0.5;
        return std::make_tuple(
            (long long)std::round(mid.x / eps),
            (long long)std::round(mid.y / eps),
            (long long)std::round(mid.z / eps)
        );
    };
    std::unordered_set<std::tuple<long long, long long, long long>, TupleHash> targetEdges;
    for (const auto& tri : targetMesh.triangles) {
        for (int j = 0; j < 3; ++j) {
            targetEdges.insert(getEdgeKey(
                targetMesh.nodes[tri.v[j]], 
                targetMesh.nodes[tri.v[(j + 1) % 3]]
            ));
        }
    }
    std::unordered_map<std::tuple<long long, long long, long long>, std::vector<size_t>, TupleHash> edgeToTris;
    for (size_t i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j) {
            edgeToTris[getEdgeKey(
                mesh.nodes[mesh.triangles[i].v[j]], 
                mesh.nodes[mesh.triangles[i].v[(j + 1) % 3]]
            )].push_back(i);
        }
    }
    std::vector<GWNNodeData> gwnData = buildGWNData(targetMesh, targetBVH);
    std::vector<bool> visited(N, false);
    std::vector<FaceClass> faceClasses(N, FaceClass::Outside);
    bool targetIsClosed = isMeshClosed(targetMesh);

    if (targetIsClosed) {
        for (size_t i = 0; i < N; ++i) {
            if (visited[i]) continue;
            double wn = evaluateGWN(mesh.centres[i], targetMesh, targetBVH, gwnData);
            FaceClass currentClass = (std::abs(wn) > 0.5) ? FaceClass::Inside : FaceClass::Outside;
            if (std::abs(std::abs(wn) - 0.5) < 0.1) {
                FaceClass exact = classifyFace(targetBVH, targetMesh, mesh.centres[i], mesh.normals[i], eps);
                if (exact == FaceClass::CoplanarSame || exact == FaceClass::CoplanarOpp) {
                    currentClass = exact;
                }
            }
            std::queue<size_t> q;
            q.push(i);
            visited[i] = true;
            faceClasses[i] = currentClass;
            while (!q.empty()) {
                size_t curr = q.front();
                q.pop();
                for (int j = 0; j < 3; ++j) {
                    auto edgeKey = getEdgeKey(
                        mesh.nodes[mesh.triangles[curr].v[j]], 
                        mesh.nodes[mesh.triangles[curr].v[(j + 1) % 3]]
                    );
                    if (targetEdges.find(edgeKey) != targetEdges.end()) {
                        continue;
                    }
                    for (size_t neighbor : edgeToTris[edgeKey]) {
                        if (!visited[neighbor]) {
                            visited[neighbor] = true;
                            faceClasses[neighbor] = currentClass;
                            q.push(neighbor);
                        }
                    }
                }
            }
        }
    }
    else{
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < N; ++i) {
            double wn = evaluateGWN(mesh.centres[i], targetMesh, targetBVH, gwnData);
            FaceClass currentClass = (std::abs(wn) > 0.5) ? FaceClass::Inside : FaceClass::Outside;
            if (std::abs(std::abs(wn) - 0.5) < 0.1) {
                FaceClass exact = classifyFace(targetBVH, targetMesh, mesh.centres[i], mesh.normals[i], eps);
                if (exact == FaceClass::CoplanarSame || exact == FaceClass::CoplanarOpp) {
                    currentClass = exact;
                }
            }
            faceClasses[i] = currentClass;
        }
    }
    for (size_t i = 0; i < N; ++i) {
        bool keep = false;
        switch (op) {
            case BoolOp::Union:
                keep = (faceClasses[i] == FaceClass::Outside || faceClasses[i] == FaceClass::CoplanarSame);
                break;
            case BoolOp::Intersect:
                keep = (faceClasses[i] == FaceClass::Inside || faceClasses[i] == FaceClass::CoplanarSame || faceClasses[i] == FaceClass::CoplanarOpp);
                break;
            case BoolOp::DifferenceBase:
                keep = (faceClasses[i] == FaceClass::Outside);
                break;
            case BoolOp::DifferenceTool:
                keep = (faceClasses[i] == FaceClass::Inside);
                break;
        }
        removeMask[i] = !keep;
    }
    return removeMask;
}
double computePolyLineArea(const PolyLine& poly) {
    Vec3 areaVec(0.0, 0.0, 0.0);
    for (const auto& seg : poly) {
        Vec3 c = cross(seg.first, seg.second);
        areaVec.x += c.x;
        areaVec.y += c.y;
        areaVec.z += c.z;
    }
    return 0.5 * std::sqrt(dot(areaVec, areaVec));
}
Connection nonConformal(const bem::TriangleMesh<3>& A, const bem::TriangleMesh<3>& B) {
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    CollisionContext ctx = detectCollisions(meshA, meshB);
    Connection conn;
    auto collectPoints = [](
        const std::unordered_map<size_t, PolyLine>& NCcoords,
        const std::unordered_map<size_t, std::vector<PolyLine>>& Ccoords,
        const MeshData& srcMesh,
        std::unordered_map<size_t, std::vector<Vec3>>& outIntersections)
    {
        for (const auto& [idx, poly] : NCcoords) {
            size_t tag = srcMesh.tags[idx];
            auto& pts = outIntersections[tag];
            for (const auto& seg : poly) {
                pts.push_back(seg.first);
                pts.push_back(seg.second);
            }
        }
        for (const auto& [idx, poly_list] : Ccoords) {
            size_t tag = srcMesh.tags[idx];
            auto& pts = outIntersections[tag];
            for (const auto& poly : poly_list) {
                for (const auto& seg : poly) {
                    pts.push_back(seg.first);
                    pts.push_back(seg.second);
                }
            }
        }
    };
    collectPoints(ctx.NCAcoords, ctx.CAcoords, ctx.meshDataA, conn.meshAIntersections);
    collectPoints(ctx.NCBcoords, ctx.CBcoords, ctx.meshDataB, conn.meshBIntersections);
    for (auto& [tag, pts] : conn.meshAIntersections) {
        SpatialGrid3D grid(ctx.eps);
        for (const auto& pt : pts) {
            grid.getOrAdd(pt);
        }
        pts = std::move(grid.getUniquePoints());
    }
    for (auto& [tag, pts] : conn.meshBIntersections) {
        SpatialGrid3D grid(ctx.eps);
        for (const auto& pt : pts) {
            grid.getOrAdd(pt);
        }
        pts = std::move(grid.getUniquePoints());
    }
    std::unordered_map<size_t, std::unordered_set<size_t>> aToBSet;
    std::unordered_map<size_t, std::unordered_set<size_t>> bToASet;
    const double areaEps = ctx.eps * ctx.eps;
    for (const auto& [idxA, targetTrisB] : ctx.Atris) {
        size_t tagA = ctx.meshDataA.tags[idxA];
        bool isCoplanar = (ctx.CAcoords.find(idxA) != ctx.CAcoords.end());
        for (size_t idxB : targetTrisB) {
            size_t tagB = ctx.meshDataB.tags[idxB];
            if (isCoplanar) {
                double area = 0.0;
                const auto& polyList = ctx.CAcoords.at(idxA);
                for (const auto& poly : polyList) {
                    area += computePolyLineArea(poly);
                }
                if (area > areaEps) {
                    conn.aToBAreas[tagA][tagB] += area;
                    aToBSet[tagA].insert(tagB);
                    bToASet[tagB].insert(tagA);
                }
            } else {
                aToBSet[tagA].insert(tagB);
                bToASet[tagB].insert(tagA);
            }
        }
    }
    for (const auto& [tagA, bSet] : aToBSet) {
        conn.aToBConnections[tagA].assign(bSet.begin(), bSet.end());
    }
    for (const auto& [tagB, aSet] : bToASet) {
        conn.bToAConnections[tagB].assign(aSet.begin(), aSet.end());
    }
    return conn;
}
CollisionContext collideAndCut(MeshData& meshA, MeshData& meshB, bool removeTouchingSurfaces) {
    CollisionContext ctx = detectCollisions(meshA, meshB);
    MeshData newMeshA = cutMesh(ctx.meshDataA, ctx.NCAcoords, ctx.CAcoords, ctx.eps, removeTouchingSurfaces);
    MeshData newMeshB = cutMesh(ctx.meshDataB, ctx.NCBcoords, ctx.CBcoords, ctx.eps, removeTouchingSurfaces);
    ctx.meshDataA = std::move(newMeshA);
    ctx.meshDataB = std::move(newMeshB);
    return ctx;
}
double meshCombine(MeshData& meshA, MeshData& meshB, bool removeTouchingSurfaces) {
    CollisionContext ctx = collideAndCut(meshA, meshB, removeTouchingSurfaces);
    meshA = std::move(ctx.meshDataA);
    meshB = std::move(ctx.meshDataB);
    return ctx.eps;
}
MeshData meshUnion(MeshData& meshA, MeshData& meshB, double* outEps) {
    CollisionContext ctx = collideAndCut(meshA, meshB, false);
    BVH bvhA = buildMeshBVH(ctx.meshDataA, ctx.eps);
    BVH bvhB = buildMeshBVH(ctx.meshDataB, ctx.eps);
    std::vector<bool> removeA = getRemovalMask(ctx.meshDataA, ctx.meshDataB, bvhB, ctx.eps, BoolOp::Union);
    std::vector<bool> removeB = getRemovalMask(ctx.meshDataB, ctx.meshDataA, bvhA, ctx.eps, BoolOp::DifferenceBase);
    filterMesh(ctx.meshDataA, removeA);
    filterMesh(ctx.meshDataB, removeB);
    if (outEps) *outEps = ctx.eps;
    return combineMeshes({ctx.meshDataA, ctx.meshDataB}, ctx.eps);
}
MeshData meshIntersect(MeshData& meshA, MeshData& meshB, double* outEps) {
        CollisionContext ctx = collideAndCut(meshA, meshB, false);
    BVH bvhA = buildMeshBVH(ctx.meshDataA, ctx.eps);
    BVH bvhB = buildMeshBVH(ctx.meshDataB, ctx.eps);
    std::vector<bool> removeA = getRemovalMask(ctx.meshDataA, ctx.meshDataB, bvhB, ctx.eps, BoolOp::Intersect);
    std::vector<bool> removeB = getRemovalMask(ctx.meshDataB, ctx.meshDataA, bvhA, ctx.eps, BoolOp::DifferenceTool);
    filterMesh(ctx.meshDataA, removeA);
    filterMesh(ctx.meshDataB, removeB);
    if (outEps) *outEps = ctx.eps;
    return combineMeshes({ctx.meshDataA, ctx.meshDataB}, ctx.eps);
}
double meshDifference(MeshData& meshA, MeshData& meshB) {
    CollisionContext ctx = collideAndCut(meshA, meshB, false);
    BVH bvhA = buildMeshBVH(ctx.meshDataA, ctx.eps);
    BVH bvhB = buildMeshBVH(ctx.meshDataB, ctx.eps);
    std::vector<bool> removeAOutsideB = getRemovalMask(ctx.meshDataA, ctx.meshDataB, bvhB, ctx.eps, BoolOp::DifferenceBase);
    std::vector<bool> removeBInsideA  = getRemovalMask(ctx.meshDataB, ctx.meshDataA, bvhA, ctx.eps, BoolOp::DifferenceTool);
    std::vector<bool> removeBOutsideA = getRemovalMask(ctx.meshDataB, ctx.meshDataA, bvhA, ctx.eps, BoolOp::DifferenceBase);
    std::vector<bool> removeAInsideB  = getRemovalMask(ctx.meshDataA, ctx.meshDataB, bvhB, ctx.eps, BoolOp::DifferenceTool);
    MeshData aShell = ctx.meshDataA;
    MeshData bCavity = ctx.meshDataB;
    filterMesh(aShell, removeAOutsideB);
    filterMesh(bCavity, removeBInsideA);
    invertWinding(bCavity);
    meshA = combineMeshes({aShell, bCavity}, ctx.eps);
    MeshData bShell = ctx.meshDataB;
    MeshData aCavity = ctx.meshDataA;
    filterMesh(bShell, removeBOutsideA);
    filterMesh(aCavity, removeAInsideB);
    invertWinding(aCavity);
    meshB = combineMeshes({bShell, aCavity}, ctx.eps);
    return ctx.eps;
}
void meshCombine(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool removeTouchingSurfaces) {
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    double eps = meshCombine(meshA, meshB, removeTouchingSurfaces);
    rebuildMesh(A, meshA, eps);
    rebuildMesh(B, meshB, eps);
}
bem::TriangleMesh<3> meshUnion(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B) {
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    double eps = 0.0;
    MeshData combined = meshUnion(meshA, meshB, &eps);
    bem::TriangleMesh<3> out;
    rebuildMesh(out, combined, eps);
    return out;
}
bem::TriangleMesh<3> meshIntersect(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B) {
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    double eps = 0.0;
    MeshData combined = meshIntersect(meshA, meshB, &eps);
    bem::TriangleMesh<3> out;
    rebuildMesh(out, combined, eps);
    return out;
}
void meshDifference(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B) {
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    double eps = meshDifference(meshA, meshB);
    rebuildMesh(A, meshA, eps);
    rebuildMesh(B, meshB, eps);
}