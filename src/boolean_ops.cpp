#include "boolean_ops.hpp"
#include "bem_interop.hpp"
#include "bvh_collisions.hpp"
#include "geom_3d.hpp"
#include "mesh_types.hpp"
#include "triangulation.hpp"
#include "classify.hpp"
#include "mesh_clean.hpp"
#include <unordered_set>
#include <cstddef>
#include <queue>
#include <utility>
void verifyMeshConformity(CollisionContext& ctx)
{
    auto checkMeshCuts = [&](const MeshData& mesh,
                             const std::unordered_map<size_t, PolyLine>& NCcoords,
                             const std::unordered_map<size_t, std::vector<PolyLine>>& Ccoords)
        -> bool
    {
        SpatialGrid3D grid(ctx.eps);
        for (const Vec3& node : mesh.nodes)
            grid.getOrAdd(node);
        const size_t meshPointCount = grid.getUniquePoints().size();
        auto isNodeInMesh = [&](const Vec3& point) -> bool
        {
            const size_t index = grid.getOrAdd(point);
            return index < meshPointCount;
        };
        auto checkPolyline = [&](const PolyLine& polyline) -> bool
        {
            for (const auto& [start, end] : polyline)
            {
                if (!isNodeInMesh(start) || !isNodeInMesh(end))
                    return false;
            }
            return true;
        };
        for (const auto& [_, polyline] : NCcoords)
        {
            if (!checkPolyline(polyline))
                return false;
        }
        for (const auto& [_, polylineList] : Ccoords)
        {
            for (const auto& polyline : polylineList)
            {
                if (!checkPolyline(polyline))
                    return false;
            }
        }
        return true;
    };
    ctx.conformal = checkMeshCuts(ctx.meshDataA, ctx.NCAcoords, ctx.CAcoords) &&
           checkMeshCuts(ctx.meshDataB, ctx.NCBcoords, ctx.CBcoords);
}
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
            FaceClass currentClass = classifyFace(targetBVH, targetMesh, mesh.centres[i], mesh.normals[i], eps);
            if (currentClass == FaceClass::Outside){
                double wn = evaluateGWN(mesh.centres[i], targetMesh, targetBVH, gwnData);
                currentClass = (std::abs(wn) > 0.5) ? FaceClass::Inside : FaceClass::Outside;
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
            FaceClass currentClass = classifyFace(targetBVH, targetMesh, mesh.centres[i], mesh.normals[i], eps);
            if (currentClass == FaceClass::Outside){
                double wn = evaluateGWN(mesh.centres[i], targetMesh, targetBVH, gwnData);
                currentClass = (std::abs(wn) > 0.5) ? FaceClass::Inside : FaceClass::Outside;
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
Connection nonConformal(const MeshData& meshA, const MeshData& meshB){
    BVH bvhA = buildMeshBVH(meshA, 0.0);
    BVH bvhB = buildMeshBVH(meshB, 0.0);
    double eps = 1e-7;
    if (!bvhA.nodes.empty() && !bvhB.nodes.empty()) {
        eps = computeMeshesEpsilon(bvhA.nodes[0].getBbox(), bvhB.nodes[0].getBbox());
    }
    bvhA = buildMeshBVH(meshA, eps);
    bvhB = buildMeshBVH(meshB, eps);
    Connection conn;
    const double areaEps = eps * eps;
    std::vector<double> triAreasA(meshA.triangles.size(), 0.0);
    for(size_t i = 0; i < meshA.triangles.size(); ++i) {
        Vec3 v0 = meshA.nodes[meshA.triangles[i].v[0]];
        Vec3 v1 = meshA.nodes[meshA.triangles[i].v[1]];
        Vec3 v2 = meshA.nodes[meshA.triangles[i].v[2]];
        Vec3 cr = cross(v1 - v0, v2 - v0);
        triAreasA[i] = 0.5 * std::sqrt(dot(cr, cr));
    }
    std::vector<double> triAreasB(meshB.triangles.size(), 0.0);
    for(size_t i = 0; i < meshB.triangles.size(); ++i) {
        Vec3 v0 = meshB.nodes[meshB.triangles[i].v[0]];
        Vec3 v1 = meshB.nodes[meshB.triangles[i].v[1]];
        Vec3 v2 = meshB.nodes[meshB.triangles[i].v[2]];
        Vec3 cr = cross(v1 - v0, v2 - v0);
        triAreasB[i] = 0.5 * std::sqrt(dot(cr, cr));
    }
    std::vector<std::pair<size_t, size_t>> stack;
    if (!bvhA.nodes.empty() && !bvhB.nodes.empty()) stack.emplace_back(0, 0);
    std::unordered_map<size_t, std::unordered_map<size_t, std::vector<Vec3>>> rawPts;
    std::unordered_map<size_t, std::unordered_map<size_t, double>> rawAreas;
    while (!stack.empty()) {
        auto [node1Idx, node2Idx] = stack.back();
        stack.pop_back();
        const Node& node1 = bvhA.nodes[node1Idx];
        const Node& node2 = bvhB.nodes[node2Idx];
        if (!boundingBoxOverlap(node1.getBbox(), node2.getBbox(), eps)) continue;
        if (node1.isLeaf() && node2.isLeaf()) {
            for (size_t i = node1.firstId; i < node1.firstId + node1.primCount; ++i) {
                size_t prim1Id = bvhA.primIds[i];
                for (size_t j = node2.firstId; j < node2.firstId + node2.primCount; ++j) {
                    size_t prim2Id = bvhB.primIds[j];
                    const TriVerts& t1 = meshA.triangles[prim1Id].getVertices(meshA.nodes);
                    const TriVerts& t2 = meshB.triangles[prim2Id].getVertices(meshB.nodes);
                    const Vec3& n1 = meshA.normals[prim1Id];
                    const Vec3& n2 = meshB.normals[prim2Id];
                    const Vec3& c1 = meshA.centres[prim1Id];
                    const Vec3& c2 = meshB.centres[prim2Id];
                    if (coplanar(n1, c1, n2, c2, eps)) {
                        std::vector<Vec3> res = findIntersectionPointsC(t1, t2, n2, c2, eps);
                        if (res.size() > 1) {
                            rawPts[prim1Id][prim2Id].insert(rawPts[prim1Id][prim2Id].end(), res.begin(), res.end());
                            PolyLine temp;
                            for (size_t k = 0; k < res.size(); k++) {
                                temp.push_back({res[k], res[(k + 1) % res.size()]});
                            }
                            rawAreas[prim1Id][prim2Id] += computePolyLineArea(temp);
                        }
                    } else {
                        if (auto optResult = findIntersectionPointsNC(t1, n1, c1, t2, n2, c2, eps)) {
                            auto [startPt, endPt] = *optResult;
                            rawPts[prim1Id][prim2Id].push_back(startPt);
                            rawPts[prim1Id][prim2Id].push_back(endPt);
                        }
                    }
                }
            }
            continue;
        }
        if (node1.isLeaf()) {
            stack.emplace_back(node1Idx, node2.firstId);
            stack.emplace_back(node1Idx, node2.firstId + 1);
        } else if (node2.isLeaf()) {
            stack.emplace_back(node1.firstId, node2Idx);
            stack.emplace_back(node1.firstId + 1, node2Idx);
        } else {
            stack.emplace_back(node1.firstId, node2.firstId);
            stack.emplace_back(node1.firstId, node2.firstId + 1);
            stack.emplace_back(node1.firstId + 1, node2.firstId);
            stack.emplace_back(node1.firstId + 1, node2.firstId + 1);
        }
    }
    for (auto& [idA, mapB] : rawPts) {
        for (auto& [idB, pts] : mapB) {
            conn.overlapsAB[idA].push_back(idB);
            conn.overlapsBA[idB].push_back(idA);
            SpatialGrid3D grid(eps);
            for (const auto& pt : pts) grid.getOrAdd(pt);
            conn.intersectionsAB[idA][idB] = grid.getUniquePoints();
            conn.intersectionsBA[idB][idA] = grid.getUniquePoints();
            double overlapArea = rawAreas[idA][idB];
            double percentAB = 0.0;
            double percentBA = 0.0;
            if (triAreasA[idA] > areaEps && overlapArea > areaEps) {
                percentAB = (overlapArea / triAreasA[idA]) * 100.0;
                percentAB = std::max(0.0, std::min(100.0, percentAB));
                percentBA = (overlapArea / triAreasB[idB]) * 100.0;
                percentBA = std::max(0.0, std::min(100.0, percentBA));
            }
            conn.overlapPercentAB[idA][idB] = percentAB;
            conn.overlapPercentBA[idB][idA] = percentBA;
        }
    }
    return conn;
}
Connection nonConformal(const bem::TriangleMesh<3>& A, const bem::TriangleMesh<3>& B) {
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    return nonConformal(meshA, meshB);
}
CollisionContext collideAndCut(MeshData& meshA, MeshData& meshB, bool removeTouchingSurfaces) {
    CollisionContext ctx = detectCollisions(meshA, meshB);
    verifyMeshConformity(ctx);
    if (ctx.conformal){
        std::cout << "Meshes are already conformal: skipping cutting stage\n";
        return ctx;
    }
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