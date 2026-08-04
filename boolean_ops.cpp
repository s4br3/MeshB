#include "boolean_ops.hpp"
#include "bem_interop.hpp"
#include "bvh_collisions.hpp"
#include "mesh_types.hpp"
#include "raycast.hpp"
#include "triangulation.hpp"
#include "mesh_clean.hpp"
#include <cstddef>
MeshData classifyAndFilterAB(const MeshData& meshA, const MeshData& meshB, double eps, BoolOp op, bool isMeshA){
    Bvh bvhB = buildMeshBVH(meshB);
    const size_t aN = meshA.triangles.size();
    std::vector<bool> removeA(aN, false);
    for (size_t i = 0; i < aN; ++i) {
        const Vec3& center = meshA.centres[i];
        const Vec3& normal = meshA.normals[i];
        FaceClass fclass = classifyFace(bvhB, meshB, center, normal, eps);
        bool keep = false;
        switch (op) {
            case BoolOp::Union:
                if (fclass == FaceClass::Outside)       keep = true;
                else if (fclass == FaceClass::Inside)  keep = false;
                else if (fclass == FaceClass::CoplanarSame) keep = isMeshA;
                else if (fclass == FaceClass::CoplanarOpp)  keep = false;
                break;
            case BoolOp::Intersect:
                if (fclass == FaceClass::Outside)       keep = false;
                else if (fclass == FaceClass::Inside)  keep = true;
                else if (fclass == FaceClass::CoplanarSame) keep = isMeshA;
                else if (fclass == FaceClass::CoplanarOpp)  keep = false;
                break;
            case BoolOp::Difference:
                if (isMeshA) {
                    if (fclass == FaceClass::Outside)       keep = true;
                    else if (fclass == FaceClass::Inside)  keep = false;
                    else if (fclass == FaceClass::CoplanarSame) keep = false;
                    else if (fclass == FaceClass::CoplanarOpp)  keep = true;
                } else {
                    if (fclass == FaceClass::Outside)       keep = false;
                    else if (fclass == FaceClass::Inside)  keep = true;
                    else if (fclass == FaceClass::CoplanarSame) keep = false;
                    else if (fclass == FaceClass::CoplanarOpp)  keep = false;
                }
                break;
        }
        removeA[i] = !keep;
    }
    MeshData out = filterMesh(meshA, removeA);
    clearUnusedNodes(out);
    return out;
}
Connection nonConformal(const bem::TriangleMesh<3>& A, const bem::TriangleMesh<3>& B)
{
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    CollisionContext ctx = detectCollisions(meshA, meshB);
    std::unordered_map<size_t, SpatialGrid3D> gridA;
    std::unordered_map<size_t, SpatialGrid3D> gridB;
    std::unordered_map<size_t, std::vector<size_t>> tagged_Atris;
    std::unordered_map<size_t, std::vector<size_t>> tagged_Btris;
    auto addSegmentsToGrid = [&](size_t idx, const PolyLine& poly, 
                                 std::unordered_map<size_t, SpatialGrid3D>& grid, 
                                 size_t tag) {
        auto [it, _] = grid.try_emplace(tag, ctx.eps);
        for (const auto& segment : poly) {
            it->second.getOrAdd(segment.first);
            it->second.getOrAdd(segment.second);
        }
    };
    for (const auto& [idx, poly] : ctx.NCAcoords) {
        addSegmentsToGrid(idx, poly, gridA, ctx.meshDataA.tags[idx]);
    }
    for (const auto& [idx, poly_list] : ctx.CAcoords) {
        for (const auto& poly : poly_list) {
            addSegmentsToGrid(idx, poly, gridA, ctx.meshDataA.tags[idx]);
        }
    }
    for (const auto& [idx, poly] : ctx.NCBcoords) {
        addSegmentsToGrid(idx, poly, gridB, ctx.meshDataB.tags[idx]);
    }
    for (const auto& [idx, poly_list] : ctx.CBcoords) {
        for (const auto& poly : poly_list) {
            addSegmentsToGrid(idx, poly, gridB, ctx.meshDataB.tags[idx]);
        }
    }
    std::unordered_map<size_t, std::vector<Vec3>> aCoords;
    std::unordered_map<size_t, std::vector<Vec3>> bCoords;
    for (const auto& [tag, grid] : gridA) {
        aCoords[tag] = grid.getUniquePoints();
    }
    for (const auto& [tag, grid] : gridB) {
        bCoords[tag] = grid.getUniquePoints();
    }
    for (const auto& [idx, list] : ctx.Atris) {
        size_t tagA = ctx.meshDataA.tags[idx];
        auto& new_list = tagged_Atris[tagA];
        for (size_t target_idx : list) {
            new_list.push_back(ctx.meshDataB.tags[target_idx]);
        }
    }
    for (const auto& [idx, list] : ctx.Btris) {
        size_t tagB = ctx.meshDataB.tags[idx];
        auto& new_list = tagged_Btris[tagB];
        for (size_t target_idx : list) {
            new_list.push_back(ctx.meshDataA.tags[target_idx]);
        }
    }
    return Connection(aCoords, bCoords, tagged_Atris, tagged_Btris);
}
CollisionContext collideAndCut(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool removeTouchingSurfaces){
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    CollisionContext ctx = detectCollisions(meshA, meshB);
    MeshData newMeshA = cutMesh(ctx.meshDataA, ctx.NCAcoords, ctx.CAcoords, ctx.eps, removeTouchingSurfaces);
    MeshData newMeshB = cutMesh(ctx.meshDataB, ctx.NCBcoords, ctx.CBcoords, ctx.eps, removeTouchingSurfaces);
    ctx.meshDataA = newMeshA;
    ctx.meshDataB = newMeshB;
    return ctx;
}
void meshCombine(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool removeTouchingSurfaces, bool cleanDegenerate) {
    CollisionContext ctx = collideAndCut(A, B, removeTouchingSurfaces);
    if (cleanDegenerate){
        cleanMesh(ctx.meshDataA, ctx.eps);
        cleanMesh(ctx.meshDataB, ctx.eps);
    }
    clearUnusedNodes(ctx.meshDataA);
    clearUnusedNodes(ctx.meshDataB);
    rebuildMesh(A, ctx.meshDataA, ctx.eps);
    rebuildMesh(B, ctx.meshDataB, ctx.eps);
}
bem::TriangleMesh<3> meshUnion(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool cleanDegenerate){
    CollisionContext ctx = collideAndCut(A, B, false);
    MeshData newMeshA = classifyAndFilterAB(ctx.meshDataA, ctx.meshDataB, ctx.eps, BoolOp::Union, true);
    MeshData newMeshB = classifyAndFilterAB(ctx.meshDataB, ctx.meshDataA, ctx.eps, BoolOp::Union, false);
    MeshData combined = combineMeshes({newMeshA, newMeshB}, ctx.eps);
    if (cleanDegenerate){
        cleanMesh(combined, ctx.eps);
    }
    clearUnusedNodes(combined);
    bem::TriangleMesh<3> out;
    rebuildMesh(out, combined, ctx.eps);
    return out;
}
bem::TriangleMesh<3> meshIntersect(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool cleanDegenerate){
    CollisionContext ctx = collideAndCut(A, B, false);
    MeshData newMeshA = classifyAndFilterAB(ctx.meshDataA, ctx.meshDataB, ctx.eps, BoolOp::Intersect, true);
    MeshData newMeshB = classifyAndFilterAB(ctx.meshDataB, ctx.meshDataA, ctx.eps, BoolOp::Intersect, false);
    MeshData combined = combineMeshes({newMeshA, newMeshB}, ctx.eps);
    if (cleanDegenerate){
        cleanMesh(combined, ctx.eps);
    }
    clearUnusedNodes(combined);
    bem::TriangleMesh<3> out;
    rebuildMesh(out, combined, ctx.eps);
    return out;
}
void meshDifference(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool cleanDegenerate){
    CollisionContext ctx = collideAndCut(A, B, false);
    MeshData newMeshA = classifyAndFilterAB(ctx.meshDataA, ctx.meshDataB, ctx.eps, BoolOp::Difference, true);
    MeshData bInA = classifyAndFilterAB(ctx.meshDataB, ctx.meshDataA, ctx.eps, BoolOp::Difference, false);
    invertWinding(bInA);
    MeshData combinedA = combineMeshes({newMeshA, bInA}, ctx.eps);
    MeshData newMeshB = classifyAndFilterAB(ctx.meshDataB, ctx.meshDataA, ctx.eps, BoolOp::Difference, true);
    MeshData aInB = classifyAndFilterAB(ctx.meshDataA, ctx.meshDataB, ctx.eps, BoolOp::Difference, false);
    invertWinding(aInB);
    MeshData combinedB = combineMeshes({newMeshB, aInB}, ctx.eps);
    if (cleanDegenerate){
        cleanMesh(combinedA, ctx.eps);
        cleanMesh(combinedB, ctx.eps);
    }
    clearUnusedNodes(combinedA);
    clearUnusedNodes(combinedB);
    rebuildMesh(A, combinedA, ctx.eps);
    rebuildMesh(B, combinedB, ctx.eps);
}