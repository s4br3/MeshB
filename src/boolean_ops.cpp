#include "boolean_ops.hpp"
#include "bem_interop.hpp"
#include "math_utils.hpp"
#include "bvh_collisions.hpp"
#include "triangulation.hpp"
#include "raycast.hpp"
#include "mesh_clean.hpp"
#include <cstddef>
std::vector<bool> getRemovalMask(const MeshData& meshA, const MeshData& meshB, double eps, BoolOp op, bool isMeshA) {
    BVH bvhB = buildMeshBVH(meshB, eps);
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
                if (fclass == FaceClass::Outside) keep = isMeshA;
                else if (fclass == FaceClass::Inside) keep = !isMeshA;
                else if (fclass == FaceClass::CoplanarSame) keep = false;
                else if (fclass == FaceClass::CoplanarSame) keep = isMeshA;
                break;
        }
        removeA[i] = !keep;
    }
    return removeA;
}
Connection nonConformal(const bem::TriangleMesh<3>& A, const bem::TriangleMesh<3>& B)
{
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    CollisionContext ctx = detectCollisions(meshA, meshB);

    auto processSide = [&](
        const std::unordered_map<size_t, PolyLine>& NCcoords,
        const std::unordered_map<size_t, std::vector<PolyLine>>& Ccoords,
        const std::unordered_map<size_t, std::vector<size_t>>& tris,
        const MeshData& srcMesh,
        const MeshData& targetMesh,
        std::unordered_map<size_t, std::vector<Vec3>>& outCoords,
        std::unordered_map<size_t, std::vector<size_t>>& outTaggedTris)
    {
        std::unordered_map<size_t, SpatialGrid3D> grids;
        auto addSegmentsToGrid = [&](size_t idx, const PolyLine& poly) {
            size_t tag = srcMesh.tags[idx];
            auto [it, _] = grids.try_emplace(tag, ctx.eps);
            for (const std::pair<Vec3, Vec3>& segment : poly) {
                it->second.getOrAdd(segment.first);
                it->second.getOrAdd(segment.second);
            }
        };
        for (const auto& [idx, poly] : NCcoords) {
            addSegmentsToGrid(idx, poly);
        }
        for (const auto& [idx, poly_list] : Ccoords) {
            for (const auto& poly : poly_list) {
                addSegmentsToGrid(idx, poly);
            }
        }
        for (const auto& [tag, grid] : grids) {
            outCoords[tag] = grid.getUniquePoints();
        }
        for (const auto& [idx, list] : tris) {
            size_t tagSrc = srcMesh.tags[idx];
            std::vector<size_t>& new_list = outTaggedTris[tagSrc];
            for (size_t target_idx : list) {
                new_list.push_back(targetMesh.tags[target_idx]);
            }
        }
    };
    std::unordered_map<size_t, std::vector<Vec3>> aCoords, bCoords;
    std::unordered_map<size_t, std::vector<size_t>> tagged_Atris, tagged_Btris;

    processSide(ctx.NCAcoords, ctx.CAcoords, ctx.Atris, ctx.meshDataA, ctx.meshDataB, aCoords, tagged_Atris);
    processSide(ctx.NCBcoords, ctx.CBcoords, ctx.Btris, ctx.meshDataB, ctx.meshDataA, bCoords, tagged_Btris);
    return Connection(aCoords, bCoords, tagged_Atris, tagged_Btris);
}
CollisionContext collideAndCut(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool removeTouchingSurfaces) {
    MeshData meshA = extractMeshData(A);
    MeshData meshB = extractMeshData(B);
    CollisionContext ctx = detectCollisions(meshA, meshB);
    MeshData newMeshA = cutMesh(ctx.meshDataA, ctx.NCAcoords, ctx.CAcoords, ctx.eps, removeTouchingSurfaces, true);
    MeshData newMeshB = cutMesh(ctx.meshDataB, ctx.NCBcoords, ctx.CBcoords, ctx.eps, removeTouchingSurfaces, false);
    ctx.meshDataA = newMeshA;
    ctx.meshDataB = newMeshB;
    return ctx;
}
void meshCombine(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool removeTouchingSurfaces, bool cleanDegenerate) {
    CollisionContext ctx = collideAndCut(A, B, removeTouchingSurfaces);
    if (cleanDegenerate) {
        cleanMesh(ctx.meshDataA, ctx.eps);
        cleanMesh(ctx.meshDataB, ctx.eps);
    }
    rebuildMesh(A, ctx.meshDataA, ctx.eps);
    rebuildMesh(B, ctx.meshDataB, ctx.eps);
}
bem::TriangleMesh<3> meshUnion(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool cleanDegenerate) {
    CollisionContext ctx = collideAndCut(A, B, false);
    std::vector<bool> removeAInB = getRemovalMask(ctx.meshDataA, ctx.meshDataB, ctx.eps, BoolOp::Union, true);
    std::vector<bool> removeBInA = getRemovalMask(ctx.meshDataB, ctx.meshDataA, ctx.eps, BoolOp::Union, false);
    filterMesh(ctx.meshDataA, removeAInB);
    filterMesh(ctx.meshDataB, removeBInA);
    MeshData combined = combineMeshes({ctx.meshDataA, ctx.meshDataB}, ctx.eps);
    if (cleanDegenerate) {
        cleanMesh(combined, ctx.eps);
    }
    bem::TriangleMesh<3> out;
    rebuildMesh(out, combined, ctx.eps);
    return out;
}
bem::TriangleMesh<3> meshIntersect(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool cleanDegenerate) {
    CollisionContext ctx = collideAndCut(A, B, false);
    std::vector<bool> removeANotInB = getRemovalMask(ctx.meshDataA, ctx.meshDataB, ctx.eps, BoolOp::Intersect, true);
    std::vector<bool> removeBNotInA = getRemovalMask(ctx.meshDataB, ctx.meshDataA, ctx.eps, BoolOp::Intersect, false);
    filterMesh(ctx.meshDataA, removeANotInB);
    filterMesh(ctx.meshDataB, removeBNotInA);
    MeshData combined = combineMeshes({ctx.meshDataA, ctx.meshDataB}, ctx.eps);
    if (cleanDegenerate) {
        cleanMesh(combined, ctx.eps);
    }
    bem::TriangleMesh<3> out;
    rebuildMesh(out, combined, ctx.eps);
    return out;
}
void meshDifference(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool cleanDegenerate) {
    CollisionContext ctx = collideAndCut(A, B, false);
    std::vector<bool> removeAInB = getRemovalMask(ctx.meshDataA, ctx.meshDataB, ctx.eps, BoolOp::Difference, true);
    std::vector<bool> removeBInA = getRemovalMask(ctx.meshDataB, ctx.meshDataA, ctx.eps, BoolOp::Difference, true);
    std::vector<bool> removeBNotInA = removeBInA;
    std::vector<bool> removeANotInB = removeAInB;
    for (auto it = removeBNotInA.begin(), e = removeBNotInA.end(); it != e; ++it) *it = !*it;
    for (auto it = removeANotInB.begin(), e = removeANotInB.end(); it != e; ++it) *it = !*it;
    MeshData bInA = ctx.meshDataB;
    MeshData aInB = ctx.meshDataA;
    filterMesh(ctx.meshDataA, removeAInB);
    filterMesh(ctx.meshDataB, removeBInA);
    // Filter interior cut boundaries
    filterMesh(bInA, removeBNotInA);
    filterMesh(aInB, removeANotInB);
    invertWinding(bInA);
    invertWinding(aInB);
    MeshData combinedA = combineMeshes({ctx.meshDataA, bInA}, ctx.eps);
    MeshData combinedB = combineMeshes({ctx.meshDataB, aInB}, ctx.eps);
    if (cleanDegenerate) {
        cleanMesh(combinedA, ctx.eps);
        cleanMesh(combinedB, ctx.eps);
    }
    rebuildMesh(A, combinedA, ctx.eps);
    rebuildMesh(B, combinedB, ctx.eps);
}