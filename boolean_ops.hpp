#pragma once
#include "mesh_types.hpp"
#include <source/geometry/mesh/triangle_mesh.hpp>
enum class BoolOp { Union, Intersect, Difference };
MeshData classifyAndFilterAB(
    const MeshData& meshA, const MeshData& meshB,
    double eps, BoolOp op, bool isMeshA);
Connection nonConformal(const bem::TriangleMesh<3>& meshA, const bem::TriangleMesh<3>& meshB);
CollisionContext collideAndCut(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool removeTouchingSurfaces = false);
void meshCombine(bem::TriangleMesh<3>& meshA, bem::TriangleMesh<3>& meshB, bool removeTouchingSurfaces = false, bool cleanDegenerate = false);
bem::TriangleMesh<3> meshUnion(bem::TriangleMesh<3>& meshA, bem::TriangleMesh<3>& meshB, bool cleanDegenerate = false);
bem::TriangleMesh<3> meshIntersect(bem::TriangleMesh<3>& meshA, bem::TriangleMesh<3>& meshB, bool cleanDegenerate = false);
void meshDifference(bem::TriangleMesh<3>& meshA, bem::TriangleMesh<3>& meshB, bool cleanDegenerate = false);
