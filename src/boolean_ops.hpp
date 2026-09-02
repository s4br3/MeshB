#pragma once
#include "mesh_types.hpp"
#include <source/geometry/mesh/triangle_mesh.hpp>

/**
* @enum BoolOp
* @brief CSG (Constructive Solid Geometry) boolean operations supported on closed surface meshes.
*/
enum class BoolOp {Union, Intersect, DifferenceBase, DifferenceTool};

/**
* @brief Overload for evaluating surface triangle classification and removal mask via a BVH tree.
* @param[in] mesh - Primary surface mesh to filter.
* @param[in] targetMesh - Target surface mesh evaluated against.
* @param[in] targetBVH - Prebuilt Bounding Volume Hierarchy tree of target mesh.
* @param[in] eps - Spatial distance tolerance.
* @param[in] op - Selected CSG operation rule.
* @return Boolean removal mask array.
*/
std::vector<bool> getRemovalMask(
    const MeshData& mesh, const MeshData& targetMesh, 
    const BVH& targetBVH, double eps, BoolOp op);

/**
* @brief Generates non-conformal mesh interface connections between two independent meshes.
* @param[in] meshA - Primary surface mesh.
* @param[in] meshB - Secondary surface mesh.
* @return Connection object containing mapping data across non-matching boundaries and intersections.
*/
Connection nonConformal(const MeshData& meshA, const MeshData& meshB);

/**
* @brief Generates non-conformal mesh interface connections between two independent meshes.
* @param[in] meshA - First openBEM surface mesh.
* @param[in] meshB - Second openBEM surface mesh.
* @return Connection object containing mapping data across non-matching boundaries and intersections.
*/
Connection nonConformal(const bem::TriangleMesh<3>& meshA, const bem::TriangleMesh<3>& meshB);

/**
* @brief Set ctx.conformal to whether the meshes stored in the collision context
* are conformal to one another.
* @param[in] ctx - Collision context to evaluate.
*/
void verifyMeshConformity(CollisionContext& ctx);

/**
* @brief Low-level pipeline computing collisions and cutting operations between two openBEM meshes.
* @param[in,out] A - First openBEM surface mesh.
* @param[in,out] B - Second openBEM surface mesh.
* @param[in] removeTouchingSurfaces Flag to strip coplanar co-facing surfaces.
* @return CollisionContext holding mapped intersection features.
*/
CollisionContext collideAndCut(bem::TriangleMesh<3>& A, bem::TriangleMesh<3>& B, bool removeTouchingSurfaces = false);
/**
* @brief Retriangulates and cuts two MeshData structures to guarantee boundary conformity.
* @param[in,out] meshA - First surface mesh modified in-place.
* @param[in,out] meshB - Second surface mesh modified in-place.
* @param[in] removeTouchingSurfaces - Option to strip coplanar co-facing surfaces.
* @return Computed spatial epsilon tolerance.
*/
double meshCombine(MeshData& meshA, MeshData& meshB, bool removeTouchingSurfaces = false);

/**
* @brief Computes CSG Union of two surface MeshData objects.
* @param[in,out] meshA - Primary surface mesh.
* @param[in,out] meshB - Secondary surface mesh.
* @param[out] outEps - Optional output pointer to store numerical tolerance.
* @return Resulting merged surface MeshData.
*/
MeshData meshUnion(MeshData& meshA, MeshData& meshB, double* outEps = nullptr);

/**
* @brief Computes CSG Intersection ($A \cap B$) of two surface MeshData objects.
* @param[in,out] meshA - Primary surface mesh.
* @param[in,out] meshB - Secondary surface mesh.
* @param[out] outEps - Optional output pointer to store numerical tolerance.
* @return Intersected MeshData object.
*/
MeshData meshIntersect(MeshData& meshA, MeshData& meshB, double* outEps = nullptr);

/**
* @brief Computes CSG Difference ($A \setminus B$) of two surface MeshData objects in-place.
* @param[in,out] meshA - Primary surface mesh modified in-place to $A \setminus B$.
* @param[in,out] meshB - Secondary surface mesh modified in-place to $B \setminus A$.
* @return Computed spatial epsilon tolerance.
*/
double meshDifference(MeshData& meshA, MeshData& meshB);
/**
* @brief Combines two meshes into a unified mesh surface without evaluating CSG interior removal.
* @param[in,out] meshA - First openBEM surface mesh at input and the modified 
(retriangulated and cut) version of the same mesh.
* @param[in,out] meshB - Second openBEM surface mesh at input and the modified
(retriangulated and cut) version of the same mesh
* @param[in] removeTouchingSurfaces - Flag to remove identical overlapping elements.
* @param[in] cleanDegenerate - Flag to remove degenerate or collapsed elements after combining.
* @details Despite `combining` the meshes, it maintains mesh A and mesh B as separate entities,
* just triangulated to ensure conformity
*/
void meshCombine(bem::TriangleMesh<3>& meshA, bem::TriangleMesh<3>& meshB, bool removeTouchingSurfaces = false);

/**
* @brief Computes CSG Union of two surface meshes.
* @param[in,out] meshA - Primary surface mesh object.
* @param[in,out] meshB - Secondary surface mesh object.
* @param[in] cleanDegenerate - Option to strip degenerate elements from final surface.
* @return Resulting merged openBEM with all the triangles on the "outside" of both initial meshes
*/
bem::TriangleMesh<3> meshUnion(bem::TriangleMesh<3>& meshA, bem::TriangleMesh<3>& meshB);

/**
* @brief Computes CSG Intersection ($A \cap B$) of two surface meshes.
* @param[in,out] meshA - Primary surface mesh.
* @param[in,out] meshB - Secondary surface mesh.
* @param[in] cleanDegenerate - Option to strip degenerate elements from final surface.
* @return Intersected openBEM mesh with all the triangles from mesh A that are inside B, and vice versa
*/
bem::TriangleMesh<3> meshIntersect(bem::TriangleMesh<3>& meshA, bem::TriangleMesh<3>& meshB);

/**
* @brief Computes symmetric CSG Difference of two surface meshes
* @param[in,out] meshA - Target surface mesh to subtract from.
* @param[in,out] meshB - Tool surface mesh to subtract.
* @param[in] cleanDegenerate - Option to strip degenerate elements from final surface.
*/
void meshDifference(bem::TriangleMesh<3>& meshA, bem::TriangleMesh<3>& meshB);