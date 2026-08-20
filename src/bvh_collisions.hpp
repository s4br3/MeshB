#pragma once
#include "vector.hpp"
#include "bvh.hpp"
#include "mesh_types.hpp"
#include <optional>

/**
* @brief Computes non-coplanar intersection line segment between two triangles, if one exists.
* @param[in] vertices1 - Array of Vec3 coordinates of first triangle's vertices.
* @param[in] n1 - Normal vector of t1.
* @param[in] c1 - Centroid of t1.
* @param[in] vertices2 - Array of Vec3 coordinates of second triangle's vertices.
* @param[in] n2 - Normal vector of t2.
* @param[in] c2 - Centroid of t2.
* @param[in] eps - Tolerance value.
* @return Optional pair of points defining intersection segment; std::nullopt if no collision.
*/
std::optional<std::pair<Vec3, Vec3>> findIntersectionPointsNC(
    const TriVerts& vertices1, const Vec3& n1, const Vec3& c1,
    const TriVerts& vertices2, const Vec3& n2, const Vec3& c2,
    const double eps);

/**
* @brief Computes intersection points between two coplanar triangles.
* @param[in] vertices1 - Array of Vec3 coordinates of first triangle's vertices.
* @param[in] vertices2 - Array of Vec3 coordinates of second triangle's vertices.
* @param[in] n2 - Shared surface normal.
* @param[in] c2 - Reference center point on plane.
* @param[in] eps - Tolerance value.
* @return Polyline defining overlapping region boundaries.
*/
const std::vector<Vec3>& findIntersectionPointsC(
    const TriVerts& vertices1,
    const TriVerts& vertices2, const Vec3& n2, const Vec3& c2,
    const double eps);

/**
* @brief Traverses two BVH trees to find all surface segment and face intersections between two meshes.
* @param[in] bvh1 - Bounding volume hierarchy of mesh 1.
* @param[in] bvh2 - Bounding volume hierarchy of mesh 2.
* @param[in] ctx - Full Collision Context struct containing both meshes and several output buffers for information.
*/
void findAllCollisions(const BVH& bvh1, const BVH& bvh2, CollisionContext& ctx);

/**
* @brief Constructs a Bounding Volume Hierarchy (BVH) tree for accelerated collision queries.
* @param[in] mesh - Target surface mesh.
* @param[in] eps - Tolerance value.
* @return Built BVH tree.
 */
BVH buildMeshBVH(const MeshData& mesh, double eps);

/**
* @brief High-level detection pipeline computing all geometric collisions between two meshes.
* @param[in] meshA - First input surface mesh.
* @param[in] meshB - Second input surface mesh.
* @return CollisionContext containing collision structures, cutting lines, and element indices.
*/
CollisionContext detectCollisions(const MeshData& meshA, const MeshData& meshB);