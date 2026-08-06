#pragma once
#include "mesh_types.hpp"
#include <bvh/v2/bvh.h>
#include <bvh/v2/node.h>
#include <optional>
#include <bvh/v2/default_builder.h>

using Node = BVH::Node<double, 3>;
using TriVerts = std::array<Vec3, 3>;
using Bvh  = BVH::Bvh<Node>;

/**
* @brief Determines if a 3D point lies inside or on the boundary of a triangle element.
* @param[in] p - Point position to test.
* @param[in] vertices - Array of Vec3 coordinates of vertices
* @param[in] eps - Geometric tolerance.
* @return True if point lies within triangle bounds.
* @details Barycentric PIT test, so it's completely separated from any issues with normal pointing direction
*/
bool pointInTriangle(const Vec3& p, const TriVerts& vertices, double eps);

/**
* @brief Calculates signed distances from the 3 vertices of triangle 1 to the supporting plane of triangle 2.
* @param[in] vertices - Array of Vec3 coordinates of vertices
* @param[in] n2 - Normal vector of second plane/triangle.
* @param[in] centre2 - Point lying on the second plane.
* @return Array of signed distances for t1's three vertices.
*/
std::array<double, 3> distancesToPlane(const TriVerts& vertices, const Vec3& n2, const Vec3& centre2);

/**
* @brief Tests whether two 3D axis-aligned bounding boxes overlap within tolerance.
* @param[in] box1 - First bounding box.
* @param[in] box2 - Second bounding box.
* @param[in] eps - Geometric tolerance.
* @return True if boxes overlap; false otherwise.
*/
bool boundingBoxOverlap(const BBox& box1, const BBox& box2, double eps);

/**
* @brief Checks if two triangles are coplanar within specified normal and distance thresholds.
* @param[in] n1 - Normal vector of first triangle.
* @param[in] centre1 - Centroid of first triangle.
* @param[in] n2 - Normal vector of second triangle.
* @param[in] centre2 - Centroid of second triangle.
* @param[in] eps - Tolerance value.
* @return True if triangles share the same plane.
*/
bool coplanar(
    const Vec3& n1, const Vec3& centre1,
    const Vec3& n2, const Vec3& centre2,
    double eps);

/**
* @brief Computes the intersection line segment formed by triangle 1 cutting plane 2.
* @param[in] vertices - Array of Vec3 coordinates of first triangle's vertices.
* @param[in] distances - Signed distances of t1's vertices to plane 2.
* @param[in] eps - Tolerance value.
* @return Pair of 3D points defining the intersection line segment or null if no such pair exists.s
*/
std::optional<std::pair<Vec3, Vec3>> segmentAOnB(
    const TriVerts& vertices,
    const std::array<double, 3>& distances,
    const double eps);

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
std::vector<Vec3> findIntersectionPointsC(
    const TriVerts& vertices1,
    const TriVerts& vertices2, const Vec3& n2, const Vec3& c2,
    const double eps);

/**
* @brief Traverses two BVH trees to find all surface segment and face intersections between two meshes.
* @param[in] bvh1 - Bounding volume hierarchy of mesh 1.
* @param[in] bvh2 - Bounding volume hierarchy of mesh 2.
* @param[in] ctx - Full Collision Context struct containing both meshes and several output buffers for information.
*/
void findAllCollisions(const Bvh& bvh1, const Bvh& bvh2, CollisionContext& ctx);

/**
* @brief Constructs a Bounding Volume Hierarchy (BVH) tree for accelerated collision queries.
* @param[in] mesh - Target surface mesh.
* @return Built BVH tree.
 */
Bvh buildMeshBVH(const MeshData& mesh);

/**
* @brief High-level detection pipeline computing all geometric collisions between two meshes.
* @param[in] meshA - First input surface mesh.
* @param[in] meshB - Second input surface mesh.
* @return CollisionContext containing collision structures, cutting lines, and element indices.
*/
CollisionContext detectCollisions(const MeshData& meshA, const MeshData& meshB);