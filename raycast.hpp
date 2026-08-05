#pragma once
#include "mesh_types.hpp"
#include "bvh_collisions.hpp"

/**
* @brief Performs ray-triangle intersection test.
* @param[in] orig - Ray origin coordinate.
* @param[in] dir - Ray direction vector.
* @param[in] t - Target triangle.
* @param[in] nodes - Global node list.
* @param[in] eps - Distance tolerance threshold.
* @return True if ray intersects triangle within positive distance.
*/
static bool intersectRayTri(
    const Vec3& orig, const Vec3& dir,
    const Triangle& t, const std::vector<Vec3>& nodes,
    double eps);

/**
* @brief Tests ray intersection with an axis-aligned bounding box (AABB).
* @param[in] orig - Ray origin position.
* @param[in] invDir - Precomputed component-wise inverse ray direction (1/dir).
* @param[in] box - Target bounding box.
* @return True if ray intersects bounding box.
*/
static bool rayBoxIntersect(const Vec3& orig, const Vec3& invDir, const BBox& box);

/**
* @brief Traverses BVH to count total ray intersections against a mesh surface.
* @param[in] bvh - Acceleration structure for target mesh.
* @param[in] mesh - Target surface mesh data.
* @param[in] orig - Ray origin point.
* @param[in] direction - Ray unit vector.
* @param[in] eps - Distance tolerance.
* @return Total number of valid ray-triangle intersections.
*/
int countRayIntersections(
    const Bvh& bvh, const MeshData& mesh,
    const Vec3& orig, const Vec3& direction,
    double eps);

/**
* @brief Determines if a 3D query point is located inside a closed surface mesh using parity ray casting.
* @param[in] bvh - Target mesh BVH tree.
* @param[in] mesh - Target surface mesh.
* @param[in] point - Query point in 3D space.
* @param[in] eps - Distance tolerance.
* @return True if point lies inside closed mesh volume.
*/
bool isInsideMesh(const Bvh& bvh, const MeshData& mesh,
    const Vec3& point,
    double eps);

/**
* @enum FaceClass
* @brief Topological classification of a mesh element relative to a target surface.
*/
enum class FaceClass { Outside, Inside, CoplanarSame, CoplanarOpp };

/**
* @brief Classifies a triangle face as inside, outside, or coplanar relative to a target surface mesh.
* @param[in] targetBvh - BVH tree of target closed domain mesh.
* @param[in] targetMesh - Target mesh dataset.
* @param[in] triCenter - Centroid of query triangle element.
* @param[in] triNormal - Normal vector of query triangle element.
* @param[in] eps - Distance tolerance.
* @return FaceClass classification result.
*/
FaceClass classifyFace(
    const Bvh& targetBvh, const MeshData& targetMesh,
    const Vec3& triCenter, const Vec3& triNormal,
    double eps);