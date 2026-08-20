#pragma once
#include "vector.hpp"
#include "mesh_types.hpp"
#include "math_utils.hpp"
#include <optional>

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
* @brief Find all cycles in a PolyLine soup ("holes").
* @param[in] edges - Group of pairs of vectors (line segments) in arbitrary order.
* @param[in] nodeGrid - Spatial node matching structure.
* @param[out] outLoops - Output buffer containing all valid polyline loops found.
*/
void findCycles(const PolyLine& edges, SpatialGrid3D& nodeGrid, std::vector<PolyLine>& outLoops);

/**
* @brief Determines if a point projected into 2D space falls inside any cycle hole region.
* @param[in] centre - 3D Query coordinate.
* @param[in] allHoles - Vector of closed boundary polyline cycles defining hole regions.
* @param[in] frame - Planar projection frame.
* @return True if point lies inside any hole boundary.
*/
bool isCentroidInHole(
    const Vec3& centre,
    const std::vector<PolyLine>& allHoles, 
    const ProjectionFrame& frame);