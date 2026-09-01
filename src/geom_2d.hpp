#pragma once
#include "mesh_types.hpp"
#include "vector.hpp"
#include "math_utils.hpp"
using EdgeKey = std::pair<size_t, size_t>;
using Edge = std::pair<Vec2, Vec2>;
/**
* @brief Standard key of any edge, with smaller index first and larger index second.
* @param [in] u - Index of one endpoint of the edge.
* @param [in] v - Index of the other endpoint of the edge.
* @return The pair of indices in ascending order.
*/
inline EdgeKey makeEdgeKey(size_t u, size_t v) {
    return {std::min(u, v), std::max(u, v)};
}

/**
* @brief Constructs an orthonormal 2D projection frame for a planar surface.
* @param[in] normal - Plane normal vector.
* @param[in] origin - Reference origin point on plane.
* @return Orthonormal ProjectionFrame object.
*/
ProjectionFrame computeSharedFrame(const Vec3& normal, const Vec3& origin);

/**
* @brief Finds whether a point P is on the line segment AB
* @param[in] P - Point to test
* @param[in] A - Start point of segment.
* @param[in] B - End point of segment.
* @param[in] eps - Tolerance value.
* return True if point is on line segment, false otherwise
*/
bool pointOnSegment(const Vec2& P,const Vec2& A,const Vec2& B,
    double eps);

/**
* @brief Computes 2D line-segment intersection point.
* @param[in] A - Start point of segment 1.
* @param[in] B - End point of segment 1.
* @param[in] C - Start point of segment 2.
* @param[in] D - End point of segment 2.
* @param [in, out] grid - 2D grid used for vertex deduplication and storing.
* @param[in, out] outs - Vector of indices of points in the 2D spatial grid.
* @param[in] eps - Distance tolerance.
* @return True if segments intersect within valid bounds.
*/
void intersect2DAllPoints(
    const Vec2& A, const Vec2& B,
    const Vec2& C, const Vec2& D,
    SpatialGrid2D& grid,
    std::vector<size_t>& outs,
    double eps);

/**
* @brief Determines if a point is inside a polygon.
* @param[in] pt - 2D query coordinate.
* @param[in] edges - Vector of edges forming a closed loop.
* @return True if point lies inside the polygon.
*/
bool isPointInsidePolygon(const Vec2& pt, const std::vector<std::pair<Vec2, Vec2>>& edges);